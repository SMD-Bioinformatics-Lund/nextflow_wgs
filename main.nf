#!/usr/bin/env nextflow

// GENERAL PATHS //
OUTDIR = params.outdir+'/'+params.subdir
CRONDIR = params.crondir

// SENTIEON CONFIGS //
K_size      = 100000000
sentieon_model = params.sentieon_model
bwa_num_shards = params.bwa_shards
bwa_shards = Channel.from( 0..bwa_num_shards-1 )
genomic_num_shards = params.genomic_shards_num

// FASTA //
genome_file = params.genome_file

PON = [F: params.GATK_PON_FEMALE, M: params.GATK_PON_MALE]

// Count lines of input csv, if more than 2(header + 1 ind) then mode is set to family //
csv = file(params.csv)
mode = csv.countLines() > 2 ? "family" : "single"
trio = csv.countLines() > 3 ? true : false
println(csv)
println("mode: "+mode)
println("trio: "+trio)
// Print commit-version of active deployment
file(params.git)
	.readLines()
	.each { println "git commit-hash: "+it }
// Print active container
container = file(params.container).toRealPath()
println("container: "+container)

workflow.onComplete {

	def msg = """\
		Pipeline execution summary
		---------------------------
		Completed at: ${workflow.complete}
		Duration    : ${workflow.duration}
		Success     : ${workflow.success}
		scriptFile  : ${workflow.scriptFile}
		workDir     : ${workflow.workDir}
		csv         : ${params.csv}
		exit status : ${workflow.exitStatus}
		errorMessage: ${workflow.errorMessage}
		errorReport :
		"""
		.stripIndent()
	def error = """\
		${workflow.errorReport}
		"""
		.stripIndent()

	base = csv.getBaseName()
	File logFile = new File("${params.crondir}/logs/${base}.complete")
	if (!logFile.getParentFile().exists()) {
		logFile.getParentFile().mkdirs()
	}
	logFile.text = msg
	logFile.append(error)
}

// Input channels for alignment, variant calling and annotation //
Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, 
				row.id, 
				(row.containsKey("bam") ? file(row.bam) : (row.containsKey("vcf") ? file(row.vcf) : file(row.read1) ) ), 
				(row.containsKey("bai") ? file(row.bai) : (row.containsKey("idx") ? file(row.idx) : file(row.read2) ) ) ) }
	.set { input_files }

fastq = Channel.create()
bam_choice = Channel.create()
//vcf_choice = Channel.create()
fastq_sharded = Channel.create()
fastq_umi = Channel.create()
annotate_only = Channel.create()

// If input-files has bam files bypass alignment, otherwise go for fastq-channels => three options for fastq, sharded bwa, normal bwa or umi trimming
input_files.view().choice(bam_choice, fastq, fastq_sharded, fastq_umi, annotate_only ) { it[2] =~ /\.bam/ ? 0 : ( it[2] =~ /\.vcf.gz/ ? 4 : (params.shardbwa ? 2 : (params.umi ? 3 : 1) )) }

annotate_only.into{
	annotate_only_vep;
	annotate_only_cadd
}

Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, row.assay) }
        .set{ meta_loqusdb_no_sv_calling }

// Input channels for various meta information //
Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.id, row.diagnosis, row.read1, row.read2) }
	.set{ qc_extra }

Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, row.id, row.sex, row.mother, row.father, row.phenotype, row.diagnosis, row.type, row.assay, row.clarity_sample_id, (row.containsKey("ffpe") ? row.ffpe : false), (row.containsKey("analysis") ? row.analysis : false) ) }
	.into { ped; yml_diag; meta_upd; meta_str }


Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, row.id, row.sex, row.type) }
	.into { meta_gatkcov; meta_exp; meta_svbed; meta_pod; meta_mutect2; meta_eklipse}


Channel
	.fromPath(params.gatkreffolders)
	.splitCsv(header:true)
	.map{ row-> tuple(row.i, row.refpart) }
	.into{ gatk_ref; gatk_postprocess }


// Check whether genome assembly is indexed //
if(genome_file) {
	bwaId = Channel
			.fromPath("${genome_file}.bwt")
			.ifEmpty { exit 1, "BWA index not found: ${genome_file}.bwt" }
}

process fastp {
	cpus 10
	tag "$id"
	time '1h'
	memory '20 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.container_fastp}"

	when:
		params.umi

	input:
		set group, val(id), r1, r2 from fastq_umi

	output:
		set group, val(id), file("${id}_R1_a_q_u_trimmed.fq.gz"), file("${id}_R2_a_q_u_trimmed.fq.gz") into fastq_trimmed
		set group, file("*versions.yml") into ch_fastp_versions

	script:
		"""
		fastp -i $r1 -I $r2 --stdout \\
			-U --umi_loc=per_read --umi_len=3 \\
			-w ${task.cpus} \\
		| fastp --stdin --interleaved_in -f 2 -F 2 \\
			-o ${id}_R1_a_q_u_trimmed.fq.gz \\
			-O ${id}_R2_a_q_u_trimmed.fq.gz \\
			-l 30 \\
			-w ${task.cpus}
		
		${fastp_version(task)}
		"""

	stub:
		"""
		touch "${id}_R1_a_q_u_trimmed.fq.gz"
		touch "${id}_R2_a_q_u_trimmed.fq.gz"

		${fastp_version(task)}
		"""
}
def fastp_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    fastp: \$(echo \$(fastp -v 2>&1) | cut -f 2 -d " ")
	END_VERSIONS
	"""
}

// Align fractions of fastq files with BWA
process bwa_align_sharded {
	cpus 50
	memory '120 GB'
	tag "$id $shard"
	time '5h'
	container = "${params.container_sentieon}"

	input:
		set val(shard), val(group), val(id), r1, r2 from bwa_shards.combine(fastq_sharded)

	output:
		set val(id), group, file("${id}_${shard}.bwa.sort.bam"), file("${id}_${shard}.bwa.sort.bam.bai") into bwa_shards_ch
		set group, file("*versions.yml") into ch_bwa_align_shareded_versions

	when:
		params.align && params.shardbwa

	script:
		"""
		sentieon bwa mem -M \\
			-R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' \\
			-K $K_size \\
			-t ${task.cpus} \\
			-p $genome_file '<sentieon fqidx extract -F $shard/$bwa_num_shards -K $K_size $r1 $r2' | sentieon util sort \\
			-r $genome_file \\
			-o ${id}_${shard}.bwa.sort.bam \\
			-t ${task.cpus} --sam2bam -i -
		
		${bwa_align_sharded_version(task)}
		"""

	stub:
		"""
		touch "${id}_${shard}.bwa.sort.bam"
		touch "${id}_${shard}.bwa.sort.bam.bai"

		${bwa_align_sharded_version(task)}
		"""
}
def bwa_align_sharded_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon util --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	    bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
	END_VERSIONS
	"""
}

// Merge the fractioned bam files
process bwa_merge_shards {
	cpus 50
	tag "$id"
	time '1h'
	memory '120 GB'
	container = "${params.container_sentieon}"

	input:
		set val(id), group, file(shard), file(shard_bai) from bwa_shards_ch.groupTuple(by: [0,1])

	output:
		set id, group, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into merged_bam_locusc
		set id, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into merged_bam_dedup
		set group, file("*versions.yml") into ch_bwa_merge_shards_versions
	when:
		params.shardbwa
	
	script:
		bams = shard.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() } .join(' ')

		"""
		sentieon util merge -o ${id}_merged.bam ${bams}

		${bwa_merge_shards_version(task)}
		"""

	stub:
		"""
		touch "${id}_merged.bam"
		touch "${id}_merged.bam.bai"

		${bwa_merge_shards_version(task)}
		"""
}
def bwa_merge_shards_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon util --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}

// ALTERNATIVE PATH: Unsharded BWA, utilize local scratch space.
process bwa_align {
	cpus 50
	memory '100 GB'
	// 64 GB peak giab //
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	tag "$id"
	container = "${params.container_sentieon}"

	input:
		set val(group), val(id), file(r1), file(r2) from fastq.mix(fastq_trimmed)

	output:
		set id, group, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into bam_locusc, bam_markdup
		// set id, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into bam_dedup remnant of distri
		set group, file("*versions.yml") into ch_bwa_align_versions

	when:
		params.align && !params.shardbwa

	script:
		"""
		sentieon bwa mem \\
			-M \\
			-K $K_size \\
			-R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' \\
			-t ${task.cpus} \\
			$genome_file $r1 $r2 \\
			| sentieon util sort \\
			-r $genome_file \\
			-o ${id}_merged.bam \\
			-t ${task.cpus} --sam2bam -i -
		
		${bwa_align_versions(task)}
		"""

	stub:
		"""
		touch "${id}_merged.bam"
		touch "${id}_merged.bam.bai"

		${bwa_align_versions(task)}
		"""

}
def bwa_align_versions(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon util --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	    bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
	END_VERSIONS
	"""
}

process markdup {
	cpus 40
	errorStrategy 'retry'
	maxErrors 5
	tag "$id"
	memory '50 GB'
	// 12gb peak giab //
	time '3h'
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'
	container = "${params.container_sentieon}"
	publishDir "${OUTDIR}/bam", mode: 'copy' , overwrite: 'true', pattern: '*_dedup.bam*'

	input:
		set id, group, file(bam), file(bai) from bam_markdup.mix(merged_bam_dedup)

	output:
		set group, id, file("${id}_dedup.bam"), file("${id}_dedup.bam.bai") into complete_bam, chanjo_bam, d4_bam, verifybamid2_bam, expansionhunter_bam, yaml_bam, cov_bam, bam_manta, bam_nator, bam_tiddit, bam_manta_panel, bam_cnvkit_panel, bam_freebayes, bam_mito, smncnc_bam, bam_gatk, depth_onco
		set id, group, file("${id}_dedup.bam"), file("${id}_dedup.bam.bai") into qc_bam, bam_melt, bam_bqsr
		set val(id), file("dedup_metrics.txt") into dedupmet_sentieonqc
		set group, file("${group}_bam.INFO") into bam_INFO
		set group, file("*versions.yml") into ch_markdup_versions

	script:
		"""
		sentieon driver \\
			--temp_dir /local/scratch/ \\
			-t ${task.cpus} \\
			-i $bam --shard 1:1-248956422 --shard 2:1-242193529 --shard 3:1-198295559 --shard 4:1-190214555 --shard 5:1-120339935 --shard 5:120339936-181538259 --shard 6:1-170805979 --shard 7:1-159345973 --shard 8:1-145138636 --shard 9:1-138394717 --shard 10:1-133797422 --shard 11:1-135086622 --shard 12:1-56232327 --shard 12:56232328-133275309 --shard 13:1-114364328 --shard 14:1-107043718 --shard 15:1-101991189 --shard 16:1-90338345 --shard 17:1-83257441 --shard 18:1-80373285 --shard 19:1-58617616 --shard 20:1-64444167 --shard 21:1-46709983 --shard 22:1-50818468 --shard X:1-124998478 --shard X:124998479-156040895 --shard Y:1-57227415 --shard M:1-16569 \\
			--algo LocusCollector \\
			--fun score_info ${id}.score
		sentieon driver \\
			--temp_dir /local/scratch/ \\
			-t ${task.cpus} \\
			-i $bam \\
			--algo Dedup --score_info ${id}.score \\
			--metrics dedup_metrics.txt \\
			--rmdup ${id}_dedup.bam
		echo "BAM	$id	/access/${params.subdir}/bam/${id}_dedup.bam" > ${group}_bam.INFO

		${markdup_versions(task)}
		"""

	stub:
		"""
		touch "${id}_dedup.bam"
		touch "${id}_dedup.bam.bai"
		touch "dedup_metrics.txt"
		touch "${group}_bam.INFO"

		${markdup_versions(task)}
		"""
}
def markdup_versions(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}

process copy_bam {

	tag "$id"
	cpus 1
	memory '2GB'
	time '1h'

	input:
		set group, id, file(bam), file(bai) from bam_choice
	
	output:
		set group, id, file("${id}_dedup.bam"), file("${id}_dedup.bam.bai") into expansionhunter_bam_choice, dnascope_bam_choice, bampath_start, cov_bam_choice, bam_manta_choice, bam_tiddit_choice, bam_mito_choice, bam_SMN_choice, bam_freebayes_choice, bam_mantapanel_choice, bam_cnvkitpanel_choice, bam_melt_choice, bam_qc_choice, dedup_dummy_choice, bam_bqsr_choice, bam_gatk_choice, verifybamid2_bam_choice
	script:
		"""
		ionice -c 2 -n 7 cp ${bam} "${id}_dedup.copy.bam"
		ionice -c 2 -n 7 cp ${bai} "${id}_dedup.copy.bam.bai"
		"""

	stub:
		"""
		touch "${id}_dedup.copy.bam"
		touch "${id}_dedup.copy.bam.bai"
		"""
}

// These processes expects ID, group instead of group, ID
// FIXME: We should really fix this to be part of the channels
def remap_bam_choice_tuple = { channel ->
		channel.map { tup -> return tuple(tup.get(1), tup.get(0), tup.get(2), tup.get(3))
	}
}
remap_bam_choice_tuple(cov_bam_choice).set { cov_bam_choice }
remap_bam_choice_tuple(bam_melt_choice).set { bam_melt_choice }
remap_bam_choice_tuple(bam_bqsr_choice).set { bam_bqsr_choice }
remap_bam_choice_tuple(bam_qc_choice).set { bam_qc_choice }

// For melt to work if started from bam-file.
process dedupdummy {
	when:
		params.run_melt
	input:
		set group, id, file(bam), file(bai) from dedup_dummy_choice
	output:
		set id, file("dummy") into dedup_dummy
	"""
	echo test > dummy
	"""
}


process bqsr {
	cpus 40
	errorStrategy 'retry'
	maxErrors 5
	tag "$id"
	memory '30 GB'
	// 12gb peak giab //
	time '5h'
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'
	container = "${params.container_sentieon}"
	publishDir "${OUTDIR}/bqsr", mode: 'copy' , overwrite: 'true', pattern: '*.table'

	input:
		set id, group, file(bam), file(bai) from bam_bqsr.mix(bam_bqsr_choice)

	output:
		set group, id, file("${id}.bqsr.table") into dnascope_bqsr
		set group, file("*versions.yml") into ch_bqsr_versions

	script:
		"""
		sentieon driver -t ${task.cpus} \\
			-r $genome_file -i $bam \\
			--algo QualCal ${id}.bqsr.table \\
			-k $params.KNOWN

		${bqsr_version(task)}
		"""

	stub:
		"""
		touch "${id}.bqsr.table"
		${bqsr_version(task)}
		"""
}
def bqsr_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}

//Collect various QC data: 
process sentieon_qc {
	cpus 52
	memory '30 GB'
	tag "$id"
	time '2h'
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'
	container = "${params.container_sentieon}"

	input:
		set id, group, file(bam), file(bai) from qc_bam.mix(bam_qc_choice)

	output:
		set id, group, file("mq_metrics.txt"), file("qd_metrics.txt"), file("gc_summary.txt"), 
			file("gc_metrics.txt"), file("aln_metrics.txt"), file("is_metrics.txt"), file("assay_metrics.txt"), 
			file("cov_metrics.txt"), file("cov_metrics.txt.sample_summary") into ch_sentieon_qc_metrics
		set group, file("*versions.yml") into ch_sentieon_qc_versions

	script:
		target = ""
		// A bit of cheating here - these are really optional arguments
		panel_command = "touch cov_metrics.txt cov_metrics.txt.sample_summary"
		cov = "WgsMetricsAlgo assay_metrics.txt"

		if (params.onco || params.exome) {
			target = "--interval $params.intervals"
			cov = "CoverageMetrics --cov_thresh 1 --cov_thresh 10 --cov_thresh 30 --cov_thresh 100 --cov_thresh 250 --cov_thresh 500 cov_metrics.txt"
			panel_command = "sentieon driver -r ${params.genome_file} -t ${task.cpus} -i ${bam} --algo HsMetricAlgo --targets_list ${params.intervals} --baits_list ${params.intervals} assay_metrics.txt"
		}
		
		"""
		sentieon driver \\
			-r $genome_file $target \\
			-t ${task.cpus} \\
			-i $bam \\
			--algo MeanQualityByCycle mq_metrics.txt \\
			--algo QualDistribution qd_metrics.txt \\
			--algo GCBias --summary gc_summary.txt gc_metrics.txt \\
			--algo AlignmentStat aln_metrics.txt \\
			--algo InsertSizeMetricAlgo is_metrics.txt \\
			--algo $cov
		$panel_command

		${sentieon_qc_version(task)}
		"""

	stub:
		"""
		touch "assay_metrics.txt"
		touch "mq_metrics.txt"
		touch "qd_metrics.txt"
		touch "gc_summary.txt"
		touch "gc_metrics.txt"
		touch "aln_metrics.txt"
		touch "is_metrics.txt"
		touch "cov_metrics.txt"
		touch "cov_metrics.txt.sample_summary"
		${sentieon_qc_version(task)}
		"""
}
def sentieon_qc_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}

process sentieon_qc_postprocess {
	cpus 2
	memory '1 GB'
	tag "$id"
	time '2h'
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'

	input:
		set id, file(dedup) from dedupmet_sentieonqc.mix(dedup_dummy)
		set id, group, file(mq_metrics), file(qd_metrics), file(gc_summary), file(gc_metrics), file(aln_metrics),
			file(is_metrics), file(assay_metrics), file(cov_metrics), file(cov_metrics_sample_summary) from ch_sentieon_qc_metrics

	output:
		set group, id, file("${id}_qc.json") into qc_cdm
		set group, id, file("${id}_qc.json") into qc_melt
	
	script:

		assay = (params.onco || params.exome) ? "panel" : "wgs"
		"""
		qc_sentieon.pl \\
			--SID ${id} \\
			--type ${assay} \\
			--align_metrics_file ${aln_metrics} \\
			--insert_file ${is_metrics} \\
			--dedup_metrics_file ${dedup} \\
			--metrics_file ${assay_metrics} \\
			--gcsummary_file ${gc_summary} \\
			--coverage_file ${cov_metrics} \\
			--coverage_file_summary ${cov_metrics_sample_summary} \\
			> ${id}_qc.json
		
		"""
}

process d4_coverage {
	cpus 16
	memory '10 GB'
	publishDir "${OUTDIR}/cov", mode: 'copy', overwrite: 'true', pattern: '*.d4'
	tag "$id"
	container = "${params.container_d4tools}"

	when:
		params.run_chanjo2

	input:
		set group, id, file(bam), file(bai) from d4_bam

	output:
		file("${id}_coverage.d4")
		set group, id, file("${id}_coverage.d4") into ch_final_d4
		set group, file("*versions.yml") into ch_d4_coverage_versions
		set group, file("${group}_d4.INFO") into d4_INFO

	script:
	"""
	d4tools create \\
		--threads ${task.cpus} \\
		"${bam}" \\
		"${id}_coverage.d4"

	echo "D4	$id	/access/${params.subdir}/cov/${id}_coverage.d4" > ${group}_d4.INFO

	${d4_coverage_version(task)}
	"""

	stub:
	"""
	touch "${id}_coverage.d4"
	touch "${group}_d4.INFO"

	${d4_coverage_version(task)}
	"""
}
def d4_coverage_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    d4tools: \$(echo \$( d4tools 2>&1 | head -1 ) | sed "s/.*version: //" | sed "s/)//" )
	END_VERSIONS
	"""
}

process verifybamid2 {
	cpus 16
	memory '10 GB'
	// publishDir "${OUTDIR}/contamination", mode: 'copy', overwrite: 'true', pattern: '*.selfSM'
	tag "$id"
	container = "${params.container_verifybamid2}"

	input:
		set group, id, file(bam), file(bai) from verifybamid2_bam.mix(verifybamid2_bam_choice)

	output:
		file("${id}.result.selfSM")
		file("${id}.result.Ancestry")
		set group, file("*versions.yml") into ch_verifybamid2_versions

	script:
		if ( params.antype == "wgs") {
			"""
			verifybamid2 \
				--SVDPrefix ${params.verifybamid2_svdprefix} \
				--Reference ${genome_file} \
				--BamFile ${bam}

				mv result.selfSM ${id}.result.selfSM
				mv result.Ancestry ${id}.result.Ancestry
			${verifybamid2_version(task)}
			"""
		}
		else {
			"""
			verifybamid2 \
				--DisableSanityCheck \
				--SVDPrefix ${params.verifybamid2_svdprefix} \
				--Reference ${genome_file} \
				--BamFile ${bam}

				mv result.selfSM ${id}.result.selfSM
				mv result.Ancestry ${id}.result.Ancestry
			${verifybamid2_version(task)}
			"""
		}


	stub:
		"""
		touch "${id}.result.selfSM"
		touch "${id}.result.Ancestry"

		${verifybamid2_version(task)}
		"""
}
def verifybamid2_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    VerifyBamID2: \$( echo \$( verifybamid2 --help 2>&1 | grep Version ) | sed "s/^.*Version://" )
	END_VERSIONS
	"""
}

// Calculate coverage for paneldepth
process depth_onco {
	cpus 2
	time '1h'
	memory '10 GB'
	publishDir "${OUTDIR}/cov", mode: 'copy', overwrite: 'true'
	tag "$id"
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'

	when:
		params.assay == "swea"

	input:	
		set group, id, file(bam), file(bai) from depth_onco

	output:
		file("${id}.lowcov.overlapping.bed") into cov_onco

	script:
		"""
		panel_depth.pl $bam $params.scoutbed > ${id}.lowcov.bed
		overlapping_genes.pl ${id}.lowcov.bed $params.gene_regions > ${id}.lowcov.overlapping.bed
		"""

	stub:
		"""
		touch "${id}.lowcov.overlapping.bed"
		"""
}

process SMNCopyNumberCaller {
	cpus 10
	memory '25GB'
	time '2h'
	publishDir "${OUTDIR}/plots/SMNcnc", mode: 'copy' , overwrite: 'true', pattern: '*.pdf*'
	tag "$id"

	when:
		params.antype == "wgs"

	input:
		set group, id, file(bam), file(bai) from smncnc_bam.mix(bam_SMN_choice)

	output:
		file("*.tsv") into smn_tsv
		set file("*.pdf"), file("*.json")
		set group, file("${group}_smn.INFO") into smn_INFO
		set group, file("*versions.yml") into ch_smn_copy_number_caller_versions

	script:
		"""
		samtools view -H $bam | \\
			sed -e 's/SN:1/SN:chr1/' | sed -e 's/SN:2/SN:chr2/' |  \\
			sed -e 's/SN:3/SN:chr3/' | sed -e 's/SN:4/SN:chr4/' |  \\
			sed -e 's/SN:5/SN:chr5/' | sed -e 's/SN:6/SN:chr6/' |  \\
			sed -e 's/SN:7/SN:chr7/' | sed -e 's/SN:8/SN:chr8/' |  \\
			sed -e 's/SN:9/SN:chr9/' | sed -e 's/SN:10/SN:chr10/' | \\
			sed -e 's/SN:11/SN:chr11/' | sed -e 's/SN:12/SN:chr12/' |  \\
			sed -e 's/SN:13/SN:chr13/' | sed -e 's/SN:14/SN:chr14/' |  \\
			sed -e 's/SN:15/SN:chr15/' | sed -e 's/SN:16/SN:chr16/' |  \\
			sed -e 's/SN:17/SN:chr17/' | sed -e 's/SN:18/SN:chr18/' |  \\
			sed -e 's/SN:19/SN:chr19/' | sed -e 's/SN:20/SN:chr20/' |  \\
			sed -e 's/SN:21/SN:chr21/' | sed -e 's/SN:22/SN:chr22/' |  \\
			sed -e 's/SN:X/SN:chrX/' | sed -e 's/SN:Y/SN:chrY/' |   \\
			sed -e 's/SN:MT/SN:chrM/' | \\
			samtools reheader - $bam > ${id}.bam
		samtools index -b ${id}.bam -@ ${task.cpus}
		echo ${id}.bam > manifest.txt
		smn_caller.py --manifest manifest.txt --genome 38 --prefix ${id} --outDir . --threads ${task.cpus}
		rm ${id}.bam
		source activate py3-env
		python /SMNCopyNumberCaller/smn_charts.py -s ${id}.json -o .
		mv ${id}.tsv ${group}_SMN.tsv
		echo "SMN ${params.accessdir}/smn/${group}_SMN.tsv" > ${group}_smn.INFO

		${smn_copy_number_caller_version(task)}
		"""

	stub:
		"""
		touch "${id}.bam"
		touch "${id}.tsv"
		touch "${id}.pdf"
		touch "${id}.json"
		touch "${group}_SMN.tsv"
		touch "${group}_smn.INFO"
	
		${smn_copy_number_caller_version(task)}
		"""
}
// collects each individual's SMNCNC-tsv and creates one tsv-file
smn_tsv
	.collectFile(keepHeader: true, storeDir: "${OUTDIR}/smn/")
def smn_copy_number_caller_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
	    smn-copy-number-caller: 1.1.2
	END_VERSIONS	
	"""
}


////////////////////////////////////////////////////////////////////////
////////////////////////// EXPANSION HUNTER ////////////////////////////
////////////////////////////////////////////////////////////////////////

// call STRs using ExpansionHunter, and plot alignments with GraphAlignmentViewer
process expansionhunter {
	tag "$group"
	cpus 2
	time '10h'
	memory '40 GB'
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'

	when:
		params.str
		
	input:
		set group, id, file(bam), file(bai), sex, type \
			from expansionhunter_bam.mix(expansionhunter_bam_choice).join(meta_exp, by: [0,1]).filter { item -> item[5] == 'proband' }

	output:
		set group, id, file("${group}.eh.vcf") into expansionhunter_vcf
		set group, id, file("${group}.eh_realigned.sort.bam"), file("${group}.eh_realigned.sort.bam.bai"), file("${group}.eh.vcf") into reviewer
		set group, file("*versions.yml") into ch_expansionhunter_versions

	script:
		"""
		source activate htslib10
		ExpansionHunter \
			--reads $bam \
			--reference $genome_file \
			--variant-catalog $params.expansionhunter_catalog \
			--output-prefix ${group}.eh
		samtools sort ${group}.eh_realigned.bam -o ${group}.eh_realigned.sort.bam
		samtools index ${group}.eh_realigned.sort.bam

		${expansionhunter_version(task)}
		"""

	stub:
		"""
		source activate htslib10
		touch "${group}.eh.vcf"
		touch "${group}.eh_realigned.sort.bam"
		touch "${group}.eh_realigned.sort.bam.bai"

		${expansionhunter_version(task)}
		"""
}
def expansionhunter_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    expansionhunter: \$(echo \$(ExpansionHunter --version 2>&1) | sed 's/.*ExpansionHunter v// ; s/]//')
	    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
	END_VERSIONS
	"""
}

// annotate expansionhunter vcf
process stranger {
	tag "$group"
	memory '1 GB'
	time '10m'
	cpus 2
	container = "${params.container_stranger}"

	input:
		set group, id, file(eh_vcf) from expansionhunter_vcf

	output:
		set group, id, file("${group}.fixinfo.eh.stranger.vcf") into expansionhunter_vcf_anno
		set group, file("*versions.yml") into ch_stranger_versions

	script:
		"""
		stranger ${eh_vcf} -f $params.expansionhunter_catalog > ${group}.eh.stranger.vcf
		grep ^# ${group}.eh.stranger.vcf > ${group}.fixinfo.eh.stranger.vcf
		grep -v ^# ${group}.eh.stranger.vcf | sed 's/ /_/g' >> ${group}.fixinfo.eh.stranger.vcf
		
		${stranger_version(task)}
		"""

	stub:
		"""
		touch "${group}.fixinfo.eh.stranger.vcf"
		${stranger_version(task)}
		"""
}
def stranger_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    stranger: \$( stranger --version )
	END_VERSIONS
	"""
}

//for i in $( ls *.svg | cut -f 2 -d "." ); do echo "STR_IMG $i /access/!{params.subdir}/plots/reviewer/!{group}/!{group}.${i}.svg" >> !{group}_rev.INFO; done
process reviewer {
	tag "$group"
	cpus 2
	time '1h'
	memory '1 GB'
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'
	errorStrategy 'ignore'
	container = "${params.container_reviewer}"
	publishDir "${OUTDIR}/plots/reviewer/${group}", mode: 'copy' , overwrite: 'true', pattern: '*.svg'
	
	input:
		set group, id, file(bam), file(bai), file(vcf) from reviewer

	output:
		file("*svg")
		//set group, file("${group}_rev.INFO") into reviewer_INFO
		set group, file("*versions.yml") into ch_reviewer_versions

	shell:
		version_str = reviewer_version(task)
		'''
		grep LocusId !{params.expansionhunter_catalog} | sed 's/[",^ ]//g' | cut -d':' -f2 | perl -na -e 'chomp; \
		system("REViewer --reads !{bam} \
		--vcf !{vcf} \
		--reference !{genome_file} \
		--catalog !{params.expansionhunter_catalog} \
		--locus $_ \
		--output-prefix !{id}");'

		echo "!{version_str}" > "!{task.process}_versions.yml"
		'''

	stub:
		version_str = reviewer_version(task)
		"""
		touch "${id}.svg"
		echo "${version_str}" > "${task.process}_versions.yml"
		"""
}
def reviewer_version(task) {
	// This docstring looks different from others as it is used within the shell '''
	// If spaces similarly to the others, this leads to additional whitespace above and below the version text
	"""${task.process}:
	    reviewer: \$(echo \$(REViewer --version 2>&1) | sed 's/^.*REViewer v//')"""
}

// split multiallelic sites in expansionhunter vcf
// FIXME: Use env variable for picard path...
process vcfbreakmulti_expansionhunter {
	publishDir "${OUTDIR}/vcf", mode: 'copy' , overwrite: 'true', pattern: '*.vcf.gz'
	tag "$group"
	time '1h'
	memory '50 GB'
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'

	input:
		set group, id, file(eh_vcf_anno) from expansionhunter_vcf_anno
		set group, id, sex, mother, father, phenotype, diagnosis, type, assay, clarity_sample_id, ffpe, analysis from meta_str.filter{ item -> item[7] == 'proband' }

	output:
		file("${group}.expansionhunter.vcf.gz") into expansionhunter_scout
		set group, file("${group}_str.INFO") into str_INFO
		set group, file("*versions.yml") into ch_vcfbreakmulti_expansionhunter_versions

	script:
		if (father == "") { father = "null" }
		if (mother == "") { mother = "null" }
		if (mode == "family") {
			"""
			java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf INPUT=${eh_vcf_anno} OUTPUT=${eh_vcf_anno}.rename.vcf NEW_SAMPLE_NAME=${id}
			vcfbreakmulti ${eh_vcf_anno}.rename.vcf > ${group}.expansionhunter.vcf.tmp
			familyfy_str.pl --vcf ${group}.expansionhunter.vcf.tmp --mother $mother --father $father --out ${group}.expansionhunter.vcf
			bgzip ${group}.expansionhunter.vcf
			tabix ${group}.expansionhunter.vcf.gz
			echo "STR	${params.accessdir}/vcf/${group}.expansionhunter.vcf.gz" > ${group}_str.INFO

			${vcfbreakmulti_expansionhunter_version(task)}
			"""
		}
		else {
			"""
			java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf INPUT=${eh_vcf_anno} OUTPUT=${eh_vcf_anno}.rename.vcf NEW_SAMPLE_NAME=${id}
			vcfbreakmulti ${eh_vcf_anno}.rename.vcf > ${group}.expansionhunter.vcf
			bgzip ${group}.expansionhunter.vcf
			tabix ${group}.expansionhunter.vcf.gz
			echo "STR	${params.accessdir}/vcf/${group}.expansionhunter.vcf.gz" > ${group}_str.INFO

			${vcfbreakmulti_expansionhunter_version(task)}
			"""

		}

	stub:
		"""
		touch "${group}.expansionhunter.vcf.gz"
		touch "${group}_str.INFO"

		${vcfbreakmulti_expansionhunter_version(task)}
		"""
}
def vcfbreakmulti_expansionhunter_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcflib: 1.0.9
	    rename-sample-in-vcf: \$(echo \$(java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf --version 2>&1) | sed 's/-SNAPSHOT//')
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	END_VERSIONS	
	"""
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
process melt_qc_val {
	tag "$id"
	time '20m'
	memory '50 MB'

	when:
		params.run_melt

	input:
		set group, id, qc from qc_melt

	output:
		set id, val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) into qc_melt_val
		set group, id, val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) into qc_cnvkit_val
	
	script:
		// Collect qc-data if possible from normal sample, if only tumor; tumor
		qc.readLines().each{
			if (it =~ /\"(ins_size_dev)\" : \"(\S+)\"/) {
				ins_dev = it =~ /\"(ins_size_dev)\" : \"(\S+)\"/
			}
			if (it =~ /\"(mean_coverage)\" : \"(\S+)\"/) {
				coverage = it =~ /\"(mean_coverage)\" : \"(\S+)\"/
			}
			if (it =~ /\"(ins_size)\" : \"(\S+)\"/) {
				ins_size = it =~ /\"(ins_size)\" : \"(\S+)\"/
			}
		}
		// might need to be defined for -resume to work "def INS_SIZE" and so on....
		INS_SIZE = ins_size[0][2]
		MEAN_DEPTH = coverage[0][2]
		COV_DEV = ins_dev[0][2]
		"""
		echo hej > hej
		"""

	stub:
		INS_SIZE = 0
		MEAN_DEPTH = 0
		COV_DEV = 0
		"""
		echo $INS_SIZE $MEAN_DEPTH $COV_DEV > qc.val
		"""
}

// MELT always give VCFs for each type of element defined in mei_list
// If none found -> 0 byte vcf. merge_melt.pl merges the three, if all empty
// it creates a vcf with only header from params.meltheader
// merge_melt.pl gives output ${id}.melt.merged.vcf
process melt {
	cpus 3
	errorStrategy 'retry'
	container = "${params.container_melt}"
	tag "$id"
	// memory seems to scale with less number of reads?
	memory '70 GB'
	time '3h'
	publishDir "${OUTDIR}/vcf", mode: 'copy' , overwrite: 'true', pattern: '*.vcf'

	input:
		set id, group, file(bam), file(bai), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) from bam_melt.mix(bam_melt_choice).join(qc_melt_val)

	when:
		params.run_melt

	output:
		set group, id, file("${id}.melt.merged.vcf") into melt_vcf_nonfiltered
		set group, file("*versions.yml") into ch_melt_versions

	script:
		"""
		java -jar /opt/MELTv2.2.2/MELT.jar Single \\
			-bamfile $bam \\
			-r 150 \\
			-h $genome_file \\
			-n /opt/MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed \\
			-z 500000 \\
			-d 50 \\
			-t $params.mei_list \\
			-w . \\
			-c $MEAN_DEPTH \\
			-e $INS_SIZE \\
			-exome
		merge_melt.pl $params.meltheader $id

		${melt_version(task)}
		"""

	stub:
		"""
		touch "${id}.melt.merged.vcf"
		${melt_version(task)}
		"""
}
def melt_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    melt: \$(echo \$(java -jar /opt/MELTv2.2.2/MELT.jar -h | grep "^MELTv" | cut -f1 -d" " | sed "s/MELTv//" ) )
	END_VERSIONS	
	"""
}

process intersect_melt {
	cpus 2
	tag "$id"
	memory '2 GB'
	time '1h'
	publishDir "${OUTDIR}/vcf", mode: 'copy' , overwrite: 'true', pattern: '*.vcf'

	input:
		set group, id, file(vcf) from melt_vcf_nonfiltered

	when:
		params.run_melt

	output:
		set group, id, file("${id}.melt.merged.intersected.vcf") into ch_melt_vcf
		set group, file("*versions.yml") into ch_intersect_melt_versions

	script:
		"""
		bedtools intersect -a $vcf -b $params.intersect_bed -header > ${id}.melt.merged.intersected.vcf
		${intersect_melt_version(task)}
		"""

	stub:
		"""
		touch "${id}.melt.merged.intersected.vcf"
		${intersect_melt_version(task)}
		"""
}
def intersect_melt_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bedtools: \$(echo \$(bedtools --version 2>&1) | sed -e "s/^.*bedtools v//" )
	END_VERSIONS
	"""
}

// When rerunning sample from bam, dnascope has to be run unsharded. this is mixed together with all other vcfs in a trio //
process dnascope {
	cpus 54
	memory '100 GB'
	// 12 GB peak giab //
	time '4h'
	tag "$id"
	container = "${params.container_sentieon}"

	when:
		params.varcall

	input:
		set group, id, bam, bai, bqsr from complete_bam.mix(dnascope_bam_choice).join(dnascope_bqsr, by: [0,1])

	output:
		set group, id, file("${id}.dnascope.gvcf.gz"), file("${id}.dnascope.gvcf.gz.tbi") into complete_vcf_choice
		set group, id, file("${id}.dnascope.gvcf.gz") into gvcf_gens_choice
		set group, file("*versions.yml") into ch_dnascope_versions

	script:
		"""
		sentieon driver \\
			-t ${task.cpus} \\
			-r $genome_file \\
			-q $bqsr \\
			-i $bam --shard 1:1-248956422 --shard 2:1-242193529 --shard 3:1-198295559 --shard 4:1-190214555 --shard 5:1-120339935 --shard 5:120339936-181538259 --shard 6:1-170805979 --shard 7:1-159345973 --shard 8:1-145138636 --shard 9:1-138394717 --shard 10:1-133797422 --shard 11:1-135086622 --shard 12:1-56232327 --shard 12:56232328-133275309 --shard 13:1-114364328 --shard 14:1-107043718 --shard 15:1-101991189 --shard 16:1-90338345 --shard 17:1-83257441 --shard 18:1-80373285 --shard 19:1-58617616 --shard 20:1-64444167 --shard 21:1-46709983 --shard 22:1-50818468 --shard X:1-124998478 --shard X:124998479-156040895 --shard Y:1-57227415 --shard M:1-16569 \\
			--algo DNAscope --emit_mode GVCF ${id}.dnascope.gvcf.gz

		${dnascope_version(task)}
		"""

	stub:
		"""
		touch "${id}.dnascope.gvcf.gz"
		touch "${id}.dnascope.gvcf.gz.tbi"

		${dnascope_version(task)}
		"""
}
def dnascope_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS	
	"""
}

process bamtoyaml {
	cpus 1
	time "5m"
	memory "2MB"

	input:
		set group, id, bam, bai from bampath_start
	
	output:
		set group, file("${group}_bamstart.INFO") into bamchoice_INFO

	script:
		"""
		echo "BAM	$id	/access/${params.subdir}/bam/${bam.getName()}" > ${group}_bamstart.INFO
		"""

	stub:
		"""
		touch "${group}_bamstart.INFO"
		"""
}


process gvcf_combine {
	cpus 16
	tag "$group"
	memory '5 GB'
	time '5h'
	container = "${params.container_sentieon}"

	input:
		set group, id, file(vcf), file(idx) from complete_vcf_choice.groupTuple()

	output:
		set group, id, file("${group}.combined.vcf"), file("${group}.combined.vcf.idx") into combined_vcf
		set group, file("*versions.yml") into ch_gvcf_combine_versions

	script:
		all_gvcfs = vcf.collect { it.toString() }.sort().join(' -v ')

		"""
		sentieon driver \\
			-t ${task.cpus} \\
			-r $genome_file \\
			--algo GVCFtyper \\
			-v $all_gvcfs ${group}.combined.vcf

		${gvcf_combine_version(task)}
		"""

	stub:
		"""
		touch "${group}.combined.vcf"
		touch "${group}.combined.vcf.idx"

		${gvcf_combine_version(task)}
		"""
}
def gvcf_combine_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}

// Create ped from input variables //
process create_ped {
	tag "$group"
	time '20m'
	publishDir "${OUTDIR}/ped", mode: 'copy' , overwrite: 'true'
	memory '1 GB'

	input:
		set group, id, sex, mother, father, phenotype, diagnosis, type, assay, clarity_sample_id, ffpe, analysis from ped.filter { item -> item[7] == 'proband' }
		
	output:
		set group, type, file("${group}_base.ped") into ped_mad, ped_peddy, ped_inher, ped_scout, ped_loqus, ped_prescore, ped_compound, ped_pod
		set group, type_ma, file("${group}_ma.ped") optional true into ped_inher_ma, ped_prescore_ma, ped_compound_ma, ped_mad_ma
		set group, type_fa, file("${group}_fa.ped") optional true into ped_inher_fa, ped_prescore_fa, ped_compound_fa, ped_mad_fa

	script:
		if ( father == "" ) {
			father = "0"
		}
		if ( mother == "" ) {
			mother = "0"
		}
		type_fa = "fa"
		type_ma = "ma"
		"""
		create_ped.pl --mother $mother --father $father --group $group --id $id --sex $sex
		"""

	stub:
		type_fa = "fa"
		type_ma = "ma"
		"""
		touch "${group}_base.ped"
		touch "${group}_ma.ped"
		touch "${group}_fa.ped"

        echo $type_fa $type_ma > type.val
		"""
}

//madeline ped, run if family mode
process madeline {
	publishDir "${OUTDIR}/ped", mode: 'copy' , overwrite: 'true', pattern: '*.xml'
	memory '1 GB'
	time '1h'
	cpus 2
	container = "${params.container_madeline}"

	input:
		set group, type, file(ped) from ped_mad.mix(ped_mad_ma,ped_mad_fa)

	output:
		file("${ped}.madeline.xml") into madeline_ped
		set group, file("${group}_madde.INFO") into madde_INFO
		set group, file("*versions.yml") into ch_madeline_versions

	when:
		mode == "family" && params.assay == "wgs"

	script:
		"""
		source activate tools
		ped_parser \\
			-t ped $ped \\
			--to_madeline \\
			-o ${ped}.madeline
		madeline2 \\
			-L "IndividualId" ${ped}.madeline \\
			-o ${ped}.madeline \\
			-x xml
		echo "MADDE	$type ${params.accessdir}/ped/${ped}.madeline.xml" > ${group}_madde.INFO

		${madeline_version(task)}
		"""

	stub:
		"""
		source activate tools
		touch "${group}_madde.INFO"
		touch "${ped}.madeline.xml"

		${madeline_version(task)}
		"""
}
def madeline_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    ped-parser: \$(echo \$(ped_parser --version 2>&1) | sed -e "s/^.*ped_parser version: //")
	    madeline: \$(echo \$(madeline2 --version 2>&1) | grep : | sed -e"s/^.*Madeline //; s/PDE : 1.*//")
	END_VERSIONS	
	"""
}

process freebayes {
	cpus 1
	time '2h'
	memory '10 GB'
	container = "${params.container_twist_myeloid}"
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'

	when: 
		params.antype == "panel"

	input:
		set group, id, file(bam), file(bai) from bam_freebayes.mix(bam_freebayes_choice)

	output:
		set group, file("${id}.pathfreebayes.lines") into freebayes_concat
		set group, file("*versions.yml") into ch_freebayes_versions

	script:
		if (params.onco) {
			"""
			freebayes -f $genome_file --pooled-continuous --pooled-discrete -t $params.intersect_bed --min-repeat-entropy 1 -F 0.03 $bam > ${id}.freebayes.vcf
			vcfbreakmulti ${id}.freebayes.vcf > ${id}.freebayes.multibreak.vcf
			bcftools norm -m-both -c w -O v -f $genome_file -o ${id}.freebayes.multibreak.norm.vcf ${id}.freebayes.multibreak.vcf
			vcfanno_linux64 -lua $params.VCFANNO_LUA $params.vcfanno ${id}.freebayes.multibreak.norm.vcf > ${id}.freebayes.multibreak.norm.anno.vcf
			grep ^# ${id}.freebayes.multibreak.norm.anno.vcf > ${id}.freebayes.multibreak.norm.anno.path.vcf
			grep -v ^# ${id}.freebayes.multibreak.norm.anno.vcf | grep -i pathogenic > ${id}.freebayes.multibreak.norm.anno.path.vcf2
			cat ${id}.freebayes.multibreak.norm.anno.path.vcf ${id}.freebayes.multibreak.norm.anno.path.vcf2 > ${id}.freebayes.multibreak.norm.anno.path.vcf3
			filter_freebayes.pl ${id}.freebayes.multibreak.norm.anno.path.vcf3 > ${id}.pathfreebayes.lines

			${freebayes_version(task)}
			"""
		}
		else {
			"""
			touch "${id}.pathfreebayes.lines"

			${freebayes_version(task)}
			"""
		}

	stub:
		"""
		touch "${id}.pathfreebayes.lines"

		${freebayes_version(task)}
		"""
}
def freebayes_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g')
	    vcflib: 1.0.9
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	    vcfanno: \$(echo \$(vcfanno_linux64 2>&1 | grep version | cut -f3 -d' ')  )
	END_VERSIONS
	"""
}

/////////////// MITOCHONDRIA SNV CALLING ///////////////
///////////////                          ///////////////

// create an MT BAM file
process fetch_MTseqs {
	cpus 2
	memory '10GB'
	time '1h'
	tag "$id"
	publishDir "${OUTDIR}/bam", mode: 'copy', overwrite: 'true', pattern: '*.bam*'

	when:
		params.antype == "wgs"

	input:
		set group, id, file(bam), file(bai) from bam_mito.mix(bam_mito_choice)

    output:
        set group, id, file ("${id}_mito.bam"), file("${id}_mito.bam.bai") into mutserve_bam, eklipse_bam, qc_mito_bam
		set group, file("${group}_mtbam.INFO") into mtBAM_INFO
		set group, file("*versions.yml") into ch_fetch_mt_seqs_versions

	script:
		"""
		sambamba view -f bam $bam M > ${id}_mito.bam
		samtools index -b ${id}_mito.bam
		echo "mtBAM	$id	/access/${params.subdir}/bam/${id}_mito.bam" > ${group}_mtbam.INFO

		${fetch_MTseqs_version(task)}
		"""

	stub:
		"""
		touch "${id}_mito.bam"
		touch "${id}_mito.bam.bai"
		touch "${group}_mtbam.INFO"

		${fetch_MTseqs_version(task)}
		"""
}
def fetch_MTseqs_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sambamba: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
	    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
	END_VERSIONS	
	"""
}


process sentieon_mitochondrial_qc {

    // Fetch mitochondrial coverage statistics
    // Calculate mean_coverage and pct_above_500x
    
    cpus 30
    memory '20 GB'
	tag "$id"
	time '2h'
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'
	container = "${params.container_sentieon}"

	when:
	    params.antype == "wgs"
    
	input:
        set group, id, file(bam), file(bai) from qc_mito_bam

	output:
    	set group, id, file("${id}_mito_coverage.tsv") into qc_mito
		set group, file("*versions.yml") into ch_sentieon_mitochondrial_qc

	script:	
		"""
		sentieon driver \\
			-r $genome_file \\
			-t ${task.cpus} \\
			-i $bam \\
			--algo CoverageMetrics \\
			--cov_thresh 500 \\
			mt_cov_metrics.txt

		head -1 mt_cov_metrics.txt.sample_interval_summary > "${id}_mito_coverage.tsv"
		grep "^M" mt_cov_metrics.txt.sample_interval_summary >> "${id}_mito_coverage.tsv"
		${sentieon_mitochondrial_qc_version(task)}
		"""

	stub:
		"""
		touch "${id}_mito_coverage.tsv"
		${sentieon_mitochondrial_qc_version(task)}
		"""
}
def sentieon_mitochondrial_qc_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}

process build_mitochondrial_qc_json {
    memory '1 GB'
    cpus 2
    tag "$id"
    time "1h"

    input:
        set group, id, file(mito_qc_file) from qc_mito
    output:
        set group, id, file("${id}_mito_qc.json") into qc_mito_json
    
	script:
		"""
		mito_tsv_to_json.py ${mito_qc_file} > "${id}_mito_qc.json"
		"""
	stub:
		"""
		touch "${id}_mito_qc.json"
		"""
}

    
// gatk FilterMutectCalls in future if FPs overwhelms tord/sofie/carro
process run_mutect2 {
	cpus 4
	memory '50 GB'
	time '1h'
	tag "$group"
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf'

	when:
		!params.onco
	
	input:
		set group, id, file(bam), file(bai) from mutserve_bam.groupTuple()

	output:
		set group, id, file("${group}.mutect2.vcf") into ms_vcfs_1, ms_vcfs_2
		set group, file("*versions.yml") into ch_run_mutect2_versions

	script:
		bams = bam.join(' -I ')

		"""
		source activate gatk4-env
		gatk Mutect2 \
		--mitochondria-mode \
		-R $params.genome_file \
		-L M \
		-I $bams \
		-O ${group}.mutect2.vcf

		${run_mutect2_version(task)}
		"""

	stub:
		"""
		source activate gatk4-env
		touch "${group}.mutect2.vcf"

		${run_mutect2_version(task)}
		"""
}
def run_mutect2_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
	END_VERSIONS
	"""
}

// split and left-align variants
process split_normalize_mito {
	cpus 2
	memory '1GB'
	time '1h'

	input:
		set group, id, file(ms_vcf) from ms_vcfs_1
		set g2, id2, sex, type from meta_mutect2.groupTuple()

	output:
		set group, file("${group}.mutect2.breakmulti.filtered5p.0genotyped.proband.vcf") into adj_vcfs
		set group, file("*versions.yml") into ch_split_normalize_mito_versions

	script:
		proband_idx = type.findIndexOf{ it == "proband" }

		"""
		grep -vP "^M\\s+955" $ms_vcf > ${ms_vcf}.fix
		bcftools norm -m-both -o ${ms_vcf}.breakmulti ${ms_vcf}.fix
		bcftools sort ${ms_vcf}.breakmulti | bgzip > ${ms_vcf}.breakmulti.fix
		tabix -p vcf ${ms_vcf}.breakmulti.fix
		bcftools norm -f $params.rCRS_fasta -o ${ms_vcf.baseName}.adjusted.vcf ${ms_vcf}.breakmulti.fix
		bcftools view -i 'FMT/AF[*]>0.05' ${ms_vcf.baseName}.adjusted.vcf -o ${group}.mutect2.breakmulti.filtered5p.vcf
		bcftools filter -S 0 --exclude 'FMT/AF[*]<0.05' ${group}.mutect2.breakmulti.filtered5p.vcf -o ${group}.mutect2.breakmulti.filtered5p.0genotyped.vcf
		filter_mutect2_mito.pl ${group}.mutect2.breakmulti.filtered5p.0genotyped.vcf ${id2[proband_idx]} > ${group}.mutect2.breakmulti.filtered5p.0genotyped.proband.vcf

		${split_normalize_mito_version(task)}
		"""

	stub:
		"""
		touch "${group}.mutect2.breakmulti.filtered5p.0genotyped.proband.vcf"
		${split_normalize_mito_version(task)}
		"""
}
def split_normalize_mito_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcflib: 1.0.9
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	END_VERSIONS
	"""
}

// use python tool HmtNote for annotating vcf
// future merging with diploid genome does not approve spaces in info-string
process run_hmtnote {
	cpus 2
	memory '5GB'
	time '1h'

	input:
		set group, file(adj_vcf) from adj_vcfs

	output:
		set group, file("${group}.fixinfo.vcf") into mito_diplod_vep
		set group, file("*versions.yml") into ch_run_hmtnote_versions

	script:
		"""
		source activate tools
		hmtnote annotate ${adj_vcf} ${group}.hmtnote --offline
		grep ^# ${group}.hmtnote > ${group}.fixinfo.vcf
		grep -v ^# ${group}.hmtnote | sed 's/ /_/g' >> ${group}.fixinfo.vcf

		${run_hmtnote_version(task)}
		"""

	stub:
		"""
		source activate tools
		touch "${group}.fixinfo.vcf"

		${run_hmtnote_version(task)}
		"""
}
def run_hmtnote_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    hmtnote: \$(echo \$(hmtnote --version 2>&1) | sed 's/^.*hmtnote, version //; s/Using.*\$//' )
	END_VERSIONS
	"""
}

// run haplogrep 2 on resulting vcf
process run_haplogrep {
	time '1h'
	memory '50 GB'
	cpus 2
	publishDir "${OUTDIR}/plots/mito", mode: 'copy', overwrite: 'true', pattern: '*.png'

	input:
		set group, id, file(ms_vcf) from ms_vcfs_2

	output:
		file("${group}.haplogrep.png")
		set group, file("${group}_haplo.INFO") into haplogrep_INFO
		set group, file("*versions.yml") into ch_run_haplogrep_versions

	shell:
		version_str = run_haplogrep_version(task)
		'''
		for sample in `bcftools query -l !{ms_vcf}`; do 
			bcftools view -c1 -Oz -s $sample -o $sample.vcf.gz !{ms_vcf}
			java  -Xmx16G -Xms16G -jar /opt/bin/haplogrep.jar classify \
			--in $sample.vcf.gz\
			--out $sample.hg2.vcf \
			--format vcf \
			--lineage 1
			dot $sample.hg2.vcf.dot -Tps2 > $sample.hg2.vcf.ps
			gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dGraphicsAlphaBits=4 -r1200 -dDownScaleFactor=3 -sOutputFile=${sample}.hg2.vcf.png ${sample}.hg2.vcf.ps
		done
		montage -mode concatenate -tile 3x1 *.png !{group}.haplogrep.png
		echo "IMG haplogrep !{params.accessdir}/plots/mito/!{group}.haplogrep.png" > !{group}_haplo.INFO
		
		echo "!{version_str}" > "!{task.process}_versions.yml"
		'''

	stub:
		version_str = run_haplogrep_version(task)
		"""
		touch "${group}.haplogrep.png"
		touch "${group}_haplo.INFO"

		echo "${version_str}" > "${task.process}_versions.yml"
		"""
}
def run_haplogrep_version(task) {
	// This docstring looks different from others as it is used within the shell '''
	// If spaces similarly to the others, this leads to additional whitespace above and below the version text
	"""${task.process}:
	    haplogrep: \$(echo \$(java -jar /opt/bin/haplogrep.jar classify 2>&1) | sed "s/htt.*Classify v// ; s/ .*//")
	    montage: \$(echo \$(gm -version 2>&1) | head -1 | sed -e "s/GraphicsMagick //" | cut -d" " -f1 )"""
}

// use eKLIPse for detecting mitochondrial deletions
process run_eklipse {

	tag "$id"
	cpus 2
	// in rare cases with samples above 50 000x this can peak at 500+ GB of VMEM. Add downsampling!
	memory '100GB'
	time '60m'
	publishDir "${OUTDIR}/plots/mito", mode: 'copy', overwrite: 'true', pattern: '*.txt'
	publishDir "${OUTDIR}/plots/mito", mode: 'copy', overwrite: 'true', pattern: '*.png'

	input:
		set group, id, file(bam), file(bai), sex, type from eklipse_bam.join(meta_eklipse, by: [0,1])
	
	output:
		set file("*.png"), file("${id}.hetplasmid_frequency.txt")
		set group, file("${id}_eklipse.INFO")  optional true into eklipse_INFO
		set group, file("*versions.yml") into ch_run_eklipse_versions

	script:
		yml_info_command = ""
		if (type == "proband") {
			yml_info_command = "echo 'IMG eklipse ${params.accessdir}/plots/mito/${id}_eklipse.png' > ${id}_eklipse.INFO"
		}
		"""
		source activate htslib10
		echo "${bam}\tsample" > infile.txt
		python /eKLIPse/eKLIPse.py \
		-in infile.txt \
		-ref /eKLIPse/data/NC_012920.1.gb
		mv eKLIPse_*/eKLIPse_deletions.csv ./${id}_deletions.csv
		mv eKLIPse_*/eKLIPse_genes.csv ./${id}_genes.csv
		mv eKLIPse_*/eKLIPse_sample.png ./${id}_eklipse.png
		hetplasmid_frequency_eKLIPse.pl --bam ${bam} --in ${id}_deletions.csv
		mv hetplasmid_frequency.txt ${id}.hetplasmid_frequency.txt
		$yml_info_command

		${run_eklipse_version(task)}
		"""

	stub:
		yml_info_command = ""
		if (type == "proband") {
			yml_info_command = "echo 'IMG eklipse ${params.accessdir}/plots/mito/${id}_eklipse.png' > ${id}_eklipse.INFO"
		}
		"""
		source activate htslib10
		touch "${id}.hetplasmid_frequency.txt"
		touch "${id}.png"
		$yml_info_command

		${run_eklipse_version(task)}
		"""
}
def run_eklipse_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    eklipse: 1.8
	END_VERSIONS	
	"""
}

//eklipseM_INFO.collectFile(name: "eklipse.INFO").set{ eklipse_INFO }

// Splitting & normalizing variants, merging with Freebayes/Mutect2, intersecting against exome/clinvar introns
process split_normalize {
	cpus 2
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	tag "$group"
	memory '50 GB'
	time '1h'
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'

	when:
		params.annotate

	input:
		set group, id, file(vcf), file(idx), file(vcfconcat) from combined_vcf.join(mito_diplod_vep.mix(freebayes_concat))

	output:
		set group, file("${group}.norm.uniq.DPAF.vcf") into split_norm, vcf_gnomad
		set group, id, file("${group}.intersected.vcf"), file("${group}.multibreak.vcf") into split_vep, split_cadd, vcf_cnvkit
		set group, file("*versions.yml") into ch_split_normalize_versions

	script:
	id = id[0]
	// rename M to MT because genmod does not recognize M
	if (params.onco || params.assay == "modycf") {
		"""
		cat $vcf $vcfconcat > ${id}.concat.freebayes.vcf
		vcfbreakmulti ${id}.concat.freebayes.vcf > ${group}.multibreak.vcf
		bcftools norm -m-both -c w -O v -f $genome_file -o ${group}.norm.vcf ${group}.multibreak.vcf
		bcftools sort ${group}.norm.vcf | vcfuniq > ${group}.norm.uniq.vcf
		wgs_DPAF_filter.pl ${group}.norm.uniq.vcf > ${group}.norm.uniq.DPAF.vcf
		bedtools intersect \
			-a ${group}.norm.uniq.DPAF.vcf \\
			-b $params.intersect_bed \\
			-u -header > ${group}.intersected.vcf
		
		${split_normalize_version(task)}
		"""
	}

	else {
		"""
		vcfbreakmulti ${vcf} > ${group}.multibreak.vcf
		bcftools norm -m-both -c w -O v -f $genome_file -o ${group}.norm.vcf ${group}.multibreak.vcf
		bcftools sort ${group}.norm.vcf | vcfuniq > ${group}.norm.uniq.vcf
		wgs_DPAF_filter.pl ${group}.norm.uniq.vcf > ${group}.norm.uniq.DPAF.vcf
		bedtools intersect \\
			-a ${group}.norm.uniq.DPAF.vcf \\
			-b $params.intersect_bed \\
			-u -header > ${group}.intersected_diploid.vcf
		java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar MergeVcfs \
		I=${group}.intersected_diploid.vcf I=$vcfconcat O=${group}.intersected.vcf
		sed 's/^M/MT/' -i ${group}.intersected.vcf
		sed 's/ID=M,length/ID=MT,length/' -i ${group}.intersected.vcf

		${split_normalize_version(task)}
		"""
	}

	stub:
		"""
		touch "${group}.norm.uniq.DPAF.vcf"
		touch "${group}.intersected.vcf"
		touch "${group}.multibreak.vcf"

		${split_normalize_version(task)}
		"""
}
def split_normalize_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcflib: 1.0.9
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	    bedtools: \$(echo \$(bedtools --version 2>&1) | sed -e "s/^.*bedtools v//" )
	    merge-vcfs: \$(echo \$(java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar MergeVcfs --version 2>&1 | sed 's/-SNAPSHOT//'))
	END_VERSIONS	
	"""
}

/////////////// Collect QC into single file ///////////////

process merge_qc_json {   
    cpus 2
    errorStrategy 'retry'
    maxErrors 5
    publishDir "${OUTDIR}/qc", mode: 'copy' , overwrite: 'true', pattern: '*.QC'
    tag "$id"
    time '1h'
	memory '1 GB'

    input:
        set group, id, file(qc) from qc_cdm.mix(qc_mito_json).groupTuple(by: [0,1])

    output:
        set id, file("${id}.QC") into qc_cdm_merged

    script:
        qc_json_files = qc.join(' ')
		"""
		merge_json_files.py ${qc_json_files} > ${id}.QC
		"""
	
	stub:
		"""
		touch "${id}.QC"
		"""
}   
    
// Load QC data into CDM (via middleman)
process qc_to_cdm {
	cpus 2
	errorStrategy 'retry'
	maxErrors 5
	publishDir "${CRONDIR}/qc", mode: 'copy' , overwrite: 'true'
	tag "$id"
	time '1h'

	when:
		!params.noupload
	
	input:
		set id, file(qc), diagnosis, r1, r2 from qc_cdm_merged.join(qc_extra)

	output:
		file("${id}.cdm") into cdm_done

	script:
		parts = r1.split('/')
		idx =  parts.findIndexOf {it ==~ /......_......_...._........../}
		rundir = parts[0..idx].join("/")
		"""
		echo "--run-folder $rundir --sample-id $id --subassay $diagnosis --assay $params.assay --qc ${OUTDIR}/qc/${id}.QC" > ${id}.cdm
		"""
}
    
    
process annotate_vep {
	container = "${params.container_vep}"
	cpus 30
	tag "$group"
	memory '50 GB'
	time '5h'

	input:
		set group, id, file(vcf), idx from split_vep.mix(annotate_only_vep)

	output:
		set group, file("${group}.vep.vcf") into vep
		set group, file("*versions.yml") into ch_annotate_vep_versions

	script:
		"""
		vep \\
			-i ${vcf} \\
			-o ${group}.vep.vcf \\
			--offline \\
			--sift b --polyphen b --ccds --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --af --af_1kg --max_af --pubmed --uniprot --mane --tsl --appris --variant_class --gene_phenotype --mirna \\
			--merged \\
			--vcf \\
			--no_stats \\
			--synonyms $params.VEP_SYNONYMS \\
			--fork ${task.cpus} \\
			--force_overwrite \\
			--fasta $params.VEP_FASTA \\
			--dir_cache $params.VEP_CACHE \\
			--dir_plugins $params.VEP_PLUGINS \\
			--distance $params.VEP_TRANSCRIPT_DISTANCE \\
			-cache \\
			--plugin CADD,$params.CADD \\
			--plugin LoFtool \\
			--plugin MaxEntScan,$params.MAXENTSCAN,SWA,NCSS \\
			--plugin dbNSFP,$params.DBNSFP,transcript_match=1,REVEL_score,REVEL_rankscore \\
			-custom $params.GNOMAD_EXOMES,gnomADe,vcf,exact,0,AF_popmax,AF,popmax \\
			-custom $params.GNOMAD_GENOMES,gnomADg,vcf,exact,0,AF_popmax,AF,popmax \\
			-custom $params.GNOMAD_MT,gnomAD_mt,vcf,exact,0,AF_hom,AF_het \\
			-custom $params.PHYLOP,phyloP100way,bigwig \\
			-custom $params.PHASTCONS,phastCons,bigwig

		${annotate_vep_version(task)}
		"""

	stub:
		"""
		touch "${group}.vep.vcf"
		${annotate_vep_version(task)}
		"""
}
def annotate_vep_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
	END_VERSIONS
	"""
}


// gene, clinvar, loqusdb, enigma(onco)
process vcfanno {
	memory '1GB'
	time '1h'
	errorStrategy 'retry'
	maxErrors 5
	cpus 2

	input:
		set group, file(vcf) from vep

	output:
		set group, file("${group}.clinvar.loqusdb.gene.vcf") into vcfanno_vcf
		set group, file("*versions.yml") into ch_vcfanno_versions

	script:
		"""
		vcfanno_linux64 -lua $params.VCFANNO_LUA $params.vcfanno $vcf > ${group}.clinvar.loqusdb.gene.vcf
		${vcfanno_version(task)}
		"""

	stub:
		"""
		touch "${group}.clinvar.loqusdb.gene.vcf"
		${vcfanno_version(task)}
		"""
}
def vcfanno_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcfanno: \$(echo \$(vcfanno_linux64 2>&1 | grep version | cut -f3 -d' ')  )
	END_VERSIONS	
	"""
}


// Extracting most severe consequence: 
// Modifying annotations by VEP-plugins, and adding to info-field: 
// Modifying CLNSIG field to allow it to be used by genmod score properly:
process modify_vcf {
	cpus 2
	tag "$group"
	memory '1 GB'
	time '1h'

	input:
		set group, file(vcf) from vcfanno_vcf

	output:
		set group, file("${group}.mod.vcf") into mod_vcf

	script:
		"""
		modify_vcf_scout.pl $vcf > ${group}.mod.vcf
		"""

	stub:
		"""
		touch "${group}.mod.vcf"
		"""
} 


// Marking splice INDELs: 
process mark_splice {
	cpus 2
	tag "$group"
	memory '1 GB'
	time '1h'

	input:
		set group, file(vcf) from mod_vcf

	output:
		set group, file("${group}.marksplice.vcf") into splice_marked

	script:
		"""
		/opt/bin/mark_spliceindels.pl $vcf > ${group}.marksplice.vcf
		"""

	stub:
		"""
		touch "${group}.marksplice.vcf"
		"""
}

// Extract all INDELs from VCF for CADD annotation
process extract_indels_for_cadd {
	cpus 2
	tag "$group"
	memory '1 GB'
	time '1h'

	input:
		set group, id, file(vcf), idx from split_cadd.mix(annotate_only_cadd)
	
	output:
		set group, file("${group}.only_indels.vcf") into indel_cadd_vep
		set group, file("*versions.yml") into ch_extract_indels_for_cadd_versions

	script:
		"""
		bcftools view $vcf -V snps -o ${group}.only_indels.vcf
		${extract_indels_for_cadd_version(task)}
		"""

	stub:
		"""
		touch "${group}.only_indels.vcf"
		${extract_indels_for_cadd_version(task)}
		""" 
}
def extract_indels_for_cadd_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS	
	"""
}

// Annotate Indels with VEP+Gnomad genomes. Filter variants below threshold
process indel_vep {
	cpus 5
	container = "${params.container_vep}"
	tag "$group"
	memory '10 GB'
	time '3h'

	input:
		set group, file(vcf) from indel_cadd_vep

	output:
		set group, file("${group}.only_indels.vep.filtered.vcf") into indel_cadd_vcf
		set group, file("*versions.yml") into ch_indel_vep_versions

	script:
		"""
		vep \\
			-i $vcf \\
			-o ${group}.only_indels.vep.vcf \\
			--offline \\
			--cache \\
			--merged \\
			--vcf \\
			--synonyms $params.VEP_SYNONYMS \\
			--fasta $params.VEP_FASTA \\
			-custom $params.GNOMAD_GENOMES,gnomADg,vcf,exact,0,AF \\
			-custom $params.GNOMAD_MT,gnomAD_mt,vcf,exact,0,AF_hom,AF_het \\
			--dir_cache $params.VEP_CACHE \\
			--force_overwrite \\
			--no_stats \\
			--fork ${task.cpus}
		filter_indels.pl ${group}.only_indels.vep.vcf > ${group}.only_indels.vep.filtered.vcf
		${indel_vep_version(task)}
		"""

	stub:
		"""
		touch "${group}.only_indels.vep.filtered.vcf"
		${indel_vep_version(task)}
		"""
}
def indel_vep_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
	END_VERSIONS	
	"""
}

// Calculate CADD scores for all indels
process calculate_indel_cadd {
	cpus 2
	container = "${params.container_cadd}"
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'
	tag "$group"
	memory '20 GB'
	time '3h'

	input:
		set group, file(vcf) from indel_cadd_vcf

	output:
		set group, file("${group}.indel_cadd.gz") into indel_cadd
		set group, file("*versions.yml") into ch_calculate_indel_cadd_versions

	script:
		"""
		/CADD-scripts/CADD.sh -c ${task.cpus} -g GRCh38 -o ${group}.indel_cadd.gz $vcf
		${calculate_indel_cadd_version(task)}
		"""

	stub:
		"""
		touch "${group}.indel_cadd.gz"
		${calculate_indel_cadd_version(task)}
		"""
}
def calculate_indel_cadd_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    cadd: \$(echo \$(/CADD-scripts/CADD.sh -v 2>&1) | sed -e "s/^.*CADD-v// ; s/ (c).*//")
	END_VERSIONS	
	"""
}

process bgzip_indel_cadd {

	tag "$group"
	cpus 4
	memory '1 GB'
	time '5m'
	container = "${params.container_bcftools}"

	input:
		set group, file(cadd_scores) from indel_cadd
	
	output:
		set group, file("${group}.cadd.gz"), file("${group}.cadd.gz.tbi") into indel_cadd_bgzip
		set group, file("*versions.yml") into ch_bgzip_indel_cadd_versions
	
	script:
		"""
		gunzip -c ${cadd_scores} > "${group}.cadd"
		bgzip -@ ${task.cpus} "${group}.cadd"
		tabix -p vcf "${group}.cadd.gz"
		${bgzip_indel_cadd_version(task)}
		"""
	
	stub:
		"""
		touch "${group}.cadd.gz"
		touch "${group}.cadd.gz.tbi"
		${bgzip_indel_cadd_version(task)}
		"""
}
def bgzip_indel_cadd_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS
	"""
}

// Add the calculated indel CADDs to the vcf
process add_cadd_scores_to_vcf {
	cpus 4
	tag "$group"
	memory '1 GB'
	time '5m'
	container = "${params.container_genmod}"

	input: 
		set group, file(vcf), file(cadd_scores), file(cadd_scores_tbi) from splice_marked.join(indel_cadd_bgzip)

	output:
		set group, file("${group}.cadd.vcf") into ma_vcf, fa_vcf, base_vcf
		set group, file("*versions.yml") into ch_add_cadd_scores_to_vcf_versions

	script:
		"""
		genmod annotate --cadd-file ${cadd_scores} ${vcf} > ${group}.cadd.vcf

		${add_cadd_scores_to_vcf_version(task)}
		"""

	stub:
		"""
		touch "${group}.cadd.vcf"
		${add_cadd_scores_to_vcf_version(task)}
		"""
}
def add_cadd_scores_to_vcf_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    genmod: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
	END_VERSIONS
	"""
}


// # Annotating variant inheritance models:
process inher_models {
	tag "$group"
	cpus 3
	memory '80 GB'
	time '1h'
	container = "${params.container_genmod}"

	input:
		set group, file(vcf), type, file(ped) from base_vcf.mix(ma_vcf, fa_vcf).join(ped_inher.mix(ped_inher_ma,ped_inher_fa))

	output:
		set group, type, file("${group}.models.vcf") into inhermod
		set group, file("*versions.yml") into ch_inher_models_versions

	script:
		"""
		genmod models $vcf -p ${task.cpus} -f $ped > ${group}.models.vcf
		${inher_models_version(task)}
		"""

	stub:
		"""
		touch "${group}.models.vcf"
		${inher_models_version(task)}
		"""
}
def inher_models_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    genmod: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
	END_VERSIONS	
	"""
}


// Scoring variants: 
// Adjusting compound scores: 
// Sorting VCF according to score: 
process genmodscore {
	tag "$group"
	cpus 2
	memory '20 GB'
	time '1h'
	container = "${params.container_genmod}"

	input:
		set group, type, file(vcf) from inhermod

	output:
		set group, type, file("${group_score}.scored.vcf") into scored_vcf
		set group, file("*versions.yml") into ch_genmodscore_versions

	script:
		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group

		if ( mode == "family" && params.antype == "wgs" ) {
			"""
			genmod score -i $group_score -c $params.rank_model -r $vcf -o ${group_score}.only_rankscore.vcf
			genmod compound \
				--threshold ${params.genmod_compound_trio_threshold} \
				--penalty ${params.genmod_compound_trio_penalty} \
				-o ${group_score}.with_compounds.vcf \
				${group_score}.only_rankscore.vcf
			sed 's/RankScore=${group}:/RankScore=${group_score}:/g' -i ${group_score}.with_compounds.vcf
			genmod sort -p -f $group_score ${group_score}.with_compounds.vcf -o ${group_score}.scored.vcf

			${genmodscore_version(task)}
			"""
		}
		else {
			"""
			genmod score -i $group_score -c $params.rank_model_s -r $vcf -o ${group_score}.only_rankscore.vcf

			# To get compounds without applying rank score penalty
			genmod compound \
				--penalty 0 \
				-o ${group_score}.with_compounds.vcf \
				${group_score}.only_rankscore.vcf

			genmod sort \
				-p \
				-f $group_score \
				-o ${group_score}.scored.vcf \
				${group_score}.with_compounds.vcf

			${genmodscore_version(task)}
			"""
		}

	stub:
		group_score = group
		"""
		touch "${group_score}.scored.vcf"
		${genmodscore_version(task)}
		"""
}
def genmodscore_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    genmod: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
	END_VERSIONS	
	"""
}

// Bgzipping and indexing VCF: 
process vcf_completion {
	cpus 16
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
	tag "$group"
	time '1h'
	memory '5 GB'

	input:
		set group, type, file(vcf) from scored_vcf

	output:
		set group, type, file("${group_score}.scored.vcf.gz"), file("${group_score}.scored.vcf.gz.tbi") into vcf_peddy, snv_sv_vcf, snv_sv_vcf_ma, snv_sv_vcf_fa, vcf_loqus
		set group, file("${group}_snv.INFO") into snv_INFO
		set group, file("*versions.yml") into ch_vcf_completion_versions

	script:
		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group

		"""
		sed 's/^MT/M/' -i $vcf
		sed 's/ID=MT,length/ID=M,length/' -i $vcf
		bgzip -@ ${task.cpus} $vcf -f
		tabix ${vcf}.gz -f
		echo "SNV	$type	${params.accessdir}/vcf/${group_score}.scored.vcf.gz" > ${group}_snv.INFO

		${vcf_completion_version(task)}
		"""

	stub:
		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group

		"""
		touch "${group_score}.scored.vcf.gz"
		touch "${group_score}.scored.vcf.gz.tbi"
		touch "${group}_snv.INFO"

		${vcf_completion_version(task)}
		"""
}
def vcf_completion_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	END_VERSIONS	
	"""
}

process peddy {

	publishDir "${OUTDIR}/ped", mode: 'copy' , overwrite: 'true', pattern: '*.ped'
	publishDir "${OUTDIR}/ped", mode: 'copy' , overwrite: 'true', pattern: '*.csv'

	cpus 4
	tag "$group"
	time '1h'
	memory '20GB'

	when:
		!params.annotate_only && params.run_peddy

	input:
		set group, type, file(vcf), file(idx), type, file(ped) from vcf_peddy.join(ped_peddy)

	output:
		set file("${group}.ped_check.csv"),file("${group}.peddy.ped"), file("${group}.sex_check.csv") into peddy_files
		set group, file("${group}_peddy.INFO") into peddy_INFO
		set group, file("*versions.yml") into ch_peddy_versions

	script:
		"""
		source activate py3-env
		python -m peddy --sites hg38 -p ${task.cpus} $vcf $ped --prefix $group
		echo "PEDDY	${params.accessdir}/ped/${group}.ped_check.csv,${params.accessdir}/ped/${group}.peddy.ped,${params.accessdir}/ped/${group}.sex_check.csv" > ${group}_peddy.INFO

		${peddy_version(task)}
		"""

	stub:
		"""
		source activate py3-env
		touch "${group}.ped_check.csv"
		touch "${group}.peddy.ped"
		touch "${group}.sex_check.csv"
		touch "${group}_peddy.INFO"

		${peddy_version(task)}
		"""
}
def peddy_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    peddy: \$(echo \$(python -m peddy --version 2>&1) | sed 's/^.*peddy, version //')
	END_VERSIONS	
	"""
}

// Extract all variants (from whole genome) with a gnomAD af > x%
process fastgnomad {
	cpus 2
	memory '40 GB'
	tag "$group"
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	time '2h'

	when:
		params.antype == "wgs"

	input:
		set group, file(vcf) from vcf_gnomad

	output:
		set group, file("${group}.SNPs.vcf") into vcf_upd, vcf_roh, vcf_pod

	script:
		"""
		gzip -c $vcf > ${vcf}.gz
		annotate -g $params.FASTGNOMAD_REF -i ${vcf}.gz > ${group}.SNPs.vcf
		"""

	stub:
		"""
		touch "${group}.SNPs.vcf"
		"""
}


// Call UPD regions from SNP vcf
process upd {
	tag "$group"
	time '1h'
	memory '1 GB'
	cpus 2

	input:
		set group, file(vcf) from vcf_upd
		set group, id, sex, mother, father, phenotype, diagnosis, type, assay, clarity_sample_id, ffpe, analysis from meta_upd.filter{ item -> item[7] == 'proband' }

	output:
		file("upd.bed") into upd_plot
		set group, file("upd.sites.bed") into upd_table
		set group, file("*versions.yml") into ch_upd_versions

	script:
		if( mode == "family" && trio == true ) {
			"""
			upd --vcf $vcf --proband $id --mother $mother --father $father --af-tag GNOMADAF regions > upd.bed
			upd --vcf $vcf --proband $id --mother $mother --father $father --af-tag GNOMADAF sites > upd.sites.bed

			${upd_version(task)}
			"""
		}
		else {
			"""
			touch "upd.bed"
			touch "upd.sites.bed"

			${upd_version(task)}
			"""
		}

	stub:
		"""
		touch "upd.bed"
		touch "upd.sites.bed"
		
		${upd_version(task)}
		"""
}
def upd_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    upd: \$(echo \$(upd --version 2>&1))
	END_VERSIONS
	"""
}


process upd_table {
	publishDir "${OUTDIR}/plots", mode: 'copy' , overwrite: 'true'
	tag "$group"
	time '1h'
	memory '1 GB'
	cpus 2

	input:
		set group, file(upd_sites) from upd_table

	output:
		file("${group}.UPDtable.xls")

	when:
		mode == "family" && trio == true

	script:
		"""
		upd_table.pl $upd_sites > ${group}.UPDtable.xls
		"""

	stub:
		"""
		touch "${group}.UPDtable.xls"
		"""
}


// Call ROH regions from SNP vcf
process roh {
	tag "$group"
	time '1h'
	memory '1 GB'
	cpus 2

	input:
		set group, file(vcf) from vcf_roh

	output:
		set group, file("roh.txt") into roh_plot
		set group, file("*versions.yml") into ch_roh_versions

	script:
		"""
		bcftools roh --rec-rate 1e-9 --AF-tag GNOMADAF ${vcf} -o roh.txt
		${roh_version(task)}
		"""

	stub:
		"""
		touch "roh.txt"
		${roh_version(task)}
		"""
}
def roh_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS	
	"""
}

// Create coverage profile using GATK
process gatkcov {
	publishDir "${OUTDIR}/cov", mode: 'copy' , overwrite: 'true', pattern: '*.tsv'
	tag "$group"
	cpus 2
	memory '80 GB'
	time '5h'

	input:
		set id, group, file(bam), file(bai), gr, sex, type from cov_bam.mix(cov_bam_choice).join(meta_gatkcov, by:1)

	output:
		set group, id, type, sex, file("${id}.standardizedCR.tsv"), file("${id}.denoisedCR.tsv") into cov_plot, cov_gens
		set group, file("*versions.yml") into ch_gatkcov_versions

	when:
		params.gatkcov

	script:
		"""
		source activate gatk4-env

		gatk CollectReadCounts \\
			-I $bam -L $params.COV_INTERVAL_LIST \\
			--interval-merging-rule OVERLAPPING_ONLY -O ${bam}.hdf5

		gatk --java-options "-Xmx30g" DenoiseReadCounts \\
			-I ${bam}.hdf5 --count-panel-of-normals ${PON[sex]} \\
			--standardized-copy-ratios ${id}.standardizedCR.tsv \\
			--denoised-copy-ratios ${id}.denoisedCR.tsv

		gatk PlotDenoisedCopyRatios \\
			--standardized-copy-ratios ${id}.standardizedCR.tsv \\
			--denoised-copy-ratios ${id}.denoisedCR.tsv \\
			--sequence-dictionary $params.GENOMEDICT \\
			--minimum-contig-length 46709983 --output . --output-prefix $id
		
		${gatkcov_version(task)}
		"""

	stub:
		"""
		source activate gatk4-env
		touch "${id}.standardizedCR.tsv"
		touch "${id}.denoisedCR.tsv"

		${gatkcov_version(task)}
		"""
}
def gatkcov_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
	END_VERSIONS	
	"""
}


// Plot ROH, UPD and coverage in a genomic overview plot
process overview_plot {
	publishDir "${OUTDIR}/plots", mode: 'copy' , overwrite: 'true', pattern: "*.png"
	tag "$group"
	time '1h'
	memory '5 GB'
	cpus 2

	input:
		file(upd) from upd_plot
		set group, file(roh) from roh_plot
		set group, id, type, sex, file(cov_stand), file(cov_denoised) from cov_plot.groupTuple()

	output:
		file("${group}.genomic_overview.png")
		set group, file("${group}_oplot.INFO") into oplot_INFO

	script:
		proband_idx = type.findIndexOf{ it == "proband" }

	script:
		"""
		genome_plotter.pl --dict $params.GENOMEDICT \\
			--sample ${id[proband_idx]} \\
			--upd $upd \\
			--roh $roh \\
			--sex ${sex[proband_idx]} \\
			--cov ${cov_denoised[proband_idx]} \\
			--out ${group}.genomic_overview.png
		echo "IMG overviewplot	${params.accessdir}/plots/${group}.genomic_overview.png" > ${group}_oplot.INFO 
		"""

	stub:
		"""
		touch "${group}.genomic_overview.png"
		touch "${group}_oplot.INFO"
		"""
}

process generate_gens_data {
	publishDir "${OUTDIR}/plot_data", mode: 'copy' , overwrite: 'true', pattern: "*.gz*"
	publishDir "${CRONDIR}/gens", mode: 'copy', overwrite: 'true', pattern: "*.gens"
	tag "$group"
	cpus 1
	time '3h'
	memory '5 GB'

	when:
		params.prepare_gens_data

	input:
		set id, group, file(gvcf), g, type, sex, file(cov_stand), file(cov_denoise) from gvcf_gens_choice.join(cov_gens, by:[1])

	output:
		set file("${id}.cov.bed.gz"), file("${id}.baf.bed.gz"), file("${id}.cov.bed.gz.tbi"), file("${id}.baf.bed.gz.tbi"), file("${id}.overview.json.gz")
		file("${id}.gens") into gens_middleman

	script:
		"""
		generate_gens_data.pl $cov_stand $gvcf $id $params.GENS_GNOMAD
		echo "gens load sample --sample-id $id --genome-build 38 --baf ${params.gens_accessdir}/${id}.baf.bed.gz --coverage ${params.gens_accessdir}/${id}.cov.bed.gz --overview-json ${params.gens_accessdir}/${id}.overview.json.gz" > ${id}.gens
		"""

	stub:
		"""
		touch "${id}.cov.bed.gz"
		touch "${id}.baf.bed.gz"
		touch "${id}.cov.bed.gz.tbi"
		touch "${id}.baf.bed.gz.tbi"
		touch "${id}.overview.json.gz"
		touch "${id}.gens"
		"""
}

// SV-calling //

// GATK panel+wgs //

process gatk_coverage {
	cpus 2
	memory '50GB'
	time '2h'
	container = "${params.container_gatk}"
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'
	tag "$id"   

	when:
		params.sv && params.gatkcnv

	input:
		set group, id, file(bam), file(bai) from bam_gatk.mix(bam_gatk_choice)

	output:
		set group, id, file("${id}.tsv") into call_ploidy, call_cnv
		set group, file("*versions.yml") into ch_gatk_coverage_versions

	script:
		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		gatk --java-options "-Xmx20g" CollectReadCounts \\
			-L $params.gatk_intervals \\
			-R $params.genome_file \\
			-imr OVERLAPPING_ONLY \\
			-I $bam \\
			--format TSV -O ${id}.tsv

		${gatk_coverage_version(task)}
		"""

	stub:
		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		touch "${id}.tsv"

		${gatk_coverage_version(task)}
		"""
}
def gatk_coverage_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
	END_VERSIONS	
	"""
}

process gatk_call_ploidy {
	cpus 10
	memory '50GB'
	time '2h'
	container = "${params.container_gatk}"
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'
	tag "$id"

	input:
		set group, id, file(tsv) from call_ploidy

	output:
		set group, id, file("ploidy.tar") into ploidy_to_cnvcall, ploidy_to_post
		set group, file("*versions.yml") into ch_gatk_call_ploidy_versions

	script:
		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		gatk --java-options "-Xmx20g" DetermineGermlineContigPloidy \\
			--model $params.ploidymodel \\
			-I $tsv \\
			-O ploidy/ \\
			--output-prefix $group
		tar -cvf ploidy.tar ploidy/

		${gatk_call_ploidy_version(task)}
		"""

	stub:
		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		touch "ploidy.tar"

		${gatk_call_ploidy_version(task)}
		"""
}
def gatk_call_ploidy_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
	END_VERSIONS	
	"""
}

process gatk_call_cnv {
	cpus 8
	memory '50GB'
	time '3h'
	container = "${params.container_gatk}"
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'
	tag "$id"

	input:
		set group, id, file(tsv), file(ploidy), i, refpart \
			from call_cnv.join(ploidy_to_cnvcall, by: [0,1]).combine(gatk_ref)

	output:
		set group, id, i, file("${group}_${i}.tar") into postprocessgatk
		set group, file("*versions.yml") into ch_gatk_call_cnv_versions

	script:
		"""
		export THEANO_FLAGS="base_compiledir=."
		set +u
		source activate gatk
		export HOME=/local/scratch
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		tar -xvf ploidy.tar
		mkdir ${group}_${i}
		gatk --java-options "-Xmx25g" GermlineCNVCaller \\
			--run-mode CASE \\
			-I $tsv \\
			--contig-ploidy-calls ploidy/${group}-calls/ \\
			--model ${refpart} \\
			--output ${group}_${i}/ \\
			--output-prefix ${group}_${i}
		tar -cvf ${group}_${i}.tar ${group}_${i}/
		
		${gatk_call_cnv_version(task)}
		"""

	stub:
		"""
		export THEANO_FLAGS="base_compiledir=."
		set +u
		source activate gatk
		export HOME=/local/scratch
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		source activate gatk
		touch "${group}_${i}.tar"

		${gatk_call_cnv_version(task)}
		"""
}
def gatk_call_cnv_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
	END_VERSIONS	
	"""
}

process postprocessgatk {
	cpus 5
	memory '50GB'
	time '3h'
	container = "${params.container_gatk}"
	publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'
	tag "$id"

	input:
		set group, id, i, file(tar), file(ploidy), shard_no, shard \
			from postprocessgatk.groupTuple(by: [0,1]).join(ploidy_to_post, by: [0,1]).combine(gatk_postprocess.groupTuple(by: [3]))


	output:
		set group, id, \
			file("genotyped-intervals-${group}-vs-cohort30.vcf.gz"), \
			file("genotyped-segments-${group}-vs-cohort30.vcf.gz"), \
			file("denoised-${group}-vs-cohort30.vcf.gz") into called_gatk
		set group, file("*versions.yml") into ch_postprocessgatk_versions

	script:
		modelshards = shard.join(' --model-shard-path ') // join each reference shard
		caseshards = []
		for (n = 1; n <= i.size(); n++) { // join each shard(n) that's been called
			tmp = group+'_'+i[n-1]+'/'+group+'_'+i[n-1]+'-calls' 
			caseshards = caseshards + tmp
		}
		caseshards = caseshards.join( ' --calls-shard-path ')
 	
	shell:
		version_str = postprocessgatk_version(task)
		'''
		THEANO_FLAGS="base_compiledir=/fs1/resources/theano"
		for model in !{tar}; do
		tar -xvf $model
		done
		tar -xvf !{ploidy}
		set +u
		source activate gatk
		export MKL_NUM_THREADS=!{task.cpus}
		export OMP_NUM_THREADS=!{task.cpus}
		gatk --java-options "-Xmx25g" PostprocessGermlineCNVCalls \
			--allosomal-contig X --allosomal-contig Y \
			--contig-ploidy-calls ploidy/!{group}-calls/ \
			--sample-index 0 \\
			--output-genotyped-intervals genotyped-intervals-!{group}-vs-cohort30.vcf.gz \
			--output-genotyped-segments genotyped-segments-!{group}-vs-cohort30.vcf.gz \
			--output-denoised-copy-ratios denoised-!{group}-vs-cohort30.vcf.gz \
			--sequence-dictionary !{params.GENOMEDICT} \
			--calls-shard-path !{caseshards} \
			--model-shard-path !{modelshards}

		echo "!{version_str}" > "!{task.process}_versions.yml"
		'''

	stub:
		version_str = postprocessgatk_version(task)
		"""
		THEANO_FLAGS="base_compiledir=/fs1/resources/theano"
		set +u
		source activate gatk
		export MKL_NUM_THREADS=!{task.cpus}
		export OMP_NUM_THREADS=!{task.cpus}
		source activate gatk
		touch "genotyped-intervals-${group}-vs-cohort30.vcf.gz"
		touch "genotyped-segments-${group}-vs-cohort30.vcf.gz"
		touch "denoised-${group}-vs-cohort30.vcf.gz"

		echo "${version_str}" > "${task.process}_versions.yml"
		"""
}
def postprocessgatk_version(task) {
	// This docstring looks different from others as it is used within the shell '''
	// If spaces similarly to the others, this leads to additional whitespace above and below the version text
	"""${task.process}:
	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')"""
}

process filter_merge_gatk {
	cpus 2
	tag "$group"
	time '2h'
	memory '1 GB'
	publishDir "${OUTDIR}/sv_vcf", mode: 'copy', overwrite: 'true'

	input:
		set group, id, file(inter), file(gatk), file(denoised) from called_gatk

	output:
		set group, id, file("${id}.gatk.filtered.merged.vcf") into merged_gatk, merged_gatk_panel

	script:
		"""
		filter_gatk.pl $gatk > ${id}.gatk.filtered.vcf
		mergeGATK.pl ${id}.gatk.filtered.vcf > ${id}.gatk.filtered.merged.vcf
		"""

	stub:
		"""
		touch "${id}.gatk.filtered.merged.vcf"
		"""
}

process manta {
	cpus = 56
	publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'
	tag "$id"
	time '10h'
	memory '150 GB'
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'

	when:
		params.sv && params.antype == "wgs"

	input:
		set group, id, file(bam), file(bai) from bam_manta.mix(bam_manta_choice)

	output:
		set group, id, file("${id}.manta.vcf.gz") into called_manta
		set group, file("*versions.yml") into ch_manta_versions

	script:
		bams = bam.join('--bam ')

		"""
		configManta.py --bam $bam --reference $genome_file --runDir .
		python runWorkflow.py -m local -j ${task.cpus}
		mv results/variants/diploidSV.vcf.gz ${id}.manta.vcf.gz
		mv results/variants/diploidSV.vcf.gz.tbi ${id}.manta.vcf.gz.tbi

		${manta_version(task)}
		"""

	stub:
		"""
		touch "${id}.manta.vcf.gz"
		${manta_version(task)}
		"""
}
def manta_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    manta: \$( configManta.py --version )
	END_VERSIONS	
	"""
}

process manta_panel {
	cpus = 20
	publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'
	tag "$id"
	time '1h'
	memory '50 GB'
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'

	when:
		params.sv && params.antype == "panel"

	input:
		set group, id, file(bam), file(bai) from bam_manta_panel.mix(bam_mantapanel_choice)

	output:
		set group, id, file("${id}.manta.vcf.gz") into called_manta_panel
		set group, file("*versions.yml") into ch_manta_panel_versions

	script:
		"""
		configManta.py --bam $bam --reference $genome_file --runDir . --exome --callRegions $params.bedgz --generateEvidenceBam
		python runWorkflow.py -m local -j ${task.cpus}
		mv results/variants/diploidSV.vcf.gz ${id}.manta.vcf.gz
		mv results/variants/diploidSV.vcf.gz.tbi ${id}.manta.vcf.gz.tbi

		${manta_panel_version(task)}
		"""

	stub:
		"""
		touch "${id}.manta.vcf.gz"
		${manta_panel_version(task)}
		"""
}
def manta_panel_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    manta: \$( configManta.py --version )
	END_VERSIONS
	"""
}


process cnvkit_panel {
	cpus = 5
	container = "${params.container_twist_myeloid}"
	publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	publishDir "${OUTDIR}/plots/", mode: 'copy', overwrite: 'true', pattern: '*.png'
	tag "$id"
	time '1h'
	memory '20 GB'
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'

	when:
		params.sv && params.antype == "panel"

	input:
		set group, id, file(bam), file(bai), file(vcf), file(multi), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) from bam_cnvkit_panel.mix(bam_cnvkitpanel_choice).join(vcf_cnvkit, by:[0,1]).join(qc_cnvkit_val, by:[0,1])
		//set id, val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) from qc_cnvkit_val.view()
		//set group, id, file(vcf) from vcf_cnvkit.view()
	
	output:
		set group, id, file("${id}.cnvkit_filtered.vcf") into called_cnvkit_panel
		file("${id}.call.cns") into unfiltered_cns
		file("${group}.genomic_overview.png")
		set group, file("${group}_oplot.INFO") into cnvkit_INFO
		set group, file("*versions.yml") into ch_cnvkit_panel_versions

	script:
		"""
		cnvkit.py batch $bam -r $params.cnvkit_reference -p 5 -d results/
		cnvkit.py call results/*.cns -v $vcf -o ${id}.call.cns
		filter_cnvkit.pl ${id}.call.cns $MEAN_DEPTH > ${id}.filtered
		cnvkit.py export vcf ${id}.filtered -i "$id" > ${id}.cnvkit_filtered.vcf
		cnvkit.py scatter -s results/*dedup.cn{s,r} -o ${group}.genomic_overview.png -v $vcf -i $id
		echo "IMG overviewplot	${params.accessdir}/plots/${group}.genomic_overview.png" > ${group}_oplot.INFO
		
		${cnvkit_panel_version(task)}
		"""

	stub:
		"""
		touch "${id}.cnvkit_filtered.vcf"
		touch "${id}.call.cns"
		touch "${group}.genomic_overview.png"
		touch "${group}_oplot.INFO"

		${cnvkit_panel_version(task)}
		"""
}
def cnvkit_panel_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
	END_VERSIONS	
	"""
}

process svdb_merge_panel {
	container = "${params.container_svdb}"
	cpus 2
	cache 'deep'
	tag "$group"
	publishDir "${OUTDIR}/sv_vcf/merged/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	time '1h'
	memory '1 GB'

	when:
		params.antype == "panel"

	input:
		set group, id, file(vcfs) from called_manta_panel.mix(called_cnvkit_panel, merged_gatk_panel).groupTuple()
				
	output:
		set group, id, file("${group}.merged.vcf") into ch_postprocess_merged_panel_sv
		set group, file("*versions.yml") into ch_svdb_merge_panel_versions

	script:
		if (vcfs.size() > 1) {
			// for each sv-caller add idx, find vcf and find priority, add in priority order! //
			// index of vcfs added from mix //
			manta_idx = vcfs.findIndexOf{ it =~ 'manta' }
			cnvkit_idx = vcfs.findIndexOf{ it =~ 'cnvkit' }
			gatk_idx = vcfs.findIndexOf{ it =~ 'gatk' }

			// find vcfs //
			manta = manta_idx >= 0 ? vcfs[manta_idx].collect {it + ':manta ' } : null
			cnvkit = cnvkit_idx >= 0 ? vcfs[cnvkit_idx].collect {it + ':cnvkit ' } : null
			gatk = gatk_idx >= 0 ? vcfs[gatk_idx].collect {it + ':gatk ' } : null
			tmp = manta + gatk + cnvkit
			tmp = tmp - null
			vcfs_svdb = tmp.join(' ')

			// find priorities //
			mantap = manta_idx >= 0 ? 'manta' : null
			gatkp = gatk_idx >= 0 ? 'gatk' : null
			cnvkitp = cnvkit_idx >= 0 ? 'cnvkit' : null
			tmpp = [mantap, gatkp, cnvkitp]
			tmpp = tmpp - null
			priority = tmpp.join(',')
			
			"""
			svdb \\
			  --merge \\
			  --vcf $vcfs_svdb \\
			  --no_intra \\
			  --pass_only \\
			  --bnd_distance 2500 \\
			  --overlap 0.7 \\
			  --priority $priority \\
			  --ins_distance 0 > ${group}.merged.tmp


			# copy callers out of INFO.set to INFO.SCOUT_CUSTOM
			add_callers_to_scout_custom.py \\
				--callers $priority \\
				--merged_vcf ${group}.merged.tmp > ${group}.merged.callers.tmp

			add_vcf_header_info_records.py \\
				--vcf ${group}.merged.callers.tmp \\
				--info SCOUT_CUSTOM . String "Custom annotations for scout" '' '' \\
				--output ${group}.merged.vcf

			${svdb_merge_panel_version(task)}
			"""
		}
		else {
			"""
			mv $vcf ${group}.merged.vcf
			${svdb_merge_panel_version(task)}
			"""
		}
	
	stub:
		"""
		touch "${group}.merged.vcf"
		${svdb_merge_panel_version(task)}
		"""
}
def svdb_merge_panel_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
	END_VERSIONS
	"""
}

process postprocess_merged_panel_sv_vcf {
	cpus 2
	tag "$group"
	publishDir "${OUTDIR}/sv_vcf/merged/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	time '1h'
	memory '1 GB'


	input:
		set group, id, file(merged_vcf) from ch_postprocess_merged_panel_sv
		set group, id, file(melt_vcf) from ch_melt_vcf

	output:
		set group, id, file("${group}.merged.bndless.genotypefix.melt.vcf") into vep_sv_panel, annotsv_panel
		set group, file("${group}.merged.bndless.genotypefix.melt.vcf") into loqusdb_sv_panel
		set group, file("*versions.yml") into ch_postprocess_merged_panel_sv_vcf_versions


	script:
		"""
		# Remove BNDs
		grep -v "BND" $merged_vcf > ${group}.merged.bndless.vcf

		# Any 0/0 GT -> 0/1, otherwise loqus will reject them.
		modify_cnv_genotypes_for_loqusdb.pl --merged_panel_sv_vcf ${group}.merged.bndless.vcf > ${group}.merged.bndless.genotypefix.vcf

		# Add MELT data to info vars:
		add_vcf_header_info_records.py \\
			--vcf ${group}.merged.bndless.genotypefix.vcf \\
			--info MELT_RANK . String "Evidence level 1-5, 5 - highest" '' '' \\
			--info MELT_QC . String "Quality of call" '' '' \\
			--output ${group}.merged.bndless.genotypefix.headers.vcf

		# Combine with MELT:
		vcf-concat  ${group}.merged.bndless.genotypefix.headers.vcf $melt_vcf | vcf-sort -c > ${group}.merged.bndless.genotypefix.melt.vcf
		${postprocess_merged_panel_sv_version(task)}
		"""

	stub:
		"""
		touch ${group}.merged.bndless.genotypefix.melt.vcf
		${postprocess_merged_panel_sv_version(task)}
		"""

}
def postprocess_merged_panel_sv_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcftools: \$(vcftools --version | cut -f 2 -d " " | tr -d "()")
	END_VERSIONS
	"""
}

process tiddit {
	cpus = 2
	publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	time '10h'
	tag "$id"
	memory '15 GB'
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'

	when:
		params.sv && params.antype == "wgs"


	input:
		set group, id, file(bam), file(bai) from bam_tiddit.mix(bam_tiddit_choice)

	output:
		set group, id, file("${id}.tiddit.filtered.vcf") into called_tiddit
		set group, file("*versions.yml") into ch_tiddit_versions

	script:
		"""
		TIDDIT.py --sv -o ${id}.tiddit --bam $bam
		grep -E \"#|PASS\" ${id}.tiddit.vcf > ${id}.tiddit.filtered.vcf
		${tiddit_version(task)}
		"""

	stub:
		"""
		touch "${id}.tiddit.filtered.vcf"
		${tiddit_version(task)}
		"""
}
def tiddit_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    tiddit: \$(echo \$(TIDDIT.py 2>&1) | sed 's/^.*TIDDIT-//; s/ .*\$//')
	END_VERSIONS	
	"""
}

process svdb_merge {
	cpus 2
	container = "${params.container_svdb}"
	tag "$group"
	publishDir "${OUTDIR}/sv_vcf/merged/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	time '2h'
	memory '1 GB'

	input:
		set group, id, file(mantaV) from called_manta.groupTuple()
		set group, id, file(tidditV) from called_tiddit.groupTuple()
		set group, id, file(gatkV) from merged_gatk.groupTuple()
		
	output:
		set group, id, file("${group}.merged.bndless.vcf") into vep_sv, annotsv_vcf
		set group, file("${group}.merged.vcf") into loqusdb_sv
		set group, file("*versions.yml") into ch_svdb_merge_versions

	script:
		if (mode == "family") {
			vcfs = []
			manta = []
			tiddit = []
			gatk = []

			/*
			 Order in which VCFs are merged matters when the merged SV
			 is annotated with final position/length, which affects
			 artefact matching in loqusdb.

			 A possibly better way to sort here would be to sort the
			 file by familial-relation (e.g. always sort proband-mother-father)
			 this would ensure the same merge-order regardless of sample-id
			 */

			mantaV = mantaV.collect { it.toString() }.sort()
			gatkV = gatkV.collect { it.toString() }.sort()
			tidditV = tidditV.collect { it.toString() }.sort()

			for (i = 1; i <= mantaV.size(); i++) {
				tmp = mantaV[i-1] + ':manta' + "${i}"
				tmp1 = tidditV[i-1] + ':tiddit' + "${i}"
				tmp2 = gatkV[i-1] + ':gatk' + "${i}"
				vcfs = vcfs + tmp + tmp1 + tmp2
				mt = 'manta' + "${i}"
				tt = 'tiddit' + "${i}"
				ct = 'gatk' + "${i}"
				manta = manta + mt
				tiddit = tiddit + tt
				gatk = gatk + ct
			}

			prio = manta + tiddit + gatk
			prio = prio.join(',')
			vcfs = vcfs.join(' ')
			"""
			svdb \\
				--merge \\
				--vcf $vcfs \\
				--no_intra \\
				--pass_only \\
				--bnd_distance 2500 \\
				--overlap 0.7 \\
				--priority $prio \\
				--ins_distance 0 > ${group}.merged.vcf

			grep -v BND ${group}.merged.vcf > ${group}.merged.bndless.vcf
			
			${svdb_merge_version(task)}
			"""
		}

		else {
			tmp = mantaV.collect {it + ':manta ' } + tidditV.collect {it + ':tiddit ' } + gatkV.collect {it + ':gatk ' }
			vcfs = tmp.join(' ')
			"""
			svdb \\
				--merge \\
				--vcf $vcfs \\
				--no_intra \\
				--pass_only \\
				--bnd_distance 2500 \\
				--overlap 0.7 \\
				--priority manta,tiddit,gatk \\
				--ins_distance 0 > ${group}.merged.vcf

			grep -v BND ${group}.merged.vcf > ${group}.merged.bndless.vcf
			
			${svdb_merge_version(task)}
			"""
		}

	stub:
		"""
		touch "${group}.merged.vcf"
		touch "${group}.merged.bndless.vcf"

		${svdb_merge_version(task)}
		"""
}
def svdb_merge_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
	END_VERSIONS
	"""
}

process dummy_svvcf_for_loqusdb {

	// add_to_loqusb won't run if no svvcf is generated
	// this process creates dummy svvcf for no-SV runs
	// assay input only exists to disable nextflow warning
	// for channels emitting with less than two input elements

	cpus 1
	tag "$group"
	memory '10 MB'
	time '10m'

	when:
		!params.sv
	input:
		set group, assay from meta_loqusdb_no_sv_calling
	output:
		set group, file("${group}.dummy.sv.vcf") into dummy_svvcf_ch

	script:
		"""
		touch ${group}.dummy.sv.vcf
		"""

	stub:
		"""
		touch ${group}.dummy.sv.vcf
		"""
}

process add_to_loqusdb {
	cpus 1
	publishDir "${CRONDIR}/loqus", mode: 'copy' , overwrite: 'true'
	tag "$group"
	memory '100 MB'
	time '25m'

	when:
		!params.noupload && !params.reanalyze

	input:
		set group, type, file(vcf), file(tbi), type, file(ped) from vcf_loqus.join(ped_loqus)
		set group, file(svvcf) from loqusdb_sv.mix(loqusdb_sv_panel, dummy_svvcf_ch)

	output:
		file("${group}*.loqus") into loqusdb_done

	script:
		"""
		sv_variants=""
		nbr_svvcf_records=\$(grep -v '^#' ${svvcf} | wc -l)

		if (( \$nbr_svvcf_records > 0 )); then
			sv_variants="--sv-variants ${params.accessdir}/sv_vcf/merged/${svvcf}"
		fi

		echo "-db $params.loqusdb load -f ${params.accessdir}/ped/${ped} --variant-file ${params.accessdir}/vcf/${vcf} \$sv_variants" > ${group}.loqus
		"""

	stub:
		"""
		touch "${group}.loqus"
		"""
}

process annotsv {
	container = "${params.container_annotsv}"
	cpus 2
	tag "$group"
	publishDir "${OUTDIR}/annotsv/", mode: 'copy', overwrite: 'true', pattern: '*.tsv'
	time '5h'
	memory '20 GB'

	input:
		set group, id, file(sv) from annotsv_vcf.mix(annotsv_panel)
			
	output:
		set group, file("${group}_annotsv.tsv") into annotsv, annotsv_ma, annotsv_fa
		set group, file("*versions.yml") into ch_annotsv_versions

	shell:
		version_str = annotsv_version(task)
		'''
		export ANNOTSV="/AnnotSV"
		/AnnotSV/bin/AnnotSV -SvinputFile !{sv} \\
			-typeOfAnnotation full \\
			-outputDir !{group} \\
			-genomeBuild GRCh38
		if [ -f !{group}/*.annotated.tsv ]; then
			mv !{group}/*.annotated.tsv !{group}_annotsv.tsv
		else
		    echo "1\n" > !{group}_annotsv.tsv
		fi
		echo "!{version_str}" > "!{task.process}_versions.yml"
		'''

	stub:
		version_str = annotsv_version(task)
		"""
		export ANNOTSV="/AnnotSV"
		touch "${group}_annotsv.tsv"

		echo "${version_str}" > "${task.process}_versions.yml"
		"""
}
def annotsv_version(task) {
	"""${task.process}:
	    annotsv: \$( echo \$(/AnnotSV/bin/AnnotSV --version) | sed -e "s/AnnotSV //g ; s/Copyright.*//" )"""
}

process vep_sv {
	cpus 10
	container = "${params.container_vep}"
	tag "$group"
	memory '50 GB'
	time '1h'
	
	input:
		set group, id, file(vcf) from vep_sv.mix(vep_sv_panel)

	output:
		set group, id, file("${group}.vep.vcf") into vep_sv_vcf
		set group, file("*versions.yml") into ch_vep_sv_versions

	script:
		"""

		# Temporary fix for VEP 111.0 annotation bug, where certain MANTA indels are being skipped by VEP
		# See: https://github.com/Ensembl/ensembl-vep/issues/1631#issuecomment-1985973568
		# Edit 2024-11-01: Fixed in VEP 112.0

		sed 's/SVTYPE=/BAZBAZ=/' $vcf > ${group}.vep111-workaround.vcf

		vep \\
			-i ${group}.vep111-workaround.vcf \\
			-o ${group}.vep.vcf \\
			--offline \\
			--merged \\
			--sift b --polyphen b --ccds --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --af --af_1kg --max_af --pubmed --uniprot --mane --tsl --appris --variant_class --gene_phenotype --mirna \\
			--synonyms $params.VEP_SYNONYMS \\
			--vcf \\
			--no_stats \\
			--fork ${task.cpus} \\
			--force_overwrite \\
			--plugin LoFtool \\
			--fasta $params.VEP_FASTA \\
			--dir_cache $params.VEP_CACHE \\
			--dir_plugins $params.VEP_PLUGINS \\
			--max_sv_size $params.VEP_MAX_SV_SIZE \\
			--distance $params.VEP_TRANSCRIPT_DISTANCE \\
			-cache \\
			--format vcf

		# Re-enable SVTYPE:
		sed -i 's/BAZBAZ=/SVTYPE=/' ${group}.vep.vcf
		${vep_sv_version(task)}
		"""

	stub:
		"""
		touch "${group}.vep.vcf"
		${vep_sv_version(task)}
		"""
}
def vep_sv_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
	END_VERSIONS
	"""
}

process postprocess_vep_sv {
	cpus = 2
	memory '10GB'
	time '1h'
	tag "$group"
	container = "${params.container_svdb}"

	input:
		set group, id, file(vcf) from vep_sv_vcf

	output:
		set group, file("${group}.vep.clean.merge.vcf") into add_omim_vcf
		set group, file("*versions.yml") into ch_postprocess_vep_versions
	
	script:
		"""
		# Filter variants with FILTER != . or PASS and variants missing CSQ field.
		postprocess_vep_vcf.py $vcf > ${group}.vep.clean.vcf
		svdb --merge --overlap 0.9 --notag --vcf ${group}.vep.clean.vcf --ins_distance 0 > ${group}.vep.clean.merge.tmp.vcf

		# --notag above will remove set from INFO:
		add_vcf_header_info_records.py \\
			--vcf ${group}.vep.clean.merge.tmp.vcf \\
			--info set 1 String "Source VCF for the merged record in SVDB" '' '' \\
			--info VARID 1 String "The variant ID of merged samples" '' '' \\
			--output ${group}.vep.clean.merge.headers.tmp.vcf

		# Prepare annotations for scout:
		normalize_caller_names_in_svdb_fields.py ${group}.vep.clean.merge.headers.tmp.vcf --callers manta gatk tiddit > ${group}.vep.clean.merge.vcf
		${postprocess_vep_sv_version(task)}
		"""
	stub:
		"""
		touch "${group}.vep.clean.merge.vcf"
		${postprocess_vep_sv_version(task)}
		"""
}
def postprocess_vep_sv_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
	END_VERSIONS
	"""
}

process add_omim {
	cpus = 2
	memory '10GB'
	time '1h'
	tag "$group"

	input:
		set group, file(vcf) from add_omim_vcf

	output:
		set group, file("${group}.vep.clean.merge.omim.vcf") into artefact_vcf

	script:
		"""
		add_omim.pl $params.OMIM_GENES < $vcf > ${group}.vep.clean.merge.omim.vcf
		"""
	stub:
		"""
		touch "${group}.vep.clean.merge.omim.vcf"
		"""
}



// Query artefact db
process artefact {
	cpus 2
	tag "$group"
	time '10h'
	memory '10 GB'
	container = "${params.container_svdb}"
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'

	input:
		set group, file(sv) from artefact_vcf

	output:
		set group, file("${group}.artefact.vcf") into manip_vcf,manip_vcf_ma,manip_vcf_fa
		set group, file("*versions.yml") into ch_artefact_versions

	script:
		// use loqusdb dump not svdb database //
		if (params.gatkcnv) {
			"""
			svdb \\
			--query --bnd_distance 25000 --overlap 0.7 --in_occ Obs --out_occ ACOUNT --in_frq Frq --out_frq AFRQ  \\
			--db $params.svdb \\
			--ins_distance 0 \\
			--query_vcf $sv > ${group}.artefact.vcf

			${artefact_version(task)}
			"""
		}
		// for oncov1-0 still use svdb database remove in future//
		else {
			"""
			svdb \\
			--sqdb $params.svdb \\
			--query \\
			--query_vcf $sv \\
			--out_occ ACOUNT \\
			--ins_distance 0 \\
			--out_frq AFRQ > ${group}.artefact.vcf

			${artefact_version(task)}
			"""
		}

	stub:
		"""
		touch "${group}.artefact.vcf"
		${artefact_version(task)}
		"""
}
def artefact_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
	END_VERSIONS
	"""
}


process prescore {
	cpus 2
	tag "$group"
	memory '10 GB'
	time '1h'

	input:
		set group, file(sv_artefact), type, file(ped), file(annotsv) from manip_vcf.mix(manip_vcf_ma,manip_vcf_fa).join(ped_prescore.mix(ped_prescore_ma,ped_prescore_fa)).join(annotsv.mix(annotsv_ma,annotsv_fa))

	output:
		set group, type, file("${group}.annotatedSV.vcf") into annotatedSV

	script:
		"""
		prescore_sv.pl \\
		--sv $sv_artefact --ped $ped --annotsv $annotsv --osv ${group}.annotatedSV.vcf
		"""

	stub:
		"""
		touch "${group}.annotatedSV.vcf"
		"""
}

process score_sv {
	tag "$group $mode"
	cpus 2
	memory '10 GB'
	time '2h'
	container = "${params.container_genmod}"

	input:
		set group, type, file(in_vcf) from annotatedSV

	output:
		set group, type, file("${group_score}.sv.scored.vcf") into ch_scored_sv
		set group, file("*versions.yml") into ch_score_sv_versions

	script:
		def model = (mode == "family" && params.antype == "wgs") ? params.svrank_model : params.svrank_model_s
		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group
		"""
		genmod score --family_id ${group_score} --score_config ${model} --rank_results --outfile "${group_score}.sv.scored.vcf" ${in_vcf}

		${score_sv_version(task)}
		"""

	stub:
		group_score = group
		"""
		touch "${group_score}.sv.scored.vcf"

		${score_sv_version(task)}
		"""
}
def score_sv_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    genmod: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
	END_VERSIONS	
	"""
}

process bgzip_scored_genmod {
	tag "$group"
	cpus 4
	memory '1 GB'
	time '5m'
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz.tbi'
	container = "${params.container_bcftools}"

	input:
		set group, type, file(scored_sv_vcf) from ch_scored_sv
	
	output:
		set group, type, file("${group_score}.sv.scored.sorted.vcf.gz"), file("${group_score}.sv.scored.sorted.vcf.gz.tbi") into sv_rescore, sv_rescore_ma, sv_rescore_fa
		set group, file("${group_score}.sv.scored.sorted.vcf.gz") into svvcf_bed, svvcf_pod
		set group, file("${group}_sv.INFO") into sv_INFO
		set group, file("*versions.yml") into ch_bgzip_scored_genmod_versions

	script:
		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group
		"""
			bcftools sort -O v -o ${group_score}.sv.scored.sorted.vcf ${scored_sv_vcf}
			bgzip -@ ${task.cpus} ${group_score}.sv.scored.sorted.vcf -f
			tabix ${group_score}.sv.scored.sorted.vcf.gz -f
			echo "SV\t$type\t${params.accessdir}/vcf/${group_score}.sv.scored.sorted.vcf.gz" > ${group}_sv.INFO

			${bgzip_score_sv_version(task)}
		"""
	stub:
		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group
		"""
			touch "${group_score}.sv.scored.sorted.vcf.gz"
			touch "${group_score}.sv.scored.sorted.vcf.gz.tbi"
			touch "${group}_sv.INFO"

			${bgzip_score_sv_version(task)}
		"""
}
def bgzip_score_sv_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS	
	"""
}

process compound_finder {
	cpus 2
	tag "$group $mode"
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
	memory '10 GB'
	time '2h'

	when:
		mode == "family" && params.assay == "wgs"

	input:
		set group, type, file(vcf), file(tbi), file(ped), file(snv), file(tbi) from sv_rescore.mix(sv_rescore_ma,sv_rescore_fa).join(ped_compound.mix(ped_compound_ma,ped_compound_fa), by: [0,1]).join(snv_sv_vcf.mix(snv_sv_vcf_ma,snv_sv_vcf_fa),by: [0,1])
		//set group, file(snv), file(tbi) from snv_sv_vcf

	output:
		set group, file("${group_score}.snv.rescored.sorted.vcf.gz"), file("${group_score}.snv.rescored.sorted.vcf.gz.tbi") into vcf_yaml
		set group, file("${group}_svp.INFO") into svcompound_INFO
		set group, file("*versions.yml") into ch_compound_finder_versions

	script:
		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group

		"""
		compound_finder.pl \\
			--sv $vcf --ped $ped --snv $snv \\
			--osv ${group_score}.sv.rescored.sorted.vcf \\
			--osnv ${group_score}.snv.rescored.sorted.vcf \\
			--skipsv
		bgzip -@ ${task.cpus} ${group_score}.snv.rescored.sorted.vcf -f
		tabix ${group_score}.snv.rescored.sorted.vcf.gz -f
		echo "SVc	$type	${params.accessdir}/vcf/${group_score}.sv.scored.sorted.vcf.gz,${params.accessdir}/vcf/${group_score}.snv.rescored.sorted.vcf.gz" > ${group}_svp.INFO
		
		${compound_finder_version(task)}
		"""

	stub:
		group_score = group
		"""
		touch "${group_score}.snv.rescored.sorted.vcf.gz"
		touch "${group_score}.snv.rescored.sorted.vcf.gz.tbi"
		touch "${group}_svp.INFO"

		${compound_finder_version(task)}
		"""
}
def compound_finder_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	END_VERSIONS	
	"""
}


process output_files {
	cpus 2
	memory '1GB'
	time '1h'

	input:
		set group, files from bam_INFO.mix(snv_INFO,sv_INFO,str_INFO,peddy_INFO,madde_INFO,svcompound_INFO,smn_INFO,bamchoice_INFO,mtBAM_INFO,oplot_INFO,haplogrep_INFO,eklipse_INFO,cnvkit_INFO,d4_INFO).groupTuple()

	output:
		set group, file("${group}.INFO") into yaml_INFO

	script:
		files = files.join( ' ' )

		"""
		cat $files > ${group}.INFO
		"""

	stub:
		"""
		touch "${group}.INFO"
		"""
}


process svvcf_to_bed {
	publishDir "${OUTDIR}/bed", mode: 'copy' , overwrite: 'true'
	tag "group"
	memory '1 GB'
	time '1h'
	cpus 2

	when:
		params.antype != "panel"

	input:
		set group, file(vcf) from svvcf_bed
		set group, id, sex, type from meta_svbed.filter { item -> item[3] == 'proband' }

	output:
		file("${group}.sv.bed")


	script:
		"""
		cnv2bed.pl --cnv $vcf --pb $id > ${group}.sv.bed
		"""

	stub:
		"""
		touch "${group}.sv.bed"
		"""
}

process plot_pod {
	container = "${params.container_pod}"
	publishDir "${OUTDIR}/pod", mode: 'copy' , overwrite: 'true'
	tag "$group"
	time '1h'
	memory '1 GB'
	cpus 2

	input:
		set group, file(snv) from vcf_pod
		set group, file(cnv), type, file(ped) from svvcf_pod.join(ped_pod)
		set group, id, sex, type from meta_pod.filter { item -> item[3] == 'proband' }		

	output:
		set file("${id}_POD_karyotype.pdf"), file("${id}_POD_results.html")

	when:
		mode == "family" && trio == true

	script:
		"""
		parental_origin_of_duplication.pl --snv $snv --cnv $cnv --proband $id --ped $ped
		"""

	stub:
		"""
		touch "${id}_POD_karyotype.pdf"
		touch "${id}_POD_results.html"
		"""
}

process create_yaml {
	publishDir "${OUTDIR}/yaml", mode: 'copy' , overwrite: 'true', pattern: '*.yaml'
	publishDir "${OUTDIR}/yaml/alt_affect", mode: 'copy' , overwrite: 'true', pattern: '*.yaml.*a'
	publishDir "${CRONDIR}/scout", mode: 'copy' , overwrite: 'true', pattern: '*.yaml'
	errorStrategy 'retry'
	maxErrors 5
	tag "$group"
	time '5m'
	memory '1 GB'

	input:
		set group, id, sex, mother, father, phenotype, diagnosis, type, assay, clarity_sample_id, ffpe, analysis, type, file(ped), file(INFO) from yml_diag.join(ped_scout).join(yaml_INFO)

	output:
		set group, file("${group}.yaml*") into yaml

	script:
		def assay = params.dev ? "dev,${analysis}" : "${assay},${analysis}"
		print(">>> params.dev: ${params.dev}")
		print(">>> Assay: ${assay}")
		print(">>> Analysis: ${analysis}")
		"""
		create_yml.pl \\
			--g "${group},${clarity_sample_id}" \\
			--d "$diagnosis" \\
			--panelsdef "$params.panelsdef" \\
			--out "${group}.yaml" \\
			--ped "$ped" \\
			--files "$INFO" \\
			--assay "$assay" \\
			--antype "$params.antype" \\
			--extra_panels "$params.extra_panels"
		"""

	stub:
		"""
		touch "${group}.yaml"
		"""
}

process combine_versions {
	publishDir "${OUTDIR}/versions", mode: 'copy', overwrite: 'true', pattern: '*.versions.yml'

	// The point of "first" here is that when a process is present in multiple instances
	// there is no need to include more than one instance of the versions
	input:
		set group, versions from ch_fastp_versions.first().mix(
			ch_bwa_align_shareded_versions.first(),
			ch_bwa_merge_shards_versions.first(),
			ch_bwa_align_versions.first(), 
			ch_markdup_versions.first(),
			ch_bqsr_versions.first(),
			ch_sentieon_qc_versions.first(),
			ch_d4_coverage_versions.first(),
			ch_verifybamid2_versions.first(),
			ch_smn_copy_number_caller_versions.first(),
			ch_expansionhunter_versions.first(),
			ch_stranger_versions.first(),
			ch_reviewer_versions.first(),
			ch_vcfbreakmulti_expansionhunter_versions.first(),
			ch_melt_versions.first(),
			ch_intersect_melt_versions.first(),
			ch_dnascope_versions.first(),
			ch_gvcf_combine_versions.first(),
			ch_madeline_versions.first(),
			ch_freebayes_versions.first(),
			ch_fetch_mt_seqs_versions.first(),
			ch_sentieon_mitochondrial_qc.first(),
			ch_run_mutect2_versions.first(),
			ch_split_normalize_mito_versions.first(),
			ch_run_hmtnote_versions.first(),
			ch_run_haplogrep_versions.first(),
			ch_run_eklipse_versions.first(),
			ch_split_normalize_versions.first(),
			ch_annotate_vep_versions.first(),
			ch_vcfanno_versions.first(),
			ch_extract_indels_for_cadd_versions.first(),
			ch_indel_vep_versions.first(),
			ch_calculate_indel_cadd_versions.first(),
			ch_add_cadd_scores_to_vcf_versions.first(),
			ch_inher_models_versions.first(),
			ch_genmodscore_versions.first(),
			ch_vcf_completion_versions.first(),
			ch_peddy_versions.first(),
			ch_upd_versions.first(),
			ch_roh_versions.first(),
			ch_gatkcov_versions.first(),
			ch_gatk_coverage_versions.first(),
			ch_gatk_call_ploidy_versions.first(),
			ch_gatk_call_cnv_versions.first(),
			ch_postprocessgatk_versions.first(),
			ch_manta_versions.first(),
			ch_manta_panel_versions.first(),
			ch_cnvkit_panel_versions.first(),
			ch_postprocess_merged_panel_sv_vcf_versions.first(),
			ch_svdb_merge_panel_versions.first(),
			ch_tiddit_versions.first(),
			ch_svdb_merge_versions.first(),
			ch_annotsv_versions.first(),
			ch_vep_sv_versions.first(),
			ch_postprocess_vep_versions.first(),
			ch_artefact_versions.first(),
 			ch_score_sv_versions.first(),
			ch_compound_finder_versions.first(),
			ch_bgzip_indel_cadd_versions.first(),
			ch_bgzip_scored_genmod_versions.first()
		).groupTuple()
	
	output:
		file("${group}.versions.yml")
	
	script:
		versions_joined = versions.sort( my_it -> my_it.name ).join(" ")
		"""
		cat $versions_joined > ${group}.versions.yml
		"""
	
	stub:
		versions_joined = versions.sort( my_it -> my_it.name ).join(" ")
		"""
		cat $versions_joined > ${group}.versions.yml
		"""
}
