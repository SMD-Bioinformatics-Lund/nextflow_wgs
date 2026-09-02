workflow MITOCHONDRIAL_ANALYSIS {

	take:
	ch_bam_bai
	ch_sample_meta
	ch_proband_meta
	ch_snv_vcf_tbi_intersected
	val_run_mito_qc
	val_run_mito_mutect2
	val_genome_fasta
	val_rcrs_fasta
	val_results_output_dir
	val_accessdir
	val_mito_bam_accessdir

	main:
	ch_versions    = channel.empty()
	ch_output_info = channel.empty()
	ch_qc_json     = channel.empty()
	ch_merged_vcf = channel.empty()

	// TODO: Split the MT BAM extraction/QC, mitochondrial SNV calling,
    // and eKLIPse deletion calling into smaller subworkflows.
	fetch_MTseqs(ch_bam_bai, val_results_output_dir, val_mito_bam_accessdir)
	ch_output_info = ch_output_info.mix(fetch_MTseqs.out.mtBAM_INFO)
	ch_versions = ch_versions.mix(fetch_MTseqs.out.versions.first())

	if (val_run_mito_qc) {
		sentieon_mitochondrial_qc(fetch_MTseqs.out.bam_bai, val_genome_fasta)
		build_mitochondrial_qc_json(sentieon_mitochondrial_qc.out.qc_tsv)

		ch_qc_json = ch_qc_json.mix(build_mitochondrial_qc_json.out.qc_json)
		ch_versions = ch_versions.mix(sentieon_mitochondrial_qc.out.versions.first())
	}

	if (val_run_mito_mutect2) {
		// TODO: Decide whether disabled mitochondrial Mutect2 should fall back to
        // ch_snv_vcf_tbi_intersected instead of emitting an empty merged_vcf channel.
		ch_mutect2_input = fetch_MTseqs.out.bam_bai.groupTuple()
		run_mutect2(ch_mutect2_input, val_genome_fasta, val_results_output_dir)

		ch_split_normalize_mito_in = run_mutect2.out.vcf
			.join(ch_proband_meta, by:[0])
			.map { group, _id, mito_snv_vcf, proband_id, meta ->
				tuple(group, proband_id, meta, mito_snv_vcf)
			}
		split_normalize_mito(ch_split_normalize_mito_in, val_rcrs_fasta)

		run_hmtnote(split_normalize_mito.out.vcf)

		ch_picard_mergevcfs_in = ch_snv_vcf_tbi_intersected
			.join(run_hmtnote.out.vcf)
			.map { group, snv_vcf, _tbi, mito_vcf -> [ group, snv_vcf, mito_vcf ] }

		picard_mergevcfs(ch_picard_mergevcfs_in)
		ch_merged_vcf = picard_mergevcfs.out.merged_vcf

		run_haplogrep(run_mutect2.out.vcf, val_results_output_dir, val_accessdir)
		ch_output_info = ch_output_info.mix(run_haplogrep.out.haplogrep_INFO)

		ch_versions = ch_versions.mix(run_hmtnote.out.versions.first())
		ch_versions = ch_versions.mix(split_normalize_mito.out.versions.first())
		ch_versions = ch_versions.mix(run_mutect2.out.versions.first())
		ch_versions = ch_versions.mix(picard_mergevcfs.out.versions.first())
		ch_versions = ch_versions.mix(run_haplogrep.out.versions.first())
	}

	run_eklipse(
		ch_sample_meta.join(fetch_MTseqs.out.bam_bai, by : [0, 1]),
		val_results_output_dir,
		val_accessdir
	)
	ch_output_info = ch_output_info.mix(run_eklipse.out.eklipse_INFO)
	ch_versions = ch_versions.mix(run_eklipse.out.versions.first())

	emit:
	merged_vcf = ch_merged_vcf
	qc_json = ch_qc_json
	output_info = ch_output_info
	versions = ch_versions
}

process picard_mergevcfs {

	tag "$group"
	cpus 2
	container "${params.container_picard}"
	time "1h"

	input:
	tuple val(group), path(snv_vcf), path(mito_snv_vcf)

	output:
	tuple val(group), path("${group}.merged.vcf.gz"), path("${group}.merged.vcf.gz.tbi"), emit: merged_vcf
	path("*versions.yml"), emit: versions

	script:
	"""
	picard MergeVcfs \
		--CREATE_INDEX \
		-I ${snv_vcf} \
		-I ${mito_snv_vcf} \
		-O ${group}.merged.vcf.gz

	${picard_mergevcfs_version(task)}
	"""

	stub:
	"""
	touch ${group}.merged.vcf.gz
	touch ${group}.merged.vcf.gz.tbi

	${picard_mergevcfs_version(task)}
	"""
}
def picard_mergevcfs_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    picard: \$( echo \$(picard MergeVcfs --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
	END_VERSIONS
	"""
}

// create an MT BAM file
process fetch_MTseqs {
	cpus 2
	memory '10GB'
	time '1h'
	tag "$id"
	// TODO: Add an explicit container for sambamba/samtools.
	publishDir "${val_results_output_dir}/bam", mode: 'copy', overwrite: true, pattern: '*.bam*'

	input:
		tuple val(group), val(id), path(bam), path(bai)
		val val_results_output_dir
		val val_mito_bam_accessdir

	output:
		tuple val(group), val(id), file ("${id}_mito.bam"), path("${id}_mito.bam.bai"), emit: bam_bai
		tuple val(group), path("${group}_mtbam.INFO"), emit: mtBAM_INFO
		path "*versions.yml", emit: versions

	script:
		"""
		sambamba view -f bam $bam M > ${id}_mito.bam
		samtools index -b ${id}_mito.bam
		echo "mtBAM	$id	${val_mito_bam_accessdir}/${id}_mito.bam" > ${group}_mtbam.INFO

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
	container  "${params.container_sentieon}"

	input:
		tuple val(group), val(id), path(bam), path(bai)
		val val_genome_fasta

	output:
		tuple val(group), val(id), path("${id}_mito_coverage.tsv"), emit: qc_tsv
		path "*versions.yml", emit: versions

	script:
		"""
		sentieon driver \\
			-r ${val_genome_fasta} \\
			-t ${task.cpus} \\
			-i $bam \\
			--algo CoverageMetrics \\
			--omit_base_output  \\
			--omit_locus_stat \\
			--omit_sample_stat \\
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
	// TODO: Add an explicit container for mito_tsv_to_json.py.

	input:
		tuple val(group), val(id), path(mito_qc_file)
	output:
		tuple val(group), val(id), path("${id}_mito_qc.json"), emit: qc_json

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
	// TODO: Add an explicit GATK container and remove the conda activate dependency.
	publishDir "${val_results_output_dir}/vcf", mode: 'copy', overwrite: true, pattern: '*.vcf'

	input:
		tuple val(group), val(id), path(bam), path(bai)
		val val_genome_fasta
		val val_results_output_dir

	output:
		tuple val(group), val(id), path("${group}.mutect2.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		bams = bam.join(' -I ')

		"""
		source activate gatk4-env
		gatk Mutect2 \
		--mitochondria-mode \
		-R ${val_genome_fasta} \
		-L M \
		-I $bams \
		-O ${group}.mutect2.vcf

		${run_mutect2_version(task)}
		"""

	stub:
		bams = bam.join(' -I ')
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
	// TODO: Add an explicit container for bcftools/tabix/filter_mutect2_mito.pl.

	input:
		tuple val(group), val(id), val(meta), path(mito_snv_vcf)
		val val_rcrs_fasta

	output:
		tuple val(group), path("${group}.mutect2.breakmulti.filtered5p.0genotyped.proband.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		# Old workaround to remove false-positive that crashes bcftools norm:
		# TODO: deal w/ this in some better way
		grep -vP "^M\\s+955" ${mito_snv_vcf} > ${mito_snv_vcf}.fix

		bcftools norm -m-both -o ${mito_snv_vcf}.breakmulti ${mito_snv_vcf}.fix
		bcftools sort ${mito_snv_vcf}.breakmulti | bgzip > ${mito_snv_vcf}.breakmulti.fix
		tabix -p vcf ${mito_snv_vcf}.breakmulti.fix
		bcftools norm -f ${val_rcrs_fasta} -o ${mito_snv_vcf.baseName}.adjusted.vcf ${mito_snv_vcf}.breakmulti.fix
		bcftools view -i 'FMT/AF[*]>0.05' ${mito_snv_vcf.baseName}.adjusted.vcf -o ${group}.mutect2.breakmulti.filtered5p.vcf
		bcftools filter -S 0 --exclude 'FMT/AF[*]<0.05' ${group}.mutect2.breakmulti.filtered5p.vcf -o ${group}.mutect2.breakmulti.filtered5p.0genotyped.vcf
		filter_mutect2_mito.pl ${group}.mutect2.breakmulti.filtered5p.0genotyped.vcf ${meta.id} > ${group}.mutect2.breakmulti.filtered5p.0genotyped.proband.vcf

		${split_normalize_mito_version(task)}
		"""

	stub:
		"""
		echo "${meta.id}" > proband.id
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
// TODO: what is this future merging issue and does it still apply?
process run_hmtnote {
	cpus 2
	memory '5GB'
	time '1h'
	// TODO: Add an explicit container for HmtNote and remove the conda activate dependency.

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.fixinfo.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		source activate tools
		hmtnote annotate ${vcf} ${group}.hmtnote --offline
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
	// TODO: Add an explicit container for bcftools, Haplogrep, graphviz, ghostscript, and montage.
	publishDir "${val_results_output_dir}/plots/mito", mode: 'copy', overwrite: true, pattern: '*.png'

	input:
		tuple val(group), val(id), path(mito_snv_vcf)
		val val_results_output_dir
		val val_accessdir

	output:
		path "${group}.haplogrep.png"
		tuple val(group), path("${group}_haplo.INFO"), emit: haplogrep_INFO
		path "*versions.yml", emit: versions

	script:
		version_str = run_haplogrep_version(task)
		"""
		for sample in \$(bcftools query -l "${mito_snv_vcf}"); do

			bcftools view -c1 -Oz -s "\$sample" -o "\${sample}.vcf.gz" "${mito_snv_vcf}"
			java  -Xmx16G -Xms16G -jar /opt/bin/haplogrep.jar classify \
			--in "\${sample}.vcf.gz" \\
			--out "\${sample}.hg2.vcf" \\
			--format vcf \\
			--lineage 1

			dot "\${sample}.hg2.vcf.dot" -Tps2 > "\${sample}.hg2.vcf.ps"

			gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dGraphicsAlphaBits=4 -r1200 -dDownScaleFactor=3 -sOutputFile=\${sample}.hg2.vcf.png \${sample}.hg2.vcf.ps

		done
		montage -mode concatenate -tile 3x1 *.png ${group}.haplogrep.png
		echo "IMG haplogrep ${val_accessdir}/plots/mito/${group}.haplogrep.png" > "${group}_haplo.INFO"

		echo "${version_str}" > "${task.process}_versions.yml"
		"""

	stub:
		version_str = run_haplogrep_version(task)
		"""
		touch "${group}.haplogrep.png"
		touch "${group}_haplo.INFO"

		echo "${version_str}" > "${task.process}_versions.yml"
		"""
}
def run_haplogrep_version(task) {
	// TODO: Reconcile this version stub with others.
	"""${task.process}:
	    haplogrep: \$(echo \$(java -jar /opt/bin/haplogrep.jar classify 2>&1) | sed "s/htt.*Classify v// ; s/ .*//")
	    montage: \$(echo \$(gm -version 2>&1) | head -1 | sed -e "s/GraphicsMagick //" | cut -d" " -f1 )"""
}

// use eKLIPse for detecting mitochondrial deletions
process run_eklipse {

	tag "$id"
	cpus 2
	// in rare cases with samples above 50 000x this can peak at 500+ GB of VMEM. Add downsampling!
	// TODO: Add an explicit eKLIPse container and remove the conda activate dependency.
	// TODO: Split deletion calling, plotting, and hetplasmid frequency post-processing into separate processes.
	memory '100GB'
	time '60m'
	publishDir "${val_results_output_dir}/plots/mito", mode: 'copy', overwrite: true, pattern: '*.txt'
	publishDir "${val_results_output_dir}/plots/mito", mode: 'copy', overwrite: true, pattern: '*.png'

	input:
		tuple val(group), val(id), val(meta), path(bam), path(bai)
		val val_results_output_dir
		val val_accessdir

	output:
		tuple path("*.png"), path("${id}.hetplasmid_frequency.txt")
		tuple val(group), path("${id}_eklipse.INFO"), emit: eklipse_INFO, optional: true
		path "*versions.yml", emit: versions

	script:
		yml_info_command = ""
		if (meta.type == "proband") {
			yml_info_command = "echo 'IMG eklipse ${val_accessdir}/plots/mito/${id}_eklipse.png' > ${id}_eklipse.INFO"
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
		if (meta.type == "proband") {
			yml_info_command = "echo 'IMG eklipse ${val_accessdir}/plots/mito/${id}_eklipse.png' > ${id}_eklipse.INFO"
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
