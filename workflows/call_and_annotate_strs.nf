workflow CALL_AND_ANNOTATE_STRS {

	take:
	ch_bam_bai                   // channel: [ val(group), val(id), path(bam), path(bai) ]
	ch_proband_meta              // channel: [ val(group), val(id), val(meta) ]
	val_analysis_mode            // string:  Analysis mode derived from sample count, either "single" or "family"
	val_genome_file              // path:    Reference FASTA.
	val_genome_fai               // path:    Reference FASTA index.
	val_expansionhunter_catalog  // path:    ExpansionHunter variant catalog JSON.
	val_accessdir                // string:  Base access path used in output metadata/INFO paths

	main:
	ch_versions    = channel.empty()
	ch_output_info = channel.empty()

	expansionhunter(
		ch_bam_bai.join(ch_proband_meta, by: [0,1]),
		val_genome_file,
		val_genome_fai,
		val_expansionhunter_catalog
	)
	stranger(
		expansionhunter.out.expansionhunter_vcf,
		val_expansionhunter_catalog
	)
	rename_sample_expansionhunter_vcf(stranger.out.vcf_annotated.join(ch_proband_meta, by:[0,1]))
	vcfbreakmulti_expansionhunter(rename_sample_expansionhunter_vcf.out.vcf)

	ch_split_expansionhunter_vcf = vcfbreakmulti_expansionhunter.out.vcf

    // Scout requires all individuals in a trio to be present in the STR VCF.
    // Here we add dummy samples to the proband STR VCF so it can be loaded
    // into scout without the actual parental calls.
	if (val_analysis_mode == "family") {
		familyfy_expansionhunter_vcf(ch_split_expansionhunter_vcf)
		ch_expansionhunter_vcf = familyfy_expansionhunter_vcf.out.vcf
		ch_versions = ch_versions.mix(familyfy_expansionhunter_vcf.out.versions.first())
	} else {
		ch_expansionhunter_vcf = ch_split_expansionhunter_vcf.map { group, id, vcf, _meta ->
			tuple(group, id, vcf)
		}
	}

	bgzip_index_expansionhunter_vcf(ch_expansionhunter_vcf, val_accessdir)
	ch_output_info = ch_output_info.mix(bgzip_index_expansionhunter_vcf.out.str_INFO)

	reviewer_loci = new groovy.json.JsonSlurper()
		.parseText(val_expansionhunter_catalog.text)
		.collect { locus_definition -> locus_definition.LocusId }
		.findAll { locus -> locus }

	ch_reviewer_input = expansionhunter.out.bam_vcf
		.flatMap { group, id, bam, bai, vcf ->
			reviewer_loci.collect { locus ->
				tuple(group, id, bam, bai, vcf, locus)
			}
		}

	reviewer(
		ch_reviewer_input,
		val_genome_file,
		val_genome_fai,
		val_expansionhunter_catalog,
		val_accessdir
	)
	ch_output_info = ch_output_info.mix(reviewer.out.reviewer_INFO)

	ch_versions = ch_versions.mix(reviewer.out.versions.first())
	ch_versions = ch_versions.mix(expansionhunter.out.versions.first())
	ch_versions = ch_versions.mix(vcfbreakmulti_expansionhunter.out.versions.first())
	ch_versions = ch_versions.mix(rename_sample_expansionhunter_vcf.out.versions.first())
	ch_versions = ch_versions.mix(bgzip_index_expansionhunter_vcf.out.versions.first())
	ch_versions = ch_versions.mix(stranger.out.versions.first())

	emit:
	vcf               = bgzip_index_expansionhunter_vcf.out.vcf // channel: [ val(group), path(vcf.gz) ]
	tbi               = bgzip_index_expansionhunter_vcf.out.tbi // channel: [ val(group), path(vcf.gz.tbi) ]
	reviewer_plot_svg = reviewer.out.svg                        // channel: path(svg)
	output_info       = ch_output_info                          // channel: [ val(group), path(INFO) ]
	versions          = ch_versions                             // channel: path(versions.yml)
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

	input:
		tuple val(group), val(id), path(bam), path(bai), val(meta)
		path(fasta)
		path(fai)
		path(expansionhunter_catalog)

	output:
		tuple val(group), val(id), path("${group}.eh.vcf"), emit: expansionhunter_vcf
		tuple val(group), val(id), path("${group}.eh_realigned.sort.bam"), path("${group}.eh_realigned.sort.bam.bai"), path("${group}.eh.vcf"), emit: bam_vcf
		path "*versions.yml", emit: versions

	script:
		"""
		source activate htslib10
		ExpansionHunter \
			--reads $bam \
			--reference ${fasta} \
			--variant-catalog ${expansionhunter_catalog} \
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
	container  "${params.container_stranger}"

	input:
		tuple val(group), val(id), path(eh_vcf)
		path(expansionhunter_catalog)

	output:
		tuple val(group), val(id), path("${group}.fixinfo.eh.stranger.vcf"), emit: vcf_annotated
		path "*versions.yml", emit: versions

	script:
		"""
		stranger ${eh_vcf} -f ${expansionhunter_catalog} > ${group}.eh.stranger.vcf
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


process reviewer {
	tag "$group:$locus"
	cpus 2
	time '1h'
	memory '1 GB'
	errorStrategy 'ignore'
	container  "${params.container_reviewer}"
	publishDir "${params.outdir}/${params.subdir}/plots/reviewer/${group}", mode: 'copy' , overwrite: true, pattern: '*.svg'

	input:
		tuple val(group), val(id), path(bam), path(bai), path(vcf), val(locus)
		path(fasta)
		path(fai)
		path(expansionhunter_catalog)
		val(accessdir)

	output:
		path("*svg"), emit: svg
		tuple val(group), path("${group}_reviewer.INFO"), emit: reviewer_INFO
		path "*versions.yml", emit: versions

	script:
		version_str = reviewer_version(task)
		"""
		REViewer --reads ${bam} \\
			--vcf ${vcf} \\
			--reference ${fasta} \\
			--catalog ${expansionhunter_catalog} \\
			--locus "${locus}" \\
			--output-prefix ${id}
		echo "STR_VARIANTS_IMG	${locus}	${accessdir}/plots/reviewer/${group}/${id}.${locus}.svg" > ${group}_reviewer.INFO

		echo "${version_str}" > "${task.process}_versions.yml"
		"""

	stub:
		version_str = reviewer_version(task)
		"""
		touch "${id}.${locus}.svg"
		echo "STR_VARIANTS_IMG	${locus}	${accessdir}/plots/reviewer/${group}/${id}.${locus}.svg" > "${group}_reviewer.INFO"
		echo "${version_str}" > "${task.process}_versions.yml"
		"""
}
def reviewer_version(task) {
	// TODO: Reconcile this version stub with others.
	"""${task.process}:
	    reviewer: \$(echo \$(REViewer --version 2>&1) | sed 's/^.*REViewer v//')"""
}

process rename_sample_expansionhunter_vcf {
	container "${params.container_picard}"
	tag "$group"
	time '1h'
	memory '5 GB'

	input:
		tuple val(group), val(id), path(eh_vcf_anno), val(meta)

	output:
		tuple val(group), val(id), path("${group}.eh.rename.vcf"), val(meta), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		picard RenameSampleInVcf \\
			INPUT=${eh_vcf_anno} \\
			OUTPUT=${group}.eh.rename.vcf \\
			NEW_SAMPLE_NAME=${id}

		${rename_sample_expansionhunter_vcf_version(task)}
		"""

	stub:
		"""
		touch "${group}.eh.rename.vcf"
		${rename_sample_expansionhunter_vcf_version(task)}
		"""
}
def rename_sample_expansionhunter_vcf_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    picard: \$( echo \$(picard RenameSampleInVcf --version 2>&1 |grep -i version | cut -f 2 -d ':') | sed 's/-SNAPSHOT//' )
	END_VERSIONS
	"""
}

// split multiallelic sites in expansionhunter vcf
process vcfbreakmulti_expansionhunter {
	container "${params.container_vcflib}"
	tag "$group"
	time '1h'
	memory '5 GB'

	input:
		tuple val(group), val(id), path(eh_vcf), val(meta)

	output:
		tuple val(group), val(id), path("${group}.expansionhunter.vcf"), val(meta), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		vcfbreakmulti ${eh_vcf} > ${group}.expansionhunter.vcf
		${vcfbreakmulti_expansionhunter_version(task)}
		"""

	stub:
		"""
		touch "${group}.expansionhunter.vcf"
		${vcfbreakmulti_expansionhunter_version(task)}
		"""
}
def vcfbreakmulti_expansionhunter_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcflib: 1.0.9
	END_VERSIONS
	"""
}

process familyfy_expansionhunter_vcf {
	container "${params.container_perl}"
	tag "$group"
	time '20m'
	memory '1 GB'

	input:
		tuple val(group), val(id), path(eh_vcf), val(meta)

	output:
		tuple val(group), val(id), path("${group}.expansionhunter.family.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		def father_arg = meta.father ? "--father '${meta.father}'" : ""
		def mother_arg = meta.mother ? "--mother '${meta.mother}'" : ""
		"""
		familyfy_str.pl --vcf ${eh_vcf} ${mother_arg} ${father_arg} --out ${group}.expansionhunter.family.vcf
		${familyfy_expansionhunter_vcf_version(task)}
		"""

	stub:
		"""
		touch "${group}.expansionhunter.family.vcf"
		${familyfy_expansionhunter_vcf_version(task)}
		"""
}
def familyfy_expansionhunter_vcf_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    perl: \$(perl -e 'printf "%vd\\n", \$^V')
	END_VERSIONS
	"""
}

process bgzip_index_expansionhunter_vcf {
	container "${params.container_bcftools}"
	publishDir "${params.outdir}/${params.subdir}/vcf", mode: 'copy' , overwrite: true, pattern: '*.vcf.gz'
	publishDir "${params.outdir}/${params.subdir}/vcf", mode: 'copy' , overwrite: true, pattern: '*.vcf.gz.tbi'
	tag "$group"
	time '1h'
	memory '5 GB'

	input:
		tuple val(group), val(id), path(eh_vcf)
		val(accessdir)

	output:
		tuple val(group), path("${group}.expansionhunter.vcf.gz"), emit: vcf
		tuple val(group), path("${group}.expansionhunter.vcf.gz.tbi"), emit: tbi
		tuple val(group), path("${group}_str.INFO"), emit: str_INFO
		path "*versions.yml", emit: versions

	script:
		"""
		cp ${eh_vcf} ${group}.expansionhunter.vcf
		bgzip ${group}.expansionhunter.vcf
		tabix ${group}.expansionhunter.vcf.gz
		echo "STR	${accessdir}/vcf/${group}.expansionhunter.vcf.gz" > ${group}_str.INFO

		${bgzip_index_expansionhunter_vcf_version(task)}
		"""

	stub:
		"""
		touch "${group}.expansionhunter.vcf.gz"
		touch "${group}.expansionhunter.vcf.gz.tbi"
		touch "${group}_str.INFO"
		${bgzip_index_expansionhunter_vcf_version(task)}
		"""
}
def bgzip_index_expansionhunter_vcf_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bgzip: \$(echo \$(bgzip --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	END_VERSIONS
	"""
}
