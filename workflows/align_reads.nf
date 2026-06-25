nextflow.enable.dsl=2

workflow ALIGN_READS {

	take:
	ch_fastq
	genome_file

	main:
	ch_genome = channel.value(file(genome_file))

	minibwa_index(ch_genome)

	ch_genome_index = ch_genome.combine(minibwa_index.out.index)
	ch_fastq_with_index = ch_fastq.combine(ch_genome_index)

	minibwa_align(ch_fastq_with_index)
	ch_sam_with_genome = minibwa_align.out.sam.combine(ch_genome)
	sam_to_sorted_bam(ch_sam_with_genome)
	markdup(sam_to_sorted_bam.out.bam_bai)

	ch_versions = minibwa_index.out.versions.first()
		.mix(minibwa_align.out.versions.first())
		.mix(sam_to_sorted_bam.out.versions.first())
		.mix(markdup.out.versions.first())

	emit:
	dedup_bam_bai = markdup.out.dedup_bam_bai
	dedup_metrics = markdup.out.dedup_metrics
	dedup_bam_INFO = markdup.out.dedup_bam_INFO
	versions = ch_versions
}

process minibwa_index {
	cpus 8
	memory '70 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	tag "${genome_file.simpleName}"
	container "${params.container_minibwa}"

	input:
		path genome_file

	output:
		tuple path("*.l2b"), path("*.mbw"), emit: index
		path "*versions.yml", emit: versions

	script:
		"""
		minibwa index -t${task.cpus} ${genome_file}

		${minibwa_index_versions(task)}
		"""

	stub:
		"""
		touch "${genome_file}.l2b"
		touch "${genome_file}.mbw"

		${minibwa_index_versions(task)}
		"""
}

def minibwa_index_versions(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    minibwa: \$(echo \$(minibwa --version 2>&1))
	END_VERSIONS
	"""
}

process minibwa_align {
	cpus 50
	memory '100 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	tag "$id"
	container "${params.container_minibwa}"

	input:
		tuple val(group), val(id), path(fastq_r1), path(fastq_r2), path(genome_file), path(index_l2b), path(index_mbw)

	output:
		tuple val(group), val(id), path("${id}.sam"), emit: sam
		path "*versions.yml", emit: versions

	script:
		"""
		minibwa map \\
			-t${task.cpus} \\
			-R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' \\
			-o ${id}.sam \\
			${genome_file} ${fastq_r1} ${fastq_r2}

		${minibwa_align_versions(task)}
		"""

	stub:
		"""
		touch "${id}.sam"

		${minibwa_align_versions(task)}
		"""
}

def minibwa_align_versions(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    minibwa: \$(echo \$(minibwa --version 2>&1))
	END_VERSIONS
	"""
}

process sam_to_sorted_bam {
	cpus 50
	memory '100 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	tag "$id"
	container "${params.container_sentieon}"

	input:
		tuple val(group), val(id), path(sam), path(genome_file)

	output:
		tuple val(group), val(id), path("${id}_merged.bam"), path("${id}_merged.bam.bai"), emit: bam_bai
		path "*versions.yml", emit: versions

	script:
		"""
		sentieon util sort \\
			-r ${genome_file} \\
			-o ${id}_merged.bam \\
			-t ${task.cpus} --sam2bam -i ${sam}

		${sam_to_sorted_bam_versions(task)}
		"""

	stub:
		"""
		touch "${id}_merged.bam"
		touch "${id}_merged.bam.bai"

		${sam_to_sorted_bam_versions(task)}
		"""
}

def sam_to_sorted_bam_versions(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon util --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}

process markdup {
	cpus 40
	errorStrategy 'retry'
	maxErrors 5
	tag "$id"
	memory '50 GB'
	time '3h'
	container "${params.container_sentieon}"
	publishDir "${params.outdir}/${params.subdir}/bam", mode: 'copy' , overwrite: true, pattern: '*_dedup.bam*'

	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${id}_dedup.bam"), path("${id}_dedup.bam.bai"), emit: dedup_bam_bai
		tuple val(group), val(id), path("dedup_metrics.txt"), emit: dedup_metrics
		tuple val(group), path("${group}_bam.INFO"), emit: dedup_bam_INFO
		path "*versions.yml", emit: versions

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
