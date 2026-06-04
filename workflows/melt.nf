
workflow MELT {

    take:
    ch_bam_bai       // ch:    [ group, id, bam, bai ]
    ch_qc_vals       // ch:    [ group, id, qc       ]
    vcf_header       // value: path(vcf_header)
    bed_intersect    // value: path(bed_intersect)

    main:
    ch_versions = channel.empty()

    ch_bam_bai
        .join(ch_qc_vals, by: [0, 1])
        .set{ ch_melt_in }

    melt(ch_melt_in)
    intersect_melt(melt.out.melt_vcf_nonfiltered)

    ch_versions = ch_versions.mix(melt.out.versions.first())
    ch_versions = ch_versions.mix(intersect_melt.out.versions.first())

    emit:
    melt_intersect_vcf = intersect_melt.out.merged_intersected_vcf
    versions           = ch_versions
    
}

// MELT always give VCFs for each type of element defined in mei_list
// If none found -> 0 byte vcf. merge_melt.pl merges the three, if all empty
// it creates a vcf with only header
// merge_melt.pl gives output ${id}.melt.merged.vcf
process melt {
	cpus 3
	errorStrategy 'retry'
	container  "${params.container_melt}"
	tag "$id"
	// memory seems to scale with less number of reads?
	memory '70 GB'
	time '3h'

	input:
		tuple val(group), val(id), path(bam), path(bai), val(qc)

	output:
		tuple val(group), val(id), path("${id}.melt.merged.vcf"), emit: melt_vcf_nonfiltered
		path "*versions.yml", emit: versions

	script:
		"""
		java -jar /opt/MELTv2.2.2/MELT.jar Single \\
			-bamfile $bam \\
			-r 150 \\
			-h ${params.genome_file} \\
			-n /opt/MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed \\
			-z 500000 \\
			-d 50 \\
			-t $params.mei_list \\
			-w . \\
			-c ${qc.mean_depth} \\
			-e ${qc.ins_size} \\
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
	publishDir "${params.results_output_dir}/vcf", mode: 'copy' , overwrite: true, pattern: '*.vcf'

	input:
		tuple val(group), val(id), path(vcf)

	output:
		tuple val(group), val(id), path("${id}.melt.merged.intersected.vcf"), emit: merged_intersected_vcf
		path "*versions.yml", emit: versions

	when:
		params.run_melt

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
