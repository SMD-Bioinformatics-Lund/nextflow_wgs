
workflow MELT {

    take:
    ch_bam_bai       // ch:    [ group, id, bam, bai   ]
    ch_mean_depth    // ch:    [ group, id, mean_depth ]
    ch_ins_size      // ch:    [ group, id, ins_size   ]
    fasta            // value: path(reference.fa)
    fai              // value: path(reference.fa.fai)
    mei_list         // value: path(mei_list.txt)
    vcf_header       // value: path(vcf_header)
    bed_intersect    // value: path(bed_intersect)

    main:
    ch_versions = channel.empty()

    ch_melt_in = ch_bam_bai
        .join(ch_mean_depth, by: [0, 1])
        .join(ch_ins_size, by: [0, 1])

    melt(ch_melt_in, fasta, fai, mei_list)
    merge_melt(melt.out.melt_vcfs, vcf_header)
    intersect_melt(merge_melt.out.melt_vcf_nonfiltered, bed_intersect)

    ch_versions = ch_versions.mix(melt.out.versions.first())
    ch_versions = ch_versions.mix(merge_melt.out.versions.first())
    ch_versions = ch_versions.mix(intersect_melt.out.versions.first())

    emit:
    melt_intersect_vcf = intersect_melt.out.merged_intersected_vcf
    versions           = ch_versions
    
}

process melt {
	cpus 3
	errorStrategy 'retry'
	container  "${params.container_melt}"
	tag "$id"
	// memory seems to scale with less number of reads?
	memory '70 GB'
	time '3h'

	input:
		tuple val(group), val(id), path(bam), path(bai), val(mean_depth), val(ins_size)
		path(fasta)
		path(fasta_fai)
		path(mei_list)
    
	output:
		tuple val(group), val(id), path("ALU.final_comp.vcf"), path("LINE1.final_comp.vcf"), path("SVA.final_comp.vcf"), emit: melt_vcfs
		path "*versions.yml", emit: versions

	script:
		"""
		java -jar /opt/MELTv2.2.2/MELT.jar Single \\
			-bamfile $bam \\
			-r 150 \\
			-h $fasta \\
			-n /opt/MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed \\
			-z 500000 \\
			-d 50 \\
			-t $mei_list \\
			-w . \\
			-c $mean_depth \\
			-e $ins_size \\
			-exome

		${melt_version(task)}
		"""

	stub:
		"""
		touch ALU.final_comp.vcf
		touch LINE1.final_comp.vcf
		touch SVA.final_comp.vcf
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

// MELT always gives VCFs for each type of element defined in mei_list.
// If none are found, merge_melt.pl creates a VCF with only header.
process merge_melt {
	cpus 1
	container  "${params.container_melt}"
	tag "$id"
	memory '2 GB'
	time '1h'
	publishDir "${params.results_output_dir}/vcf", mode: 'copy' , overwrite: true, pattern: '*.vcf'

	input:
		tuple val(group), val(id), path(alu_vcf), path(line1_vcf), path(sva_vcf)
		path vcf_header

	output:
		tuple val(group), val(id), path("${id}.melt.merged.vcf"), emit: melt_vcf_nonfiltered
		path "*versions.yml", emit: versions

	script:
		"""
		merge_melt.pl \\
			--vcf-header $vcf_header \\
			--id $id \\
			--alu $alu_vcf \\
			--line1 $line1_vcf \\
			--sva $sva_vcf \\
			--out ${id}.melt.merged.vcf

		${merge_melt_version(task)}
		"""

	stub:
		"""
		touch "${id}.melt.merged.vcf"
		${merge_melt_version(task)}
		"""
}
def merge_melt_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcftools: \$(vcf-concat --version 2>&1 | head -n 1 | sed 's/^.*vcf-concat //')
	END_VERSIONS
	"""
}

process intersect_melt {
	cpus 2
	tag "$id"
	memory '2 GB'
	time '1h'
	publishDir "${params.results_output_dir}/vcf", mode: 'copy' , overwrite: true, pattern: '*.vcf'
    container "${params.container_bedtools"}

	input:
		tuple val(group), val(id), path(vcf)
		path bed_intersect

	output:
		tuple val(group), val(id), path("${id}.melt.merged.intersected.vcf"), emit: merged_intersected_vcf
		path "*versions.yml", emit: versions

	when:
		params.run_melt

	script:
		"""
		bedtools intersect -a $vcf -b $bed_intersect -header > ${id}.melt.merged.intersected.vcf
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
