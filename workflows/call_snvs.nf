workflow CALL_SNVS {

    take:
    ch_bam_bai                                // channel: [ val(group), val(id), path(bam), path(bai) ]
    ch_genome_file                            // channel: [ path(reference_fasta)                     ]
    ch_genome_index                           // channel: [ path(reference_fasta_fai)                 ]
    ch_bqsr_known_polymorphic_sites_vcf       // channel: [ path(bqsr_known_polymorphic_sites_vcf)       ]
    ch_bqsr_known_polymorphic_sites_vcf_index // channel: [ path(bqsr_known_polymorphic_sites_vcf_index) ]
    ch_intersect_bed                          // channel: [ path(intersect_bed)                       ]
    ch_vcfanno                                // channel: [ path(vcfanno_config)                      ]
    ch_vcfanno_lua                            // channel: [ path(vcfanno_lua)                         ]
    val_run_freebayes                         // bool:    Whether Freebayes should be run

    main:
    ch_versions = channel.empty()

    bqsr(
        ch_bam_bai,
        ch_genome_file,
        ch_genome_index,
        ch_bqsr_known_polymorphic_sites_vcf,
        ch_bqsr_known_polymorphic_sites_vcf_index
    )

    ch_dnascope_in = ch_bam_bai.join(bqsr.out.bqsr_table, by: [0,1])
    
    dnascope(
        ch_dnascope_in,
        ch_genome_file,
        ch_genome_index
    )
    gvcf_combine(
        dnascope.out.gvcf_tbi.groupTuple(),
        ch_genome_file,
        ch_genome_index
    )

    ch_versions = ch_versions.mix(bqsr.out.versions.first())
    ch_versions = ch_versions.mix(dnascope.out.versions.first())
    ch_versions = ch_versions.mix(gvcf_combine.out.versions.first())

    ch_snv_vcf_tbi = channel.empty()
    
    if(val_run_freebayes) {
        freebayes(
            ch_bam_bai,
            ch_genome_file,
            ch_genome_index,
            ch_intersect_bed,
            ch_vcfanno,
            ch_vcfanno_lua
        )
		ch_versions = ch_versions.mix(freebayes.out.versions.first())
        
        ch_concat_gvcf_freebayes_in = freebayes.out.freebayes_variants
            .join(
                gvcf_combine.out.combined_vcf
            )
            .map {
                group, freebayes_vcf, _id, gvcf, _gvcf_tbi ->
                [ group, gvcf, freebayes_vcf ]
            }

        concat_gvcf_freebayes(ch_concat_gvcf_freebayes_in)
        ch_snv_vcf_tbi = concat_gvcf_freebayes.out.vcf_tbi
        ch_versions = ch_versions.mix(concat_gvcf_freebayes.out.versions.first())
        
    } else {
        ch_snv_vcf_tbi = gvcf_combine.out.combined_vcf
            .map { group, _id, vcf, tbi -> [ group, vcf, tbi ] }
    }

    emit:
    group_vcf_tbi   = ch_snv_vcf_tbi        // channel: [ group, id, vcf, tbi  ]
    sample_gvcf_tbi = dnascope.out.gvcf_tbi // channel: [ group, id, gvcf, tbi ]
    versions        = ch_versions           // channel: [ path(versions.yml)   ]
    
}

process bqsr {
	cpus 40
	errorStrategy 'retry'
	maxErrors 5
	tag "$id"
	memory '30 GB'
	// 12gb peak giab //
	time '5h'
	container  "${params.container_sentieon}"
	publishDir "${params.outdir}/${params.subdir}/bqsr", mode: 'copy' , overwrite: true, pattern: '*.table'

	input:
		tuple val(group), val(id), path(bam), path(bai)
		path(genome_file)
		path(genome_index)
		path(bqsr_known_polymorphic_sites_vcf)
		path(bqsr_known_polymorphic_sites_vcf_index)

	output:
		tuple val(group), val(id), path("${id}.bqsr.table"), emit: bqsr_table
		path "*versions.yml", emit: versions

	script:
		"""
		sentieon driver -t ${task.cpus} \\
			-r $genome_file -i $bam \\
			--algo QualCal ${id}.bqsr.table \\
			-k $bqsr_known_polymorphic_sites_vcf

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

process dnascope {
	cpus 54
	memory '100 GB'
	// 12 GB peak giab //
	time '4h'
	tag "$id"
	container  "${params.container_sentieon}"

	input:
		tuple val(group), val(id), path(bam), path(bai), path(bqsr)
		path(genome_file)
		path(genome_index)

	output:
		tuple val(group), val(id), path("${id}.dnascope.gvcf.gz"), path("${id}.dnascope.gvcf.gz.tbi"), emit: gvcf_tbi
		path "*versions.yml", emit: versions

	// TODO: Move shards out to some config file:
	script:
		"""
		sentieon driver \\
			-t ${task.cpus} \\
			-r ${genome_file} \\
	        -q ${bqsr} \\
	        -i ${bam} \\
			--shard 1:1-248956422  \\
			--shard 2:1-242193529  \\
			--shard 3:1-198295559  \\
			--shard 4:1-190214555  \\
			--shard 5:1-120339935  \\
			--shard 5:120339936-181538259  \\
			--shard 6:1-170805979  \\
			--shard 7:1-159345973  \\
			--shard 8:1-145138636  \\
			--shard 9:1-138394717  \\
			--shard 10:1-133797422  \\
			--shard 11:1-135086622  \\
			--shard 12:1-56232327  \\
			--shard 12:56232328-133275309  \\
			--shard 13:1-114364328  \\
			--shard 14:1-107043718  \\
			--shard 15:1-101991189  \\
			--shard 16:1-90338345  \\
			--shard 17:1-83257441  \\
			--shard 18:1-80373285  \\
			--shard 19:1-58617616  \\
			--shard 20:1-64444167  \\
			--shard 21:1-46709983  \\
			--shard 22:1-50818468  \\
			--shard X:1-124998478  \\
			--shard X:124998479-156040895  \\
			--shard Y:1-57227415  \\
			--shard M:1-16569 \\
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
process gvcf_combine {
	cpus 16
	tag "$group"
	memory '5 GB'
	time '5h'
	container  "${params.container_sentieon}"

	input:
	tuple val(group), val(id), path(gvcfs), path(gvcf_idxs)
    path(genome_file)
    path(genome_index)

	output: // Off to split_normalize, together with other stuff
		tuple val(group), val(id), path("${group}.combined.vcf.gz"), path("${group}.combined.vcf.gz.tbi"), emit: combined_vcf
		path "*versions.yml", emit: versions

	script:
		all_gvcfs = gvcfs.collect { gvcf -> gvcf.toString() }.sort().join(' -v ')
		"""
		sentieon driver \\
			-t ${task.cpus} \\
			-r ${genome_file} \\
			--algo GVCFtyper \\
			-v $all_gvcfs ${group}.combined.vcf.gz

		${gvcf_combine_version(task)}
		"""

	stub:
		all_gvcfs = gvcfs.collect { gvcf -> gvcf.toString() }.sort().join(' -v ')
		"""
		touch "${group}.combined.vcf.gz"
		touch "${group}.combined.vcf.gz.tbi"

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
process freebayes {
	cpus 1
	time '2h'
	memory '10 GB'
	container  "${params.container_twist_myeloid}"

	input:
		tuple val(group), val(id), path(bam), path(bai)
		path(genome_file)
		path(genome_index)
		path(bed_intersect)
		path(vcfanno)
		path(vcfanno_lua)

	output:
		tuple val(group), path("${id}.pathfreebayes.vcf_no_header.tsv.gz"), emit: freebayes_variants
		path "*versions.yml", emit: versions

	script:
		"""
		freebayes -f ${genome_file} --pooled-continuous --pooled-discrete -t $bed_intersect --min-repeat-entropy 1 -F 0.03 $bam > ${id}.freebayes.vcf
		vcfbreakmulti ${id}.freebayes.vcf > ${id}.freebayes.multibreak.vcf
		bcftools norm -m-both -c w -O v -f ${genome_file} -o ${id}.freebayes.multibreak.norm.vcf ${id}.freebayes.multibreak.vcf
		vcfanno_linux64 -lua $vcfanno_lua $vcfanno ${id}.freebayes.multibreak.norm.vcf > ${id}.freebayes.multibreak.norm.anno.vcf
		grep ^# ${id}.freebayes.multibreak.norm.anno.vcf > ${id}.freebayes.multibreak.norm.anno.path.vcf
		grep -v ^# ${id}.freebayes.multibreak.norm.anno.vcf | grep -i pathogenic > ${id}.freebayes.multibreak.norm.anno.path.vcf2
		cat ${id}.freebayes.multibreak.norm.anno.path.vcf ${id}.freebayes.multibreak.norm.anno.path.vcf2 > ${id}.freebayes.multibreak.norm.anno.path.vcf3
		filter_freebayes.pl ${id}.freebayes.multibreak.norm.anno.path.vcf3 | bgzip -c > "${id}.pathfreebayes.vcf_no_header.tsv.gz"

		${freebayes_version(task)}
		"""
	stub:
		"""
		touch "${id}.pathfreebayes.vcf_no_header.tsv.gz"

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

process concat_gvcf_freebayes {
    cpus 2
	tag "$group"
	memory '10 GB'
	time '1h'

    container "${params.container_bcftools}"
    
    input:
	tuple val(group), path(gvcf), path(freebayes_headerless_vcf)
    output:
    tuple val(group), path("${group}.sorted.vcf.gz"), path("${group}.sorted.vcf.gz.tbi"), emit: vcf_tbi
    path "*versions.yml", emit: versions
    
    script:
    """
    zcat ${gvcf} ${freebayes_headerless_vcf} |\
        bcftools sort -o ${group}.sorted.vcf.gz -O z --write-index=tbi
    ${concat_gvcf_freebayes_version(task)}
    """

    stub:
    """
    touch ${group}.sorted.vcf.gz
    touch ${group}.sorted.vcf.gz.tbi
    ${concat_gvcf_freebayes_version(task)}
    """
}
def concat_gvcf_freebayes_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS
	"""
}
