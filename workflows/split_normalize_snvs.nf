workflow SPLIT_NORMALIZE_SNVS {

    take:
    ch_family_snv_vcf_idx // channel: [val(group_id), path(family_vcf), path(family_tbi)]
    bed_intersect         // value:   [path(bed_intersect)]

    main:
    ch_versions = channel.empty()

	params.results_output_dir = params.outdir + '/' + params.subdir
    
    vcflib_vcfbreakmulti(ch_family_snv_vcf_idx)
    bcftools_norm_sort(vcflib_vcfbreakmulti.out.vcf_tbi)
    vcflib_vcfuniq(bcftools_norm_sort.out.vcf_tbi)
    wgs_dpaf_filter(vcflib_vcfuniq.out.vcf_tbi)
    bedtools_intersect(wgs_dpaf_filter.out.vcf_tbi, bed_intersect)
    bgzip_tabix(bedtools_intersect.out.intersected_vcf)
    
    ch_versions = ch_versions.mix(vcflib_vcfbreakmulti.out.versions.first())
    ch_versions = ch_versions.mix(bcftools_norm_sort.out.versions.first())
    ch_versions = ch_versions.mix(vcflib_vcfuniq.out.versions.first())
    ch_versions = ch_versions.mix(wgs_dpaf_filter.out.versions.first())
    ch_versions = ch_versions.mix(bedtools_intersect.out.versions.first())
    ch_versions = ch_versions.mix(bgzip_tabix.out.versions.first())

    emit:
    vcf_tbi_full            = wgs_dpaf_filter.out.vcf_tbi // channel: [ val(group), path(vcf), path(tbi)]
    vcf_tbi_intersected     = bgzip_tabix.out.vcf_tbi     // channel: [ val(group), path(vcf), path(tbi)]
    versions                = ch_versions                 // channel: [ path(versions) ]
}

process vcflib_vcfbreakmulti {
    cpus 2
	tag "$group"
	memory '10 GB'
	time '1h'

    container "${params.container_vcflib}"
    
    input:
		tuple val(group), path(vcf), path(idx)
    output:
	    tuple val(group), path("${group}.multibreak.vcf.gz"), path("${group}.multibreak.vcf.gz.tbi"), emit: vcf_tbi
        path "*versions.yml", emit: versions

    script:
    """
    vcfbreakmulti ${vcf} | bgzip -c > ${group}.multibreak.vcf.gz
    tabix -p vcf ${group}.multibreak.vcf.gz

    ${vcflib_vcfbreakmulti_version(task)}
    """

    stub:
    """
    touch "${group}.multibreak.vcf.gz"
    touch "${group}.multibreak.vcf.gz.tbi"
    ${vcflib_vcfbreakmulti_version(task)}
    """
}
def vcflib_vcfbreakmulti_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcflib: 1.0.9
	END_VERSIONS
	"""
}

process bcftools_norm_sort {
    cpus 2
	tag "$group"
	memory '50 GB'
	time '1h'

    container "${params.container_bcftools}"
    
    input:
		tuple val(group), path(vcf), path(idx)
    output:
	    tuple val(group), path("${group}.norm.vcf.gz"), path("${group}.norm.vcf.gz.tbi"), emit: vcf_tbi
        path "*versions.yml", emit: versions

    // TODO:  Look into replacing downstreams vcflib_vcfuniq wby moving sort step to top
    //        and adding --rm-dup flag to norm
    script:
    """
    bcftools norm -m-both -c w -f ${params.genome_file} ${vcf} |\
        bcftools sort -O z --write-index=tbi -o ${group}.norm.vcf.gz

    ${bcftools_norm_version(task)}
    """

    stub:
    """
    touch "${group}.norm.vcf.gz"
    touch "${group}.norm.vcf.gz.tbi"
    ${bcftools_norm_version(task)}
    """
}
def bcftools_norm_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS
	"""
}

process vcflib_vcfuniq {
    cpus 2
	tag "$group"
	memory '50 GB'
	time '1h'

    container "${params.container_vcflib}"
    
    input:
		tuple val(group), path(vcf), path(idx)
    output:
	tuple val(group), path("${group}.uniq.vcf.gz"), path("${group}.uniq.vcf.gz.tbi"), emit: vcf_tbi
    path "*versions.yml", emit: versions

    script:
    """
    vcfuniq ${vcf} | bgzip -c > "${group}.uniq.vcf.gz"
    tabix -p vcf "${group}.uniq.vcf.gz" 

    ${vcflib_vcfuniq_version(task)}
    """

    stub:
    """
    touch "${group}.uniq.vcf.gz"
    touch "${group}.uniq.vcf.gz.tbi"
    ${vcflib_vcfuniq_version(task)}
    """
}
def vcflib_vcfuniq_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcflib: 1.0.9
	END_VERSIONS
	"""
}

process wgs_dpaf_filter {
    cpus 2
	publishDir "${params.results_output_dir}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	tag "$group"
	memory '50 GB'
	time '1h'

    container "${params.container_perl}"
    
    input:
		tuple val(group), path(vcf), path(idx)
    output:
    tuple val(group), path("${group}.DPAF.vcf.gz"), path("${group}.DPAF.vcf.gz.tbi"), emit: vcf_tbi
    path "*versions.yml", emit: versions
    
    script:
    """
    wgs_DPAF_filter.pl ${vcf} | bgzip -c > ${group}.DPAF.vcf.gz
    tabix -p vcf ${group}.DPAF.vcf.gz

    ${wgs_dpaf_filter_version(task)}
    """

    stub:
    """
    touch "${group}.DPAF.vcf.gz"
    touch "${group}.DPAF.vcf.gz.tbi"
    ${wgs_dpaf_filter_version(task)}
    """
}
def wgs_dpaf_filter_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bgzip: \$(echo \$(bgzip --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	END_VERSIONS
	"""
}
process bedtools_intersect {
	cpus 2
	tag "$group"
	memory '50 GB'
	time '1h'

    container "${params.container_bedtools}"
    
	input:
		tuple val(group), path(vcf), path(idx) // is ids supposed to be tuple?
	    path(intersect_bed)

	output:
		tuple val(group), path("${group}.intersected.vcf"), emit: intersected_vcf
		path "*versions.yml", emit: versions

    script:
    """
    bedtools intersect \\
	        -a ${vcf} \\
			-b ${intersect_bed} \\
			-u -header > ${group}.intersected.vcf

    ${bedtool_intersect_version(task)}
    """

    stub:
    """
    touch "${group}.intersected.vcf"
    ${bedtool_intersect_version(task)}
    """
}
def bedtool_intersect_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bedtools: \$(echo \$(bedtools --version 2>&1) | sed -e "s/^.*bedtools v//" )
	END_VERSIONS
	"""
}
process bgzip_tabix {
	cpus 2
	tag "$group"
	memory '50 GB'
	time '1h'

    container "${params.container_perl}"
    
	input:
		tuple val(group), path(vcf) // is ids supposed to be tuple?

	output:
	tuple val(group), path("${group}.bgzip.vcf.gz"), path("${group}.bgzip.vcf.gz.tbi"), emit: vcf_tbi
	path "*versions.yml", emit: versions

    script:
    """
    bgzip -c ${vcf} > "${group}.bgzip.vcf.gz"
    tabix -p vcf ${group}.bgzip.vcf.gz

    ${bgzip_tabix_version(task)}
    """

    stub:
    """
    touch "${group}.bgzip.vcf.gz"
    touch "${group}.bgzip.vcf.gz.tbi"
    ${bgzip_tabix_version(task)}
    """
}
def bgzip_tabix_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bgzip: \$(echo \$(bgzip --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	END_VERSIONS
	"""
}
