process gnomad_sv_bedtools {
    tag "$group"
    cpus 1
    memory '4 GB'
    time '1h'
    // The existing base image provides Python 3 and bedtools.
    container "${params.container}"

    input:
    tuple val(group), val(id), path(vcf)
    path database
    path synonyms

    output:
    tuple val(group), val(id), path("${group}.gnomad.vcf"), emit: vcf
    path "*versions.yml", emit: versions

    script:
    """
    annotate_gnomad_sv_bedtools.py --vcf "${vcf}" --database "${database}" --synonyms "${synonyms}" --output "${group}.gnomad.vcf"
    ${gnomad_sv_bedtools_version(task)}
    """

    stub:
    """
    touch ${group}.gnomad.vcf
    ${gnomad_sv_bedtools_version(task)}
    """
}

def gnomad_sv_bedtools_version(task) {
    """
    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
