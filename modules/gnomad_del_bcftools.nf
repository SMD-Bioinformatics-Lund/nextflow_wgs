workflow gnomad_del_bcftools {
    take:
    vcf
    database
    synonyms

    main:
    gnomad_del_intervals(vcf, synonyms)
    gnomad_del_annotate(gnomad_del_intervals.out.intervals, database)
    gnomad_del_restore(gnomad_del_annotate.out.annotated)

    emit:
    vcf = gnomad_del_restore.out.vcf
    versions = gnomad_del_annotate.out.versions
}

process gnomad_del_intervals {
    tag "$group"
    cpus 1
    memory '2 GB'
    time '1h'
    container "${params.container}"

    input:
    tuple val(group), val(id), path(vcf)
    path synonyms

    output:
    tuple val(group), val(id), path(vcf), path("${group}.del.intervals.vcf"), emit: intervals

    script:
    """
    gnomad_del.py query --vcf "${vcf}" --synonyms "${synonyms}" --output "${group}.del.intervals.vcf"
    """

    stub:
    """
    touch ${group}.del.intervals.vcf
    """
}

process gnomad_del_annotate {
    tag "$group"
    cpus 1
    memory '4 GB'
    time '1h'
    container "${params.container_bcftools}"

    input:
    tuple val(group), val(id), path(vcf), path(intervals)
    path database

    output:
    tuple val(group), val(id), path(vcf), path("${group}.del.matches.vcf"), emit: annotated
    path "*versions.yml", emit: versions

    script:
    // The helper applies an exact integer 70% check after candidate retrieval.
    // This tolerance prevents single-precision rounding from losing boundary hits.
    """
    bcftools annotate --no-version \\
        --annotations "${database}" \\
        --columns CHROM,FROM,TO,GNOMAD_DEL_MATCHES \\
        --merge-logic GNOMAD_DEL_MATCHES:append \\
        --min-overlap 0.699999:0.699999 \\
        --output-type v --output "${group}.del.matches.vcf" "${intervals}"
    ${gnomad_del_bcftools_version(task)}
    """

    stub:
    """
    touch ${group}.del.matches.vcf
    ${gnomad_del_bcftools_version(task)}
    """
}

process gnomad_del_restore {
    tag "$group"
    cpus 1
    memory '2 GB'
    time '1h'
    container "${params.container}"

    input:
    tuple val(group), val(id), path(vcf), path(matches)

    output:
    tuple val(group), val(id), path("${group}.gnomad.vcf"), emit: vcf

    script:
    """
    gnomad_del.py restore --vcf "${vcf}" --matches "${matches}" --output "${group}.gnomad.vcf"
    """

    stub:
    """
    touch ${group}.gnomad.vcf
    """
}

def gnomad_del_bcftools_version(task) {
    """
    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
