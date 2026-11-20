workflow gnomad_sv_bcftools {
    take:
    vcf
    database
    synonyms

    main:
    gnomad_sv_intervals(vcf, synonyms)
    gnomad_sv_annotate(gnomad_sv_intervals.out.intervals, database)
    gnomad_sv_restore(gnomad_sv_annotate.out.annotated)

    emit:
    vcf = gnomad_sv_restore.out.vcf
    versions = gnomad_sv_annotate.out.versions
}

process gnomad_sv_intervals {
    tag "$group"
    cpus 1
    memory '2 GB'
    time '1h'
    container "${params.container}"

    input:
    tuple val(group), val(id), path(vcf)
    path synonyms

    output:
    tuple val(group), val(id), path(vcf), path("${group}.*.intervals.vcf"), emit: intervals

    script:
    """
    gnomad_sv.py query --vcf "${vcf}" --synonyms "${synonyms}" --output "${group}"
    """

    stub:
    """
    touch ${group}.overlap.intervals.vcf ${group}.ins_site.intervals.vcf ${group}.ins_length.intervals.vcf
    """
}

process gnomad_sv_annotate {
    tag "$group"
    cpus 1
    memory '4 GB'
    time '1h'
    container "${params.container_bcftools}"

    input:
    tuple val(group), val(id), path(vcf), path(intervals)
    path database

    output:
    tuple val(group), val(id), path(vcf), path("${group}.*.matches.vcf"), emit: annotated
    path "*versions.yml", emit: versions

    script:
    // The helper applies an exact integer 70% check after candidate retrieval.
    // This tolerance prevents single-precision rounding from losing boundary hits.
    """
    for mode in overlap ins_site ins_length; do
        overlap=0.699999:0.699999
        if [ "\$mode" = ins_site ]; then overlap=0:0; fi
        bcftools annotate --no-version \\
            --annotations "${database}/\$mode.bed.gz" \\
            --columns CHROM,FROM,TO,GNOMAD_SV_MATCHES \\
            --merge-logic GNOMAD_SV_MATCHES:append \\
            --min-overlap "\$overlap" \\
            --output-type v --output "${group}.\$mode.matches.vcf" "${group}.\$mode.intervals.vcf"
    done
    ${gnomad_sv_bcftools_version(task)}
    """

    stub:
    """
    touch ${group}.overlap.matches.vcf ${group}.ins_site.matches.vcf ${group}.ins_length.matches.vcf
    ${gnomad_sv_bcftools_version(task)}
    """
}

process gnomad_sv_restore {
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
    gnomad_sv.py restore --vcf "${vcf}" --output "${group}.gnomad.vcf" \\
        --matches "${group}.overlap.matches.vcf" "${group}.ins_site.matches.vcf" "${group}.ins_length.matches.vcf"
    """

    stub:
    """
    touch ${group}.gnomad.vcf
    """
}

def gnomad_sv_bcftools_version(task) {
    """
    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
