// Population frequencies use SVDB before VEP, as in Nallo's ANNOTATE_SVS.
process GNOMAD_SV {
    tag "$group"
    cpus 2
    memory '10 GB'
    time '10h'
    container "${params.container_svdb_2_12_2}"

    input:
    tuple val(group), val(id), path(vcf)
    path(database)

    output:
    tuple val(group), val(id), path("${group}.gnomad.vcf"), emit: vcf
    path "*versions.yml", emit: versions

    script:
    """
    # N_HET is already retained in the slim reference. It is copied as a count;
    # matching uses geometry/type, and AF is copied directly from the best match.
    # VEP used 70% reciprocal overlap with no breakpoint-distance cutoff.
    # SVDB uses intersection/span instead. A 50 kb endpoint-distance limit
    # allows larger shifts while keeping candidate lookups local.
    svdb \\
        --query \\
        --db ${database} \\
        --in_occ N_HET --in_frq AF \\
        --out_occ gnomad_svN_HET --out_frq gnomad_svAF \\
        --overlap 0.6 --bnd_distance 50000 --data_profile cohort \\
        --query_vcf ${vcf} > ${group}.gnomad.vcf

    # Ensure annotation does not lose query records.
    if ! awk -F '\\t' 'NR==FNR {if (!/^#/) expected++; next} !/^#/ {if (NF<8) invalid=1; observed++} END {exit invalid || expected!=observed}' ${vcf} ${group}.gnomad.vcf; then
        echo 'SVDB did not retain every query record. Check the query VCF and SVDB log.' >&2
        exit 1
    fi

    ${gnomad_sv_version(task)}
    """

    stub:
    """
    touch ${group}.gnomad.vcf
    ${gnomad_sv_version(task)}
    """
}

def gnomad_sv_version(task) {
    """
    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
        svdb: \$(python3 -c 'from importlib.metadata import version; print(version("svdb"))')
    END_VERSIONS
    """
}
