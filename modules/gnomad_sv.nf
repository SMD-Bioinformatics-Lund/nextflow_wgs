// Population frequencies use SVDB before VEP, as in Nallo's ANNOTATE_SVS.
process gnomad_sv {
    tag "$group"
    cpus 2
    memory '10 GB'
    time '10h'
    container "${params.container_svdb}"

    input:
    tuple val(group), val(id), path(vcf)
    path(database)

    output:
    tuple val(group), val(id), path("${group}.gnomad.vcf"), emit: vcf
    path "*versions.yml", emit: versions

    script:
    """
    svdb \\
        --query \\
        --db ${database} \\
        --in_occ AC --in_frq AF \\
        --out_occ gnomad_svAC --out_frq gnomad_svAF \\
        --overlap 0.6 --bnd_distance 10000 --ins_distance 50 \\
        --query_vcf ${vcf} > ${group}.gnomad.vcf

    # SVDB 2.8 can exit successfully with only a header if source AC/AF is missing.
    if ! awk -F '\\t' 'NR==FNR {if (!/^#/) expected++; next} !/^#/ {if (NF<8) invalid=1; observed++} END {exit invalid || expected!=observed}' ${vcf} ${group}.gnomad.vcf; then
        echo 'SVDB did not retain every query record. Check source AC/AF; regenerate older gnomAD slim resources.' >&2
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
        svdb: \$(svdb 2>&1 | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/')
    END_VERSIONS
    """
}
