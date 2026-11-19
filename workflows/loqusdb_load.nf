include { vcfHasVariants } from './util.nf'

workflow LOQUSDB_LOAD {
    take:
    ch_ped                  // [group, type, ped]
    ch_snv                  // [group, type, vcf, index]
    ch_sv                   // [group, vcf or []]
    ch_loqusdb_database     // [group, database]
    val_accessdir

    main:
    ch_sv_for_load = ch_sv.map { group, vcf ->
        def has_variants = vcf && vcf.exists() && vcf.size() > 0 && vcfHasVariants(vcf)
        tuple(group, has_variants ? vcf : [])
    }

    ch_load_input = ch_ped
        .join(ch_snv, by: [0, 1])
        .join(ch_sv_for_load, by: 0)
        .join(ch_loqusdb_database, by: 0)
        .map { group, _type, ped, snv, _index, sv, database ->
            tuple(group, ped, snv, sv, database)
        }

    add_to_loqusdb(ch_load_input, val_accessdir)

    emit:
    loqusdb_done = add_to_loqusdb.out.loqusdb_done
}

process add_to_loqusdb {
    cpus 1
    publishDir "${params.crondir}/loqus", mode: 'copy', overwrite: true
    tag "$group"
    memory '100 MB'
    time '25m'

    input:
    tuple val(group), path(ped), path(vcf_snvs), path(vcf_svs), val(database)
    val accessdir

    output:
    path("${group}*.loqus"), emit: loqusdb_done

    script:
    def sv_variants_arg = vcf_svs ? "--sv-variants ${accessdir}/sv_vcf/merged/${vcf_svs}" : ""

    """
    echo "-db ${database} -g GRCh38 load -f ${accessdir}/ped/${ped} --variant-file ${accessdir}/vcf/${vcf_snvs} ${sv_variants_arg}" > "${group}.loqus"
    """

    stub:
    """
    touch "${group}.loqus"
    """
}
