#!/usr/bin/env nextflow

include { annotate_vep; vcfanno; modify_vcf; mark_splice; inher_models; genmodscore } from './workflows/annotate_snvs.nf'

params.variant = null
params.query_genotype = '0/1'
params.query_sex = '0'
params.score_outdir = 'variant-score'

process PREPARE_VARIANT_QUERY {
    container params.container_python
    cpus 1
    memory '1 GB'
    time '5m'

    input:
    path 'variant.input'
    path 'reference.fa'
    path 'reference.fa.fai'
    path 'prepare_variant_query.py'
    val genotype
    val sex

    output:
    path 'snv.vcf', emit: snv
    path 'query.ped', emit: ped

    script:
    """
    python3 prepare_variant_query.py \
        --variant variant.input --reference reference.fa \
        --genotype '${genotype}' --sex '${sex}' --output-dir .
    """
}

process REPORT_VARIANT_SCORE {
    container params.container_python
    cpus 1
    memory '1 GB'
    time '5m'
    publishDir params.score_outdir, mode: 'copy', overwrite: true

    input:
    path 'ranked.vcf'
    path 'prepare_variant_query.py'
    path 'report_variant_score.py'

    output:
    path 'score.tsv'
    path 'scored.vcf'

    script:
    """
    python3 report_variant_score.py ranked.vcf > score.tsv
    cp ranked.vcf scored.vcf
    """
}

workflow {
    if (!params.variant || !params.genome_file) {
        error 'Provide --variant and configure the GRCh38 genome_file (for example with the Hopper wgs profile).'
    }
    def genotype = params.query_genotype.toString()
    def sex = params.query_sex.toString()
    if (!(genotype in ['0/1', '1/1']) || !(sex in ['0', '1', '2'])) {
        error 'Use --query_genotype 0/1 or 1/1 and --query_sex 0, 1 or 2.'
    }

    PREPARE_VARIANT_QUERY(
        channel.of(file(params.variant, checkIfExists: true)),
        file(params.genome_file, checkIfExists: true),
        file("${params.genome_file}.fai", checkIfExists: true),
        file("${projectDir}/bin/prepare_variant_query.py"),
        genotype,
        sex
    )
    ped = PREPARE_VARIANT_QUERY.out.ped

    snv_input = PREPARE_VARIANT_QUERY.out.snv.map { tuple('query', it) }
    annotate_vep(snv_input)
    vcfanno(annotate_vep.out.vcf)
    modify_vcf(vcfanno.out.vcf)
    mark_splice(modify_vcf.out.vcf)
    snv_models_input = mark_splice.out.splice_marked.combine(ped)
        .map { group, vcf, pedigree -> tuple(group, vcf, 'single', pedigree) }
    inher_models(snv_models_input)
    genmodscore(inher_models.out.vcf, false)

    REPORT_VARIANT_SCORE(
        genmodscore.out.scored_vcf.map { group, type, vcf -> vcf },
        file("${projectDir}/bin/prepare_variant_query.py"),
        file("${projectDir}/bin/report_variant_score.py")
    )
}
