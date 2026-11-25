#!/usr/bin/env nextflow

include { annotate_vep; vcfanno; modify_vcf; mark_splice; inher_models; genmodscore } from './workflows/annotate_snvs.nf'

params.variant = null
params.hgvs = null
params.coord = null
params.query_genotype = '0/1'
params.query_sex = '0'
params.score_outdir = 'variant-score'

process RESOLVE_HGVS {
    container params.container_vep
    cpus 1
    memory '4 GB'
    time '30m'

    input:
    val hgvs

    output:
    path 'resolved.vcf', emit: vcf

    script:
    """
    vep --input_data '${hgvs}' --format hgvs --assembly GRCh38 \\
        --cache --merged --dir_cache "${params.VEP_CACHE}" \\
        --fasta "${params.VEP_FASTA}" --synonyms "${params.VEP_SYNONYMS}" \\
        --check_ref --vcf --no_stats --output_file resolved.vcf
    """
}

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

process DECLARE_HMTVAR_INFO {
    cpus 1
    memory '1 GB'
    time '5m'

    input:
    tuple val(group), val(type), path(vcf)

    output:
    tuple val(group), val(type), path('rank-input.vcf'), emit: vcf

    script:
    """
    # The WGS rank model expects the HmtNote header even when an SNV has no HmtVar value.
    awk '
        /^##INFO=<ID=Pathogenicity[,>]/ { declared = 1 }
        /^#CHROM[[:space:]]/ && !declared {
            print "##INFO=<ID=Pathogenicity,Number=1,Type=String,Description=\\"Pathogenicity predicted by HmtVar\\">"
        }
        { print }
    ' "${vcf}" > rank-input.vcf
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
    if ([params.variant, params.hgvs, params.coord].findAll { it != null }.size() != 1) {
        error 'Provide exactly one of --variant (a VCF file), --hgvs (an HGVS SNV), or --coord (CHROM:POS:REF:ALT).'
    }
    if (!params.genome_file) {
        error 'Configure the GRCh38 genome_file (for example with the Hopper wgs profile).'
    }
    def genotype = params.query_genotype.toString()
    def sex = params.query_sex.toString()
    if (!(genotype in ['0/1', '1/1']) || !(sex in ['0', '1', '2'])) {
        error 'Use --query_genotype 0/1 or 1/1 and --query_sex 0, 1 or 2.'
    }

    if (params.coord != null) {
        def coordinate = params.coord.toString().trim()
        if (!(coordinate ==~ /[A-Za-z0-9_.-]+:[1-9][0-9]*:[ACGTacgt]:[ACGTacgt]/)) {
            error 'Use --coord CHROM:POS:REF:ALT with a 1-based GRCh38 position and single-base alleles.'
        }
        def (chrom, pos, ref, alt) = coordinate.split(':')
        def row = [chrom, pos, '.', ref, alt, '.', '.', '.'].join('\t')
        query_vcf = channel.of(row).collectFile(name: 'coordinate.vcf', newLine: true)
    }
    else if (params.hgvs != null) {
        def hgvs = params.hgvs.toString().trim()
        def genomic_snv = /(?:NC_[0-9]+\.[0-9]+|[0-9]+|X|Y|MT):g\.[1-9][0-9]*[ACGT]>[ACGT]/
        def transcript_snv = /(?:ENST[0-9]+|[NX][MR]_[0-9]+)(?:\.[0-9]+)?:[cn]\.[*-]?[1-9][0-9]*(?:[+-][1-9][0-9]*)?[ACGT]>[ACGT]/
        if (!(hgvs ==~ genomic_snv || hgvs ==~ transcript_snv)) {
            error 'Use a genomic or transcript HGVS SNV with a chromosome or transcript reference, such as NC_000011.10:g.5227002T>A or NM_000518.5:c.20A>T.'
        }
        RESOLVE_HGVS(channel.of(hgvs))
        query_vcf = RESOLVE_HGVS.out.vcf
    }
    else {
        query_vcf = channel.of(file(params.variant, checkIfExists: true))
    }

    PREPARE_VARIANT_QUERY(
        query_vcf,
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
    DECLARE_HMTVAR_INFO(inher_models.out.vcf)
    genmodscore(DECLARE_HMTVAR_INFO.out.vcf, false)

    REPORT_VARIANT_SCORE(
        genmodscore.out.scored_vcf.map { group, type, vcf -> vcf },
        file("${projectDir}/bin/prepare_variant_query.py"),
        file("${projectDir}/bin/report_variant_score.py")
    )
}
