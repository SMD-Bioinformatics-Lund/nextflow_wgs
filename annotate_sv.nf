#!/usr/bin/env nextflow

include { vep_sv             } from './modules/vep_sv.nf'
include { gnomad_sv_bcftools } from './modules/gnomad_sv_bcftools.nf'
include { bgzip_index_vcf    } from './workflows/util.nf'

nextflow.enable.dsl = 2

params.vcf = null
params.sample_id = null
params.outdir = 'results'

workflow {
    if (!params.vcf || !params.sample_id) {
        error 'Provide --vcf and --sample_id to annotate an SV VCF'
    }

    def sample_id = params.sample_id.toString()
    if (!(sample_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
        error 'Sample ID must start with a letter or digit and contain only letters, digits, dots, underscores, or hyphens'
    }

    ch_input = channel.of(tuple(sample_id, sample_id, file(params.vcf, checkIfExists: true)))
    vep_sv(ch_input)
    gnomad_sv_bcftools(
        vep_sv.out.vep_sv_vcf,
        file(params.GNOMAD_SV_INTERVALS, checkIfExists: true),
        file(params.VEP_SYNONYMS, checkIfExists: true)
    )
    bgzip_index_vcf(gnomad_sv_bcftools.out.vcf, "${params.outdir}/vcf")
}
