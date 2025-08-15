def bcftools_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS
	"""
}

process IDSNP_CALL {
    label 'process_single'
    tag "${id}"
    container "${params.container_bcftools}"

    input:
        tuple val(group), val(id), path(bam), path(bai)
        val idsnp_params

    output:
        tuple val(group), val(id), path("*final.vcf"), emit: vcf
        path("*versions.yml"), emit: versions

    script:
        def prefix  = "${id}"
        def max_depth = 1000
        def min_map_qual = 10
        """
        bcftools mpileup \\
            -Ou \\
            -R "${idsnp_params.idsnp_bed}" \\
            -f "${idsnp_params.genome}" \\
            -d "${max_depth}" \\
            -q "${min_map_qual}" \\
            "${bam}" | \\
        bcftools call \\
            -A \\
            -C alleles \\
            -T "${idsnp_params.idsnp_tsv}" \\
            -m \\
            -Ov \\
            > "${prefix}.raw.vcf"
        
        bcftools annotate \\
            -a "${idsnp_params.variant_rsids_bed}" \\
            -c "CHROM,FROM,TO,ID" \\
            -h "${idsnp_params.header}" \\
            -o "${prefix}.final.vcf" "${prefix}.raw.vcf"

        ${bcftools_version(task)}
        """

    stub:
        def prefix = "${id}"
        """
        touch ${prefix}.final.vcf

        ${bcftools_version(task)}
        """
}

process IDSNP_VCF_TO_JSON {
    label 'process_single'
    tag "${id}"
    container "${params.container_python}"
	publishDir "${params.outdir}/${params.subdir}/qc", mode: 'copy' , overwrite: 'true', pattern: '*.json'

    input:
        tuple val(group), val(id), path(vcf)
    
    output:
        tuple val(group), val(id), path("*.json"), emit: json
    
    script:
    def prefix = "${id}"
    """
    genotype_to_json.py "${vcf}" "${prefix}.genotypes.json"
    """

    stub:
    def prefix = "${id}"
    """
    touch "${prefix}.genotypes.json"
    """
}
