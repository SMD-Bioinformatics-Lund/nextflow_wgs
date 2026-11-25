process vep_sv {
	cpus 10
	container  "${params.container_vep}"
	tag "$group"
	memory '50 GB'
	time '1h'

	input:
		tuple val(group), val(id), path(vcf)

	output:
		tuple val(group), val(id), path("${group}.vep.vcf"), emit: vep_sv_vcf
		path "*versions.yml", emit: versions

	script:
		"""
		vep \\
			-i ${vcf} \\
			-o ${group}.vep.vcf \\
			--offline \\
			--merged \\
			--sift b --polyphen b --ccds --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --af --af_1kg --max_af --pubmed --uniprot --mane --tsl --appris --variant_class --gene_phenotype --mirna \\
			--synonyms $params.VEP_SYNONYMS \\
			--vcf \\
			--no_stats \\
			--fork ${task.cpus} \\
			--force_overwrite \\
			--plugin LoFtool \\
			--fasta $params.VEP_FASTA \\
			--dir_cache $params.VEP_CACHE \\
			--dir_plugins $params.VEP_PLUGINS \\
			--max_sv_size $params.VEP_MAX_SV_SIZE \\
			--distance $params.VEP_TRANSCRIPT_DISTANCE \\
			-cache \\
			--format vcf

		${vep_sv_version(task)}
		"""

	stub:
		"""
		touch "${group}.vep.vcf"
		${vep_sv_version(task)}
		"""
}
def vep_sv_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
	END_VERSIONS
	"""
}
