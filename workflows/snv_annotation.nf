#!/usr/bin/env nextflow

workflow ANNOTATE_SNV_INDELS {


	take:
	ch_snv_indels_vcf
	ch_ped

	main:

	// TODO: Better system so these two do not have to be redefined here.
	params.results_output_dir = params.outdir + '/' + params.subdir
	params.mode = file(params.csv).countLines() > 2 ? "family" : "single"

	ch_versions = Channel.empty()
	ch_output_info = Channel.empty()

	annotate_vep(ch_snv_indels_vcf)
	vcfanno(annotate_vep.out.vcf)
	modify_vcf(vcfanno.out.vcf)
	mark_splice(modify_vcf.out.vcf)

	// INDELS //
	extract_indels_for_cadd(ch_snv_indels_vcf)
	indel_vep(extract_indels_for_cadd.out.vcf)
	calculate_indel_cadd(indel_vep.out.vcf)
	bgzip_indel_cadd(calculate_indel_cadd.out.cadd_gz)


	add_cadd_scores_to_vcf(mark_splice.out.splice_marked.join(bgzip_indel_cadd.out.cadd_tbi))


	// INHERITANCE MODELS //
	ch_inher_models_input = add_cadd_scores_to_vcf.out.vcf
		.cross(ch_ped)
		.map { vcf_tuple, ped_tuple ->
				def group = vcf_tuple[0]
				def vcf = vcf_tuple[1]
				def type = ped_tuple[1]
				def ped = ped_tuple[2]
				tuple(group, vcf, type, ped)
		}

	inher_models(ch_inher_models_input)

	// SCORE VARIANTS //
	genmodscore(inher_models.out.vcf)
	vcf_completion(genmodscore.out.scored_vcf)
	ch_output_info = ch_output_info.mix(vcf_completion.out.snv_INFO)

	// VERSIONS
	ch_versions = ch_versions.mix(annotate_vep.out.versions.first())
	ch_versions = ch_versions.mix(vcfanno.out.versions.first())
	ch_versions = ch_versions.mix(extract_indels_for_cadd.out.versions.first())
	ch_versions = ch_versions.mix(indel_vep.out.versions.first())
	ch_versions = ch_versions.mix(calculate_indel_cadd.out.versions.first())
	ch_versions = ch_versions.mix(bgzip_indel_cadd.out.versions.first())
	ch_versions = ch_versions.mix(add_cadd_scores_to_vcf.out.versions.first())


	emit:
	annotated_snv_vcf = vcf_completion.out.vcf_tbi
	versions = ch_versions
	output_info = ch_output_info

}


process annotate_vep {
	container  "${params.container_vep}"
	cpus 30
	tag "$group"
	memory '50 GB'
	time '5h'

	input:
		tuple val(group), val(id), path(vcf)

	output:
		tuple val(group), path("${group}.vep.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		vep \\
			-i ${vcf} \\
			-o ${group}.vep.vcf \\
			--offline \\
			--sift b --polyphen b --ccds --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --af --af_1kg --max_af --pubmed --uniprot --mane --tsl --appris --variant_class --gene_phenotype --mirna \\
			--merged \\
			--vcf \\
			--no_stats \\
			--synonyms $params.VEP_SYNONYMS \\
			--fork ${task.cpus} \\
			--force_overwrite \\
			--fasta $params.VEP_FASTA \\
			--dir_cache $params.VEP_CACHE \\
			--dir_plugins $params.VEP_PLUGINS \\
			--distance $params.VEP_TRANSCRIPT_DISTANCE \\
			-cache \\
			--plugin CADD,$params.CADD \\
			--plugin LoFtool \\
			--plugin MaxEntScan,$params.MAXENTSCAN,SWA,NCSS \\
			--plugin dbNSFP,$params.DBNSFP,transcript_match=1,REVEL_score,REVEL_rankscore,BayesDel_addAF_score,BayesDel_addAF_pred,BayesDel_noAF_score,BayesDel_noAF_pred,GERP++_RS \\
			-custom $params.GNOMAD_EXOMES,gnomADe,vcf,exact,0,AF_grpmax,AF,grpmax \\
			-custom $params.GNOMAD_GENOMES,gnomADg,vcf,exact,0,AF_grpmax,AF,grpmax \\
			-custom $params.GNOMAD_MT,gnomAD_mt,vcf,exact,0,AF_hom,AF_het \\
			-custom $params.PHYLOP,phyloP100way,bigwig \\
			-custom $params.PHASTCONS,phastCons,bigwig


		${annotate_vep_version(task)}
		"""

	stub:
		"""
		touch "${group}.vep.vcf"
		${annotate_vep_version(task)}
		"""
}
def annotate_vep_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
	END_VERSIONS
	"""
}


// gene, clinvar, loqusdb, enigma(onco)
process vcfanno {
	memory '1GB'
	time '1h'
	errorStrategy 'retry'
	maxErrors 5
	cpus 2

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.clinvar.loqusdb.gene.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		vcfanno_linux64 -lua $params.VCFANNO_LUA $params.vcfanno $vcf > ${group}.clinvar.loqusdb.gene.vcf
		${vcfanno_version(task)}
		"""

	stub:
		"""
		touch "${group}.clinvar.loqusdb.gene.vcf"
		${vcfanno_version(task)}
		"""
}
def vcfanno_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcfanno: \$(echo \$(vcfanno_linux64 2>&1 | grep version | cut -f3 -d' ')  )
	END_VERSIONS
	"""
}


// Extracting most severe consequence:
// Modifying annotations by VEP-plugins, and adding to info-field:
// Modifying CLNSIG field to allow it to be used by genmod score properly:
// TODO: give process better name
process modify_vcf {
	cpus 2
	tag "$group"
	memory '1 GB'
	time '1h'

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.mod.vcf"), emit: vcf

	script:
		"""
		modify_vcf_scout.pl $vcf > ${group}.mod.vcf
		"""

	stub:
		"""
		touch "${group}.mod.vcf"
		"""
}


// Marking splice INDELs:
process mark_splice {
	cpus 2
	tag "$group"
	memory '1 GB'
	time '1h'

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.marksplice.vcf"), emit: splice_marked

	script:
		"""
		/opt/bin/mark_spliceindels.pl $vcf > ${group}.marksplice.vcf
		"""

	stub:
		"""
		touch "${group}.marksplice.vcf"
		"""
}

// Extract all INDELs
process extract_indels_for_cadd {
	cpus 2
	tag "$group"
	memory '1 GB'
	time '1h'

	input:
		tuple val(group), val(id), path(vcf)

	output:
		tuple val(group), path("${group}.only_indels.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		bcftools view $vcf -V snps -o ${group}.only_indels.vcf
		${extract_indels_for_cadd_version(task)}
		"""

	stub:
		"""
		touch "${group}.only_indels.vcf"
		${extract_indels_for_cadd_version(task)}
		"""
}
def extract_indels_for_cadd_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS
	"""
}

// Annotate Indels with VEP+Gnomad genomes. Filter variants below threshold
process indel_vep {
	cpus 5
	container  "${params.container_vep}"
	tag "$group"
	memory '10 GB'
	time '3h'

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.only_indels.vep.filtered.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		vep \\
			-i $vcf \\
			-o ${group}.only_indels.vep.vcf \\
			--offline \\
			--cache \\
			--merged \\
			--vcf \\
			--synonyms $params.VEP_SYNONYMS \\
			--fasta $params.VEP_FASTA \\
			-custom $params.GNOMAD_GENOMES,gnomADg,vcf,exact,0,AF \\
			-custom $params.GNOMAD_MT,gnomAD_mt,vcf,exact,0,AF_hom,AF_het \\
			--dir_cache $params.VEP_CACHE \\
			--force_overwrite \\
			--no_stats \\
			--fork ${task.cpus}
		filter_indels.pl ${group}.only_indels.vep.vcf > ${group}.only_indels.vep.filtered.vcf
		${indel_vep_version(task)}
		"""

	stub:
		"""
		touch "${group}.only_indels.vep.filtered.vcf"
		${indel_vep_version(task)}
		"""
}
def indel_vep_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
	END_VERSIONS
	"""
}

// Calculate CADD scores for all indels
process calculate_indel_cadd {
	cpus 2
	container  "${params.container_cadd}"
	tag "$group"
	memory '20 GB'
	time '3h'

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.indel_cadd.gz"), emit: cadd_gz
		path "*versions.yml", emit: versions

	script:
		"""
		/CADD-scripts/CADD.sh -c ${task.cpus} -g GRCh38 -o ${group}.indel_cadd.gz $vcf
		${calculate_indel_cadd_version(task)}
		"""

	stub:
		"""
		touch "${group}.indel_cadd.gz"
		${calculate_indel_cadd_version(task)}
		"""
}
def calculate_indel_cadd_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    cadd: \$(echo \$(/CADD-scripts/CADD.sh -v 2>&1) | sed -e "s/^.*CADD-v// ; s/ (c).*//")
	END_VERSIONS
	"""
}

process bgzip_indel_cadd {

	tag "$group"
	cpus 4
	memory '1 GB'
	time '5m'
	container  "${params.container_bcftools}"

	input:
		tuple val(group), path(cadd_scores)

	output:
		tuple val(group), path("${group}.cadd.gz"), path("${group}.cadd.gz.tbi"), emit: cadd_tbi
		path "*versions.yml", emit: versions

	script:
		"""
		gunzip -c ${cadd_scores} > "${group}.cadd"
		bgzip -@ ${task.cpus} "${group}.cadd"
		tabix -p vcf "${group}.cadd.gz"
		${bgzip_indel_cadd_version(task)}
		"""

	stub:
		"""
		touch "${group}.cadd.gz"
		touch "${group}.cadd.gz.tbi"
		${bgzip_indel_cadd_version(task)}
		"""
}
def bgzip_indel_cadd_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS
	"""
}

// Add the calculated indel CADDs to the vcf
process add_cadd_scores_to_vcf {
	cpus 4
	tag "$group"
	memory '1 GB'
	time '5m'
	container  "${params.container_genmod}"

	input:
		tuple val(group), path(vcf), path(cadd_scores), path(cadd_scores_tbi)

	output:
		tuple val(group), path("${group}.cadd.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		genmod annotate --cadd-file ${cadd_scores} ${vcf} > ${group}.cadd.vcf

		${add_cadd_scores_to_vcf_version(task)}
		"""

	stub:
		"""
		touch "${group}.cadd.vcf"
		${add_cadd_scores_to_vcf_version(task)}
		"""
}
def add_cadd_scores_to_vcf_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    genmod: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
	END_VERSIONS
	"""
}


// # Annotating variant inheritance models:
process inher_models {
	tag "$group"
	cpus 3
	memory '80 GB'
	time '1h'
	container  "${params.container_genmod}"

	input:
		tuple val(group), path(vcf), val(type), path(ped)

	output:
		tuple val(group), val(type), path("${group}.models.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		genmod models $vcf -p ${task.cpus} -f $ped > ${group}.models.vcf
		${inher_models_version(task)}
		"""

	stub:
		"""
		touch "${group}.models.vcf"
		${inher_models_version(task)}
		"""
}
def inher_models_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    genmod: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
	END_VERSIONS
	"""
}


// Scoring variants:
// Adjusting compound scores:
// Sorting VCF according to score:
process genmodscore {
	tag "$group"
	cpus 2
	memory '20 GB'
	time '1h'
	container  "${params.container_genmod}"

	input:
		tuple val(group), val(type), path(vcf)

	output:
		tuple val(group), val(type), path("${group_score}.scored.vcf"), emit: scored_vcf
		path "*versions.yml", emit: versions

	script:
		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group

		if ( params.mode == "family" && params.antype == "wgs" ) {
			"""
			genmod score -i $group_score -c $params.rank_model -r $vcf -o ${group_score}.only_rankscore.vcf
			genmod compound \
				--threshold ${params.genmod_compound_trio_threshold} \
				--penalty ${params.genmod_compound_trio_penalty} \
				-o ${group_score}.with_compounds.vcf \
				${group_score}.only_rankscore.vcf
			sed 's/RankScore=${group}:/RankScore=${group_score}:/g' -i ${group_score}.with_compounds.vcf
			genmod sort -p -f $group_score ${group_score}.with_compounds.vcf -o ${group_score}.scored.vcf

			${genmodscore_version(task)}
			"""
		}
		else {
			"""
			genmod score -i $group_score -c $params.rank_model_s -r $vcf -o ${group_score}.only_rankscore.vcf

			# To get compounds without applying rank score penalty
			genmod compound \
				--penalty 0 \
				-o ${group_score}.with_compounds.vcf \
				${group_score}.only_rankscore.vcf

			genmod sort \
				-p \
				-f $group_score \
				-o ${group_score}.scored.vcf \
				${group_score}.with_compounds.vcf

			${genmodscore_version(task)}
			"""
		}

	stub:
		group_score = group
		"""
		touch "${group_score}.scored.vcf"
		${genmodscore_version(task)}
		"""
}
def genmodscore_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    genmod: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
	END_VERSIONS
	"""
}

// Bgzipping and indexing VCF:
process vcf_completion {
	cpus 16
	publishDir "${params.results_output_dir}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
	tag "$group"
	time '1h'
	memory '5 GB'

	input:
		tuple val(group), val(type), path(vcf)

	output:
		tuple val(group), val(type), path("${group_score}.scored.vcf.gz"), path("${group_score}.scored.vcf.gz.tbi"), emit: vcf_tbi
		tuple val(group), path("${group}_snv.INFO"), emit: snv_INFO
		path "*versions.yml", emit: versions

	script:
		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group

		"""
		sed 's/^MT/M/' -i $vcf
		sed 's/ID=MT,length/ID=M,length/' -i $vcf
		bgzip -@ ${task.cpus} $vcf -f
		tabix ${vcf}.gz -f
		echo "SNV	$type	${params.accessdir}/vcf/${group_score}.scored.vcf.gz" > ${group}_snv.INFO

		${vcf_completion_version(task)}
		"""

	stub:
		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group

	"""
		touch "${group_score}.scored.vcf.gz"
		touch "${group_score}.scored.vcf.gz.tbi"
		touch "${group}_snv.INFO"

		${vcf_completion_version(task)}
		"""
}
def vcf_completion_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	END_VERSIONS
	"""
}
