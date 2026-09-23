workflow PEDDY_QC {
    take:
    ch_annotated_vcf        // [group, type, vcf, index]
    ch_ped                  // [group, type, ped]
    ch_samplesheet          // Samplesheet rows, including group, id and sequencing_run.
    val_results_output_dir
    val_accessdir
    val_cdm_assay

    main:
    ch_peddy_input = ch_annotated_vcf
        .filter { group, type, vcf, index -> type == "proband" }
        .join(ch_ped, by: [0, 1])

    peddy(ch_peddy_input, val_accessdir)

    ch_sample_runs = ch_samplesheet
        .map { row -> tuple(row.group, row.id, row.sequencing_run) }
        .groupTuple()

    peddy2cdm(
        peddy.out.peddy_files.join(ch_sample_runs, by: [0]),
        val_results_output_dir,
        val_cdm_assay
    )

    emit:
    peddy_files = peddy.out.peddy_files // [group, ped_check.csv, peddy.ped, sex_check.csv]
    output_info = peddy.out.peddy_INFO
    versions = peddy.out.versions.first()
    json = peddy2cdm.out.json
    cdm = peddy2cdm.out.cdm
}

process peddy {

	publishDir "${params.outdir}/${params.subdir}/ped", mode: 'copy' , overwrite: true, pattern: '*.ped'
	publishDir "${params.outdir}/${params.subdir}/ped", mode: 'copy' , overwrite: true, pattern: '*.csv'

	cpus 4
	tag "$group"
	time '1h'
	memory '20GB'

	input:
		tuple val(group), val(type), path(vcf), path(idx), path(ped)
		val accessdir

	output:
		tuple val(group), path("${group}.ped_check.csv"),path("${group}.peddy.ped"), path("${group}.sex_check.csv"), emit: peddy_files
		tuple val(group), path("${group}_peddy.INFO"), emit: peddy_INFO
		path "*versions.yml", emit: versions

	script:
		"""
		source activate py3-env
		python -m peddy --sites hg38 -p ${task.cpus} $vcf $ped --prefix $group
		echo "PEDDY	${accessdir}/ped/${group}.ped_check.csv,${accessdir}/ped/${group}.peddy.ped,${accessdir}/ped/${group}.sex_check.csv" > ${group}_peddy.INFO

		${peddy_version(task)}
		"""

	stub:
		"""
		source activate py3-env
		touch "${group}.ped_check.csv"
		touch "${group}.peddy.ped"
		touch "${group}.sex_check.csv"
		touch "${group}_peddy.INFO"

		${peddy_version(task)}
		"""
}
def peddy_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    peddy: \$(echo \$(python -m peddy --version 2>&1) | sed 's/^.*peddy, version //')
	END_VERSIONS
	"""
}

process peddy2cdm {
	cpus 2
	memory '20 MB'
	tag "$group"
	publishDir "${params.outdir}/${params.subdir}/qc", mode: 'copy', overwrite: true, pattern: '*.json'
	publishDir "${params.crondir}/peddy", mode: 'copy' , overwrite: true, pattern: '*.peddy2cdm'
	container "${params.container_pysam_cmdvcf}"
	time '20m'

	input:
		tuple val(group), path(ped_check),path(peddy_ped), path(sex_check), val(id), val(sequencing_run)
		val results_output_dir
		val cdm_assay

	output:
		tuple val(group), path("*peddy.json"), emit: json
		tuple val(group), path("*peddy2cdm"), emit: cdm

	script:
		def sample_arg = [id, sequencing_run]
			.transpose()
			.collect { sample_id, run_id -> "${sample_id}:${run_id}" }
			.join(' --sample ')
			
		"""
		peddy2cdm.py \
		--ped $ped_check \
		--sex $sex_check \
		--sample $sample_arg \
		--cdmassay $cdm_assay \
		--results_dir ${results_output_dir}/qc
		"""
		

	stub:
		"""
		touch "${group}_peddy.json"
		touch "${group}.peddy2cdm"
	    """

}
