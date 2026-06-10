workflow QC_TO_CDM {

	take:
	ch_qc_json             // ch:     [val(group), val(id), path(qc_jsons)]
	ch_sample_meta         // ch:     [val(group), val(id), val(meta)]
	val_results_output_dir // string: Full result base directory under which pipeline results are published.
	val_cdm_assay          // string: CDM assay name used when creating QC cron files.
	val_skip_cdm_cron      // bool:   Whether to skip creating CDM QC cron files.

	main:
	merge_qc_json(ch_qc_json)
	ch_cdm_done = channel.empty()

	if (!val_skip_cdm_cron) {
		qc_to_cdm(
			ch_sample_meta.join(merge_qc_json.out.qc_cdm_merged, by: [0,1]),
			val_results_output_dir,
			val_cdm_assay
		)
		ch_cdm_done = qc_to_cdm.out.cdm_done
	}

	emit:
	cdm_done = ch_cdm_done
}


/////////////// Collect QC, emit: single file ///////////////

process merge_qc_json {
    cpus 2
    errorStrategy 'retry'
    maxErrors 5
    publishDir "${params.outdir}/${params.subdir}/qc", mode: 'copy' , overwrite: true, pattern: '*.QC'
    tag "$id"
    time '1h'
	memory '1 GB'

    input:
        tuple val(group), val(id), path(qc)

    output:
    	tuple val(group), val(id), path("${id}.QC"), emit: qc_cdm_merged

    script:
        qc_json_files = qc.join(' ')
		"""
		merge_json_files.py ${qc_json_files} > ${id}.QC
		"""

	stub:
		"""
		touch "${id}.QC"
		"""
}

// Load QC data, emit: CDM (via middleman)
process qc_to_cdm {
	cpus 2
	errorStrategy 'retry'
	maxErrors 5
	publishDir "${params.crondir}/qc", mode: 'copy' , overwrite: true
	tag "$id"
	time '1h'

	input:
		tuple val(group), val(id), val(meta), path(qc_json)
		val(results_output_dir)
		val(cdm_assay)

	output:
		path("${id}.cdmpy"), emit: cdm_done

	script:
		"""
		echo "--sequencing-run ${meta.sequencing_run} --sample-id ${id} --assay ${cdm_assay} --subassay ${meta.diagnosis} --qc ${results_output_dir}/qc/${id}.QC --lims-id ${meta.clarity_sample_id}" > ${id}.cdmpy
		"""
}
