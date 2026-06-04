nextflow.enable.dsl=2

def sampleMeta(row) {
	[
		group: row.group,
		id: row.id,
		type: row.type,
		sex: row.sex,
		mother: row.mother,
		father: row.father,
		diagnosis: row.diagnosis,
		phenotype: row.phenotype,
		assay: row.assay,
		platform: row.platform,
		clarity_sample_id: row.clarity_sample_id,
		sequencing_run: row.sequencing_run,
		n_reads: row.containsKey("n_reads") ? row.n_reads : null,
		analysis: row.containsKey("analysis") ? row.analysis : false,
		ffpe: row.containsKey("ffpe") ? row.ffpe : false,
		priority: row.containsKey("priority") ? row.priority : null
	]
}

def caseMetaFromProband(meta) {
	[
		group: meta.group,
		proband_id: meta.id,
		proband_sex: meta.sex,
		diagnosis: meta.diagnosis,
		phenotype: meta.phenotype,
		assay: meta.assay,
		clarity_sample_id: meta.clarity_sample_id,
		analysis: meta.analysis,
		scout_case_status: meta.priority == "highest" ? "prioritized" : ""
	]
}

workflow PREPARE_META_CHANNELS {

	take:
	ch_samplesheet

	main:

	ch_sample_meta = ch_samplesheet.map { row ->
		tuple(sampleMeta(row))
	}

	ch_proband_meta = ch_samplesheet
		.filter { row ->
			row.type == "proband"
		}
		.map { row ->
			tuple(sampleMeta(row))
		}

	ch_case_meta = ch_samplesheet
		.filter { row ->
			row.type == "proband"
		}
		.map { row ->
			tuple(caseMetaFromProband(sampleMeta(row)))
		}

	ch_fastq_meta = ch_samplesheet
		.filter { row ->
			row.read1.endsWith("q.gz") && row.read2.endsWith("q.gz")
		}
		.map { row ->
			tuple(sampleMeta(row), row.read1, row.read2)
		}

	ch_bam_meta = ch_samplesheet
		.filter { row ->
			row.read1.endsWith("bam")
		}
		.map { row ->
			tuple(sampleMeta(row), row.read1, row.read2)
		}

	ch_vcf_meta = ch_samplesheet
		.filter { row ->
			row.read1.endsWith(".vcf") || row.read1.endsWith(".vcf.gz")
		}
		.map { row ->
			tuple(sampleMeta(row), row.read1)
		}

	emit:
		sample_meta = ch_sample_meta
		proband_meta = ch_proband_meta
		case_meta = ch_case_meta
		fastq = ch_fastq_meta
		bam = ch_bam_meta
		vcf = ch_vcf_meta
}