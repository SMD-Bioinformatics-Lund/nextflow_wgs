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
		priority: row.containsKey("priority") ? (row.priority == "highest" ? "prioritized" : "") : ""
	]
}

workflow PREPARE_INPUT_AND_META_CHANNELS {

	take:
	ch_samplesheet

	main:

	ch_sample_meta = ch_samplesheet.map { row ->
		def meta = sampleMeta(row)
		tuple(meta.group, meta.id, meta)
	}

	ch_proband_meta = ch_samplesheet
		.filter { row ->
			row.type == "proband"
		}
		.map { row ->
			def meta = sampleMeta(row)
			tuple ( meta.group, meta.id, meta)
		}

	ch_fastq_meta = ch_samplesheet
		.filter { row ->
			row.read1.endsWith("q.gz") && row.read2.endsWith("q.gz")
		}
		.map { row ->
			def meta = sampleMeta(row)
			tuple(meta.group, meta.id, row.read1, row.read2)
		}

	ch_bam_meta = ch_samplesheet
		.filter { row ->
			row.read1.endsWith("bam")
		}
		.map { row ->
			def meta = sampleMeta(row)
			tuple(meta.group, meta.id, row.read1, row.read2)
		}

	ch_vcf_meta = ch_samplesheet
		.filter { row ->
			row.read1.endsWith(".vcf") || row.read1.endsWith(".vcf.gz")
		}
		.map { row ->
			def meta = sampleMeta(row)
			tuple(meta.group, row.read1)
		}

	emit:
		sample_meta = ch_sample_meta
		proband_meta = ch_proband_meta
		fastq = ch_fastq_meta
		bam = ch_bam_meta
		vcf = ch_vcf_meta
}