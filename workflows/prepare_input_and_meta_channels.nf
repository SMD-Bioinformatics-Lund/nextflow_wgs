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

def readCsvRows(csvPath) {
	def csv = file(csvPath)
	def lines = csv.readLines().findAll { line -> line.trim() }
	def header = lines[0].split(',', -1)*.trim()
	header[0] = header[0].replaceFirst(/^\uFEFF/, '')

	lines.tail().collect { line ->
		def values = line.split(',', -1)*.trim()
		[header, values].transpose().collectEntries()
	}
}

def gatkRefKey(platform, sex) {
	"${platform?.toString()?.trim()?.toLowerCase()}_${sex?.toString()?.trim()?.toUpperCase()}"
}

/*
 * Parse GATK-CNV reference metadata from one profile-defined CSV.
 *
 * Panel profiles should provide one row per sex, even when the same reference
 * paths are used for M and F. Empty pon values are allowed for profiles that do
 * not run the WGS denoised coverage process.
 *
 * CSV content is validated outside the pipeline; validate_params.nf only checks
 * that the CSV path exists and is non-empty.
 */
def gatkRefConfigsFromCsv(gatkRefCsv) {
	def rows = readCsvRows(gatkRefCsv)
	rows.collectEntries { row ->
		def key = gatkRefKey(row.platform, row.sex)
		[key: [
			key: key,
			platform: row.platform,
			sex: row.sex,
			intervals: row.intervals,
			ploidy_model: row.ploidy_model,
			ref_folders: row.ref_folders,
			pon: row.pon
		]]
	}
}

def gatkRefConfigForSample(row, gatkRefConfigs) {
	def key = gatkRefKey(row.platform, row.sex)
	def ref = gatkRefConfigs[key]
	if (!ref) {
		error "No GATK reference found for sample ${row.id} with platform='${row.platform}' and sex='${row.sex}'. Expected key '${key}' in GATK reference CSV. Available keys: ${gatkRefConfigs.keySet().join(', ')}"
	}
	return ref
}

/*
 * Expand each selected reference into GermlineCNVCaller model shards.
 * The ref_folders file is expected to contain columns 'i' and 'refpart'.
 */
def gatkRefShardsFromConfigs(gatkRefConfigs) {
	gatkRefConfigs.collectMany { key, ref ->
		def rows = readCsvRows(ref.ref_folders)
		rows.collect { row ->
			tuple(key, row.i, row.refpart)
		}
	}
}

workflow PREPARE_INPUT_AND_META_CHANNELS {

	take:
	ch_samplesheet
	val_gatk_ref_csv

	main:
	def gatk_ref_configs = val_gatk_ref_csv ? gatkRefConfigsFromCsv(val_gatk_ref_csv) : [:]

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

	ch_gatk_ref_meta = val_gatk_ref_csv ? ch_samplesheet.map { row ->
		def gatk_ref = gatkRefConfigForSample(row, gatk_ref_configs)
		tuple(row.group, row.id, gatk_ref)
	} : channel.empty()

	ch_gatk_ref_shards = val_gatk_ref_csv ? channel.fromList(gatkRefShardsFromConfigs(gatk_ref_configs)) : channel.empty()

	emit:
		sample_meta = ch_sample_meta
		proband_meta = ch_proband_meta
		fastq = ch_fastq_meta
		bam = ch_bam_meta
		vcf = ch_vcf_meta
		gatk_ref_meta = ch_gatk_ref_meta
		gatk_ref_shards = ch_gatk_ref_shards
}
