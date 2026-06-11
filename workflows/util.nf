nextflow.enable.dsl = 2


/*
 * Returns true if VCF contains at least one non-header, non-blank line 
 * Returns false if VCF does not exist OR contains NO non-header records 
 */
def vcfHasVariants(Path vcf) {
    
    if (vcf.size() < 1) {
        return false
    }

    def found = false

    // Detect  gzipped files
    def inputStream = java.nio.file.Files.newInputStream(vcf)
    if (vcf.name.endsWith('.gz')) {
        inputStream = new java.util.zip.GZIPInputStream(inputStream)
    }

    inputStream.withReader { reader ->
        reader.eachLine { line ->
            def s = line.trim()
            if (s && !s.startsWith('#')) {
                found = true
                return
            }
        }
    }

    return found
}

/*
 * Build a single GATK-CNV reference config from the legacy flat params.
 *
 * Panel profiles such as onco and constitutional only have one GATK-CNV
 * reference setup, so their config defines params.gatk_intervals,
 * params.ploidymodel and params.gatkreffolders directly. gatkRefConfigs()
 * wraps this into the same map shape used by profiles with multiple
 * reference sets.
 */
def gatkDefaultRefConfig(workflowParams) {
	[
		intervals: workflowParams.gatk_intervals,
		ploidy_model: workflowParams.ploidymodel,
		ref_folders: workflowParams.gatkreffolders,
		pon: [
			F: workflowParams.gatk_pon_female,
			M: workflowParams.gatk_pon_male
		]
	]
}

/*
 * Return all available GATK-CNV reference configs.
 *
 * WGS defines params.gatk_ref_map because it can select different references
 * from sample metadata, for example illumina vs illuminax. Profiles without
 * params.gatk_ref_map fall back to one default config keyed as 'illumina'.
 */
def gatkRefConfigs(workflowParams) {
	if (workflowParams.containsKey('gatk_ref_map')) {
		return workflowParams.gatk_ref_map
	}

	return [illumina: gatkDefaultRefConfig(workflowParams)]
}

/*
 * Choose which GATK-CNV reference config a sample should use.
 *
 * The lookup prefers explicit metadata-derived keys, then falls back to the
 * default 'illumina' key. This lets WGS route samples by platform while keeping
 * single-reference panel profiles compatible with the same downstream flow.
 */
def gatkRefKey(meta, workflowParams) {
	def refs = gatkRefConfigs(workflowParams)
	def platform = meta.platform?.toString()
	def sex = meta.sex?.toString()
	def candidates = [
		platform && sex ? "${platform}_${sex}" : null,
		platform,
		sex,
		'illumina'
	].findAll { candidate -> candidate }

	def key = candidates.find { candidate -> refs.containsKey(candidate) }
	if (!key) {
		key = candidates.collect { candidate ->
			refs.keySet().find { refKey -> refKey.toString().equalsIgnoreCase(candidate) }
		}.find { matchedKey -> matchedKey }
	}
	if (!key) {
		error "No GATK reference config found for sample ${meta.id} with platform='${platform}' and sex='${sex}'. Tried keys: ${candidates.join(', ')}"
	}

	return key
}

/*
 * Return the selected GATK-CNV reference config and include the chosen key.
 * The key is carried through channels so sharded CNV calls and postprocessing
 * are joined only against reference shards from the same config.
 */
def gatkRefConfig(meta, workflowParams) {
	def key = gatkRefKey(meta, workflowParams)
	def ref = gatkRefConfigs(workflowParams)[key]
	return ref + [key: key]
}

/*
 * Expand each GATK-CNV reference config into model shards.
 *
 * params.gatkreffolders / ref_folders points to a small CSV-like file with at
 * least columns 'i' and 'refpart'. Each row becomes one GermlineCNVCaller
 * shard: 'i' is the shard label used in output names, and 'refpart' is the
 * model shard folder passed to GATK with --model.
 */
def gatkRefShards(workflowParams) {
	gatkRefConfigs(workflowParams).collectMany { key, ref ->
		['intervals', 'ploidy_model', 'ref_folders'].each { field ->
			if (!ref.containsKey(field)) {
				error "GATK reference config '${key}' is missing '${field}'"
			}
		}

		def refCsv = file(ref.ref_folders)
		def lines = refCsv.readLines().findAll { line -> line.trim() }
		if (!lines) {
			error "GATK reference shard file is empty for config '${key}': ${ref.ref_folders}"
		}

		def header = lines[0].split(',', -1)*.trim()
		def iIdx = header.indexOf('i')
		def refpartIdx = header.indexOf('refpart')
		if (iIdx < 0 || refpartIdx < 0) {
			error "GATK reference shard file for config '${key}' must contain columns 'i' and 'refpart': ${ref.ref_folders}"
		}

		lines.tail().collect { line ->
			def fields = line.split(',', -1)*.trim()
			tuple(key, fields[iIdx], fields[refpartIdx])
		}
	}
}

process bgzip_index_vcf {
	cpus 16
	publishDir "${publish_dir}", mode: 'copy', overwrite: true, pattern: '*.vcf.gz*'
	tag "$group"
	time '1h'
	memory '5 GB'

	input:
		tuple val(group), val(type), path(vcf)
		val publish_dir

	output:
		tuple val(group), val(type), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: compressed_indexed_vcf
		path "*versions.yml", emit: versions

	script:
		"""
		bgzip -@ ${task.cpus} $vcf -f
		tabix ${vcf}.gz -f

        ${bgzip_index_vcf_version(task)}
		"""

	stub:
	"""
		touch "${group}.vcf.gz"
		touch "${group}.vcf.gz.tbi"

        ${bgzip_index_vcf_version(task)}
		"""
}
def bgzip_index_vcf_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	END_VERSIONS
	"""
}
