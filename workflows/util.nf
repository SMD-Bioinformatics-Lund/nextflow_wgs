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

def gatkDefaultRefConfig() {
	[
		intervals: params.gatk_intervals,
		ploidy_model: params.ploidymodel,
		ref_folders: params.gatkreffolders,
		pon: [
			F: params.gatk_pon_female,
			M: params.gatk_pon_male
		]
	]
}

def gatkRefConfigs() {
	if (params.containsKey('gatk_ref_map')) {
		return params.gatk_ref_map
	}

	return ["default": gatkDefaultRefConfig()]
}

def gatkRefKey(meta) {
	def refs = gatkRefConfigs()
	def platform = meta.platform?.toString()
	def sex = meta.sex?.toString()
	def candidates = [
		platform && sex ? "${platform}_${sex}" : null,
		platform,
		sex,
		'default'
	].findAll { it }

	def key = candidates.find { refs.containsKey(it) }
	if (!key) {
		key = candidates.collect { candidate ->
			refs.keySet().find { refKey -> refKey.toString().equalsIgnoreCase(candidate) }
		}.find { it }
	}
	if (!key) {
		error "No GATK reference config found for sample ${meta.id} with platform='${platform}' and sex='${sex}'. Tried keys: ${candidates.join(', ')}"
	}

	return key
}

def gatkRefConfig(meta) {
	def key = gatkRefKey(meta)
	def ref = gatkRefConfigs()[key]
	return ref + [key: key]
}

def gatkRefShards() {
	gatkRefConfigs().collectMany { key, ref ->
		['intervals', 'ploidy_model', 'ref_folders'].each { field ->
			if (!ref.containsKey(field)) {
				error "GATK reference config '${key}' is missing '${field}'"
			}
		}

		def refCsv = file(ref.ref_folders)
		def lines = refCsv.readLines().findAll { it.trim() }
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
	publishDir "${params.outdir}/${params.subdir}/vcf", mode: 'copy', overwrite: true, pattern: '*.vcf.gz*'
	tag "$group"
	time '1h'
	memory '5 GB'

	input:
		tuple val(group), val(type), path(vcf)

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
