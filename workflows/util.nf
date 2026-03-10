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

process bgzip_index_vcf {
	cpus 16
	publishDir "${params.results_output_dir}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
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
