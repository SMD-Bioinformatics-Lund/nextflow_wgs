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
