nextflow.enable.dsl=2

workflow VALIDATE_PARAMETERS {

	take:
	params_to_validate

	main:

	log.info "Validating file and directory parameters..."

	params_to_validate.each { key, value ->
		
		if (!params.containsKey(key)) {
			log.info "Parameter '${key}' is listed in params_to_validate but is not defined in params."
		}

		if (!(value instanceof String)) {
			error "ERROR: Parameter '${key}' is not a string path."
		}

		def f = file(value)

		if (!f.exists()) {
			error "ERROR: Param '${key}' points to missing path: ${value}"
		}

		if (f.isFile() && f.size() == 0) {
			error "ERROR: Param '${key}' points to empty file: ${value}"
	}

	}

	log.info "All configured parameter paths validated ✓"

	emit:
	params_to_validate = params_to_validate
}

