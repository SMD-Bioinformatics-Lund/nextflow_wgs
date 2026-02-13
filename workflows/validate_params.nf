nextflow.enable.dsl=2

workflow VALIDATE_PARAMETERS {

	take:
	params_to_validate

	main:

	log.info "Validating file and directory parameters..."
	// Check that all required validation keys exist in params
	params_to_validate.unique().each { key ->
		if (!params.containsKey(key)) {
			log.info "Parameter '${key}' is listed in params_to_validate but is not defined in params."
		}
	}
	params.each { key, value ->
		
		// iterate over all params and only check defined ones
		// profiles may have unique params and this makes it easier
		if (!params_to_validate.contains(key))
			return

		if (!(value instanceof String)) {
			error "ERROR: Parameter '${key}' is not a string path."
		}

		// some params are set to 'PH' in some profiles while being a file path in others
		if (value.startsWith('PH'))
			return

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

