nextflow.enable.dsl=2

workflow VALIDATE_PARAMETERS {

	take:
	global_parameters // params_to_validate

	main:

	log.info "Validating file and directory parameters..."

    params.each { key, value ->

        if (!global_parameters.contains(key))
            return

        if (!(value instanceof String)) {
            error "ERROR: Parameter '${key}' is not a string path."
        }

        if (!value.startsWith('/'))
            return

        if (value == null)
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
	global_parameters = global_parameters
}

