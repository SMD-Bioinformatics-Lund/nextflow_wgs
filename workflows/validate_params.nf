nextflow.enable.dsl=2

workflow VALIDATE_PARAMETERS {

	take:
	global_parameters // params_to_validate

	main:

	log.info "Validating file and directory parameters..."

    global_parameters.each { key ->

        if (!params.containsKey(key)) {
            error "ERROR: Parameter '${key}' is listed as a paramater to validate but is not defined as an actual parameter."
        }

        def value = params[key]

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

        log.info "✓ ${key} OK"
    }

    log.info "All configured parameter paths validated ✓"

	emit:
	global_parameters = global_parameters
}

