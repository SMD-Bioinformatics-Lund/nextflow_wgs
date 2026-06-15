nextflow.enable.dsl=2

workflow VALIDATE_PARAMETERS {

	take:
	params_to_validate

	main:

	log.info "Validating file and directory parameters..."

	params_to_validate.each { key ->
		
		if (!params.containsKey(key)) {
			log.info "Parameter '${key}' is listed in params_to_validate but is not defined in params."
		} else {
			validatePathValue(key, params[key])
		}

	}

	log.info "All configured parameter paths validated ✓"

	emit:
	params_to_validate = params_to_validate
}

def validatePathValue(label, value) {
	if (value instanceof Map) {
		value.each { key, nestedValue ->
			validatePathValue("${label}.${key}", nestedValue)
		}
		return
	}

	if (value instanceof List) {
		value.eachWithIndex { nestedValue, idx ->
			validatePathValue("${label}[${idx}]", nestedValue)
		}
		return
	}

	if (!(value instanceof String)) {
		error "ERROR: Parameter '${label}' : '${value}' is not a string path."
	}

	def f = file(value)

	if (!f.exists()) {
		error "ERROR: Param '${label}' points to missing path: ${value}"
	}

	if (f.isFile() && f.size() == 0) {
		error "ERROR: Param '${label}' points to empty file: ${value}"
	}
}