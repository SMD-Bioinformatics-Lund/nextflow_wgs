nextflow.enable.dsl=2

workflow VALIDATE_SAMPLES_CSV {

	take:
	samples_csv

	main:

	def csv = file(samples_csv)
	def errorMessages = []

	if (!csv.exists()) {
		def msg = "Samples CSV does not exist: ${csv}"
		errorMessages.add(msg)
	}

	def lines = csv.readLines()

	if (lines.size() < 2) {
		def msg = "Samples CSV is empty or has no data rows"
		errorMessages.add(msg)
	}

	def REQUIRED_COLUMNS = [
		'clarity_sample_id',
		'id',
		'type',
		'assay',
		'sex',
		'diagnosis',
		'phenotype',
		'group',
		'father',
		'mother',
		'platform',
		'read1',
		'read2',
		'analysis'
	]

	def header = lines[0].split(',', -1)*.trim()
	def missing = REQUIRED_COLUMNS - header

	if (missing) {
		def msg = "Missing required columns: ${missing.join(', ')}"
		errorMessages.add(msg)
	}	

	def rows = []
	lines.drop(1).eachWithIndex { line, i ->
		def vals = line.split(',', -1)
		if (vals.size() != header.size()) {
			def msg = "Line ${idx+2} has ${vals.size()} columns, expected ${header.size()}"
			errorMessages.add(msg)
		}
		rows << [header, vals].transpose().collectEntries()
	}

	def probandsByGroup = rows
	.findAll { it.type == 'proband' }
	.groupBy { it.group }

	probandsByGroup.each { groupId, probands ->
		if (probands.size() > 1) {
			def msg = "Group '${groupId}' has multiple probands: ${probands*.id.join(', ')}"
			errorMessages.add(msg)
		}
	}

	def ids = rows*.id
	def dup = ids.findAll { id -> ids.count(id) > 1 }.unique()
	if (dup) {
		def msg = "Duplicate sample IDs: ${dup.join(', ')}"
		errorMessages.add(msg)
	}

	rows.each {
		if (!(it.sex in ['M','F'])) {
			def msg = "Invalid sex '${it.sex}' for sample ${it.id}"
			errorMessages.add(msg)
		}
	}

	def probands = rows.findAll { it.type == 'proband' }
	if (!probands) {
		def msg = "At least one sample must have type=proband"
		errorMessages.add(msg)
	}

	def byId = rows.collectEntries { [(it.id): it] }
	def groupBySample = rows.collectEntries { [(it.id): it.group] }

	rows.each { r ->
		['read1','read2'].each { fq ->
			if (!file(r[fq]).exists()) {
				def msg = "FASTQ not found: ${r[fq]} (sample ${r.id})"
				errorMessages.add(msg)
			}
		}
	}

	def parentUsage = [:].withDefault { [] }
	probands.each { p ->
		def fatherId = p.father
		def motherId = p.mother
		def groupId = p.group

		/* single sample */
		if ( rows.size() == 1 || (!fatherId && !motherId ))
			return

		/* partial trio not allowed */
		if (!fatherId || !motherId) {
			def msg = "Proband ${p.id} must define both parents or none"
			errorMessages.add(msg)
		}
		/* parent existence */
		if (!byId.containsKey(fatherId)) {
			def msg = "Father '${fatherId}' not found in CSV (proband ${p.id})"
			errorMessages.add(msg)
		}

		if (!byId.containsKey(motherId)) {
			def msg = "Mother '${motherId}' not found in CSV (proband ${p.id})"
			errorMessages.add(msg)
		}

		def father = byId[fatherId]
		def mother = byId[motherId]

		/* sex checks */
		if (father.sex != 'M') {
			def msg = "Father '${fatherId}' must have sex=M (found '${father.sex}')"
			errorMessages.add(msg)
		}

		if (mother.sex != 'F') {
			def msg = "Mother '${motherId}' must have sex=F (found '${mother.sex}')"
			errorMessages.add(msg)
		}
		
		/* parents must be in same group */
		if (groupBySample[fatherId] != groupId) {
			def msg = "Father '${fatherId}' is in group '${groupBySample[fatherId]}', expected '${groupId}' (proband ${p.id})"
			errorMessages.add(msg)
		}

		if (groupBySample[motherId] != groupId) {
			def msg = "Mother '${motherId}' is in group '${groupBySample[motherId]}', expected '${groupId}' (proband ${p.id})"
			errorMessages.add(msg)
		}

		/* track parent reuse */
		parentUsage[fatherId] << groupId
		parentUsage[motherId] << groupId
	}

	
	if (errorMessages) {
		error """
			CSV validation failed with ${errorMessages.size()} error(s):

			${errorMessages.join('\n')}
		"""
	}
	else {
		log.info "✔ CSV validation passed (${rows.size()} samples)"
	}
	emit:
	validated_csv = samples_csv
	validation_errors = Channel.value(errorMessages)
}

