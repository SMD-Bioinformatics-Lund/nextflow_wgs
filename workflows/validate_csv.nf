nextflow.enable.dsl=2

workflow VALIDATE_SAMPLES_CSV {

	take:
	samples_csv

	main:

	def csv = file(samples_csv)
	def errorMessages = []

	if (!csv.exists())
		errorMessages << "Samples CSV does not exist: ${csv}"
		error "Samples CSV does not exist: ${csv}"

	def lines = csv.readLines()

	if (lines.size() < 2)
		errorMessages << "CSV is empty or has no data rows: ${csv}"
		error "Samples CSV is empty or has no data rows"

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

	if (missing)
		errorMessages << "Missing columns: ${missing.join(', ')}"
		error "Missing required columns: ${missing.join(', ')}"

	def rows = []
	lines.drop(1).eachWithIndex { line, i ->
		def vals = line.split(',', -1)
		if (vals.size() != header.size())
			def msg = "Line ${idx+2} has ${vals.size()} columns, expected ${header.size()}"
			errorMessages << msg
			error msg

		rows << [header, vals].transpose().collectEntries()
	}

    rows.each { r ->
        if (!(r.sex in ['M','F'])) {
            def msg = "Invalid sex '${r.sex}' for sample ${r.id}"
            errorMessages << msg
            error msg
        }
    }

	def probandsByGroup = rows
	.findAll { it.type == 'proband' }
	.groupBy { it.group }

	probandsByGroup.each { groupId, probands ->
		if (probands.size() > 1)
			errorMessages << "Group '${groupId}' has multiple probands: ${probands*.id.join(', ')}"
			error "Group '${groupId}' has multiple probands: ${probands*.id.join(', ')}"
	}

	def ids = rows*.id
	def dup = ids.findAll { id -> ids.count(id) > 1 }.unique()
	if (dup)
		errorMessages << "Duplicate sample IDs: ${dup.join(', ')}"
		error "Duplicate sample IDs: ${dup.join(', ')}"

	rows.each {
		if (!(it.sex in ['M','F']))
			errorMessages << "Invalid sex '${it.sex}' for sample ${it.id}"
			error "Invalid sex '${it.sex}' for sample ${it.id}"
	}

	def probands = rows.findAll { it.type == 'proband' }
	if (!probands)
		errorMessages << "At least one sample must have type=proband"
		error "At least one sample must have type=proband"

	def byId = rows.collectEntries { [(it.id): it] }
	def groupBySample = rows.collectEntries { [(it.id): it.group] }

	rows.each { r ->
		['read1','read2'].each { fq ->
			if (!file(r[fq]).exists())
				errorMessages << "FASTQ not found: ${r[fq]} (sample ${r.id})"
				error "FASTQ not found: ${r[fq]} (sample ${r.id})"
		}
	}

	def parentUsage = [:].withDefault { [] }
	probands.each { p ->
		def fatherId = p.father
		def motherId = p.mother
		def groupId = p.group

		/* single sample */
		if (!fatherId && !motherId)
			return

		/* partial trio not allowed */
		if (!fatherId || !motherId)
			errorMessages << "Proband ${p.id} must define both parents or none"
			error "Proband ${p.id} must define both parents or none"

		/* parent existence */
		if (!byId.containsKey(fatherId))
			errorMessages << "Father '${fatherId}' not found in CSV (proband ${p.id})"
			error "Father '${fatherId}' not found in CSV (proband ${p.id})"

		if (!byId.containsKey(motherId))
			errorMessages << "Mother '${motherId}' not found in CSV (proband ${p.id})"
			error "Mother '${motherId}' not found in CSV (proband ${p.id})"

		def father = byId[fatherId]
		def mother = byId[motherId]

		/* sex checks */
		if (father.sex != 'M')
			errorMessages << "Father '${fatherId}' must have sex=M (found '${father.sex}')"
			error "Father '${fatherId}' must have sex=M (found '${father.sex}')"

		if (mother.sex != 'F')
			errorMessages << "Mother '${motherId}' must have sex=F (found '${mother.sex}')"
			error "Mother '${motherId}' must have sex=F (found '${mother.sex}')"
		
		/* parents must be in same group */
		if (groupBySample[fatherId] != groupId)
			errorMessages << "Father '${fatherId}' is in group '${groupBySample[fatherId]}', expected '${groupId}' (proband ${p.id})"
			error "Father '${fatherId}' is in group '${groupBySample[fatherId]}', expected '${groupId}' (proband ${p.id})"

		if (groupBySample[motherId] != groupId)
			errorMessages << "Mother '${motherId}' is in group '${groupBySample[motherId]}', expected '${groupId}' (proband ${p.id})"
			error "Mother '${motherId}' is in group '${groupBySample[motherId]}', expected '${groupId}' (proband ${p.id})"

		/* track parent reuse */
		parentUsage[fatherId] << groupId
		parentUsage[motherId] << groupId
	}

	log.info "✔ CSV validation passed (${rows.size()} samples)"

	emit:
	validated_csv = samples_csv
	validation_errors = errorMessages
}

