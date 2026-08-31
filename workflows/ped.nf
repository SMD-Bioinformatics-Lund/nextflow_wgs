workflow PED {

	take:
	ch_proband_meta            // channel: [val(group), val(id), val(meta)]
	val_mode                   // string:  Analysis mode.
	val_create_alt_affect_ped  // bool:    Whether family runs should create alternate affected-parent PEDs.
	val_accessdir              // string:  Access base directory used in Scout INFO outputs.

	main:
	ch_versions = channel.empty()
	ch_output_info = channel.empty()

	create_ped(ch_proband_meta)
	ch_ped_base = create_ped.out.ped_base

	ch_ped_trio_affected_permutations = channel.empty()
	ch_ped_trio_affected_permutations = ch_ped_trio_affected_permutations.mix(ch_ped_base)

	if (val_mode == "family" && val_create_alt_affect_ped) {
		ch_ped_trio_affected_permutations = ch_ped_trio_affected_permutations
			.mix(create_ped.out.ped_fa)
			.mix(create_ped.out.ped_ma)

		madeline(
			ch_ped_trio_affected_permutations,
			val_accessdir
		)

		ch_versions = ch_versions.mix(madeline.out.versions.first())
		ch_output_info = ch_output_info.mix(madeline.out.madde_INFO)
	}

	emit:
	ped_base                       = ch_ped_base
	ped_trio_affected_permutations = ch_ped_trio_affected_permutations
	versions                       = ch_versions
	output_info                    = ch_output_info
}


// Create ped
process create_ped {
	tag "${meta.group}"
	time '20m'
	publishDir "${params.outdir}/${params.subdir}/ped", mode: 'copy' , overwrite: true
	memory '1 GB'
    container "${params.container_perl}"
	input:
		tuple val(group), val(id), val(meta)

	output:
		tuple val(meta.group), val(meta.type), path("${meta.group}_base.ped"), emit: ped_base
		tuple val(meta.group), val("ma"), path("${meta.group}_ma.ped"), emit: ped_ma, optional: true
		tuple val(meta.group), val("fa"), path("${meta.group}_fa.ped"), emit: ped_fa, optional: true

	script:
        def father = meta.father == "" ? "0" : meta.father
        def mother = meta.mother == "" ? "0" : meta.mother
        """
	    create_ped.pl --mother ${mother} --father ${father} --group ${meta.group} --id ${meta.id} --sex ${meta.sex}
	    """

	stub:
		"""
		touch "${meta.group}_base.ped"
		touch "${meta.group}_ma.ped"
		touch "${meta.group}_fa.ped"
		"""
}

// Madeline ped, run if family mode
process madeline {
	publishDir "${params.outdir}/${params.subdir}/ped", mode: 'copy' , overwrite: true, pattern: '*.xml'
	memory '1 GB'
	time '1h'
	cpus 2
	container  "${params.container_madeline}"

	input:
		tuple val(group), val(type), path(ped)
		val(accessdir)

	output:
		path("${ped}.madeline.xml")
		tuple val(group), path("${group}_madde.INFO"), emit: madde_INFO
		path "*versions.yml", emit: versions

	script:
		"""
		source activate tools
		ped_parser \\
			-t ped $ped \\
			--to_madeline \\
			-o ${ped}.madeline
		madeline2 \\
			-L "IndividualId" ${ped}.madeline \\
			-o ${ped}.madeline \\
			-x xml
		echo "MADDE	$type ${accessdir}/ped/${ped}.madeline.xml" > ${group}_madde.INFO

		${madeline_version(task)}
		"""

	stub:
		"""
		source activate tools
		touch "${group}_madde.INFO"
		touch "${ped}.madeline.xml"

		${madeline_version(task)}
		"""
}
def madeline_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    ped-parser: \$(echo \$(ped_parser --version 2>&1) | sed -e "s/^.*ped_parser version: //")
	    madeline: \$(echo \$(madeline2 --version 2>&1) | grep : | sed -e"s/^.*Madeline //; s/PDE : 1.*//")
	END_VERSIONS
	"""
}
