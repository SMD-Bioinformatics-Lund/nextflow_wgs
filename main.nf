#!/usr/bin/env nextflow

include { SNV_ANNOTATE } from './workflows/annotate_snvs.nf'
include { IDSNP_CALL } from './modules/idsnp.nf'
include { IDSNP_VCF_TO_JSON } from './modules/idsnp.nf'

nextflow.enable.dsl=2


workflow {

	// Print startup and conf output dirs and modes.

	// TODO: Params assignment inside workflow block is a temp solution:
	//       outdir and subdir need to be combined in new var, re-setting
	//       params.outdir won't work.
	params.results_output_dir = params.outdir + '/' + params.subdir

	// TODO: Pass these to processes in meta?
	params.mode = file(params.csv).countLines() > 2 ? "family" : "single"
	params.trio = file(params.csv).countLines() > 3 ? true : false

	// Check whether genome assembly is indexed //
	// TODO: Move to some pre-processing workflow:
	if(params.genome_file) {
		_bwaId = Channel
			.fromPath("${params.genome_file}.bwt")
			.ifEmpty { exit 1, "BWA index not found: ${params.genome_file}.bwt" }
	}

	// Count lines of input csv, if more than 2(header + 1 ind) then mode is set to family //
	log.info("Input CSV: " + params.csv)
	log.info("mode: " + params.mode)
	log.info("trio analysis: " + params.trio)
	log.info("Results output dir: " + params.results_output_dir)
	log.info("Results subdir: " + params.subdir)
	log.info("CRON output dir: " + params.crondir)

	// Print commit-version of active deployment
	// TODO: stuff this one into versions too, for good measure.
	file(params.git)
		.readLines()
		.each { println "git commit-hash: " + it }

	// Print active container
	log.info("container: " + file(params.container).toRealPath())

	ch_versions = Channel.empty()

	Channel
		.fromPath(params.csv)
		.splitCsv(header: true)
		.set { ch_samplesheet }

	NEXTFLOW_WGS(ch_samplesheet)

	ch_versions = ch_versions.mix(NEXTFLOW_WGS.out.versions).collect{ it }

	combine_versions(
		ch_samplesheet
			.first()
			.map{ row -> row.group },
		ch_versions
	)

}

// TODO: needs to be moved into process/workflow to silence lsp error:
workflow.onComplete {

		def completed_at = "${workflow.complete}"

		def msg = """\
		Pipeline execution summary
		---------------------------
		Completed at: ${completed_at}
		Duration    : ${workflow.duration}
		Success     : ${workflow.success}
		scriptFile  : ${workflow.scriptFile}
		workDir     : ${workflow.workDir}
		csv         : ${params.csv}
		exit status : ${workflow.exitStatus}
		errorMessage: ${workflow.errorMessage}
		errorReport :
			"""
			.stripIndent()
		def error = """\
			${workflow.errorReport}
		"""
			.stripIndent()

		def base = file(params.csv).getBaseName()
		File logFile = new File("${params.crondir}/logs/${base}.complete")
		if (!logFile.getParentFile().exists()) {
			logFile.getParentFile().mkdirs()
		}
		logFile.text = msg

	if (error) {
		def error_report = """\
		errorReport :
			${error}
			"""
		logFile.append(error_report)
	}
}


workflow NEXTFLOW_WGS {

	take:
	ch_samplesheet

	main:
	// Output channels:
	ch_versions    = Channel.empty() // Gather software versions
	ch_output_info = Channel.empty() // Gather data for .INFO
	ch_qc_json     = Channel.empty() // Gather and merge QC JSONs per sample

	// CHANNEL PREP //
	// TODO: Better solution for this. Assume shomehow that everything non-bam/non-vcf is fq.
	ch_fastq = ch_samplesheet
		.filter {
			row -> row.read1.endsWith("q.gz") && row.read2.endsWith("q.gz")
		}
		.map { row ->
			def group = row.group
			def id = row.id
			def fastq_r1 = row.read1
			def fastq_r2 = row.read2
			tuple(group, id, fastq_r1, fastq_r2) // TODO: filter non fq
	}


	ch_vcf_annotation_only = ch_samplesheet
		.filter {
			row -> row.read1.endsWith(".vcf") || row.read1.endsWith(".vcf.gz")
		}
		.map { row ->
			def group = row.group
			def id = row.id
			def vcf = row.read1
			tuple(group, id, vcf)
		}

	// GATK Ref:
	ch_gatk_ref = Channel
		.fromPath(params.gatkreffolders)
		.splitCsv(header:true)
		.map{ row-> tuple(row.i, row.refpart) }

	// TODO: all the meta channels below should be consolidated into some meta object?

	// meta sent to split_normalize_mito, eklipse, bgzip_scored_genmod and used
	// to trigger loqusdb sv dummy
	ch_meta = ch_samplesheet.map{ row ->
		tuple(
			row.group,
			row.id,
			row.sex,
			row.type)
	}

	ch_gatkcov_meta = ch_meta

	ch_split_normalize_meta = ch_samplesheet
		.filter { row ->
			row.type == "proband"
		}
		.map{ row ->
			tuple(
				row.group,
				row.id,
				row.sex,
				row.type
			)
		}

	ch_expansionhunter_meta = ch_split_normalize_meta
	ch_svvcf_to_bed_meta  = ch_split_normalize_meta

	// meta sent to qc_to_cdm
	ch_qc_extra = ch_samplesheet
		.map { row ->
			def group = row.group
			def id = row.id
			def diagnosis = row.diagnosis
			def fastq_r1 = row.read1
			def fastq_r2 = row.read2
			tuple(group, id, diagnosis, fastq_r1, fastq_r2)
		}

	// meta sent to stranger
	ch_stranger_meta = ch_samplesheet
		.filter { row ->
			row.type == "proband"
		}
		.map { row -> // TODO: Rename this channel
			tuple(
				row.group,
				row.id,
				row.sex,
				row.mother,
				row.father,
				row.phenotype,
				row.diagnosis,
				row.type,
				row.assay,
				row.clarity_sample_id,
				(row.containsKey("ffpe") ? row.ffpe : false),
				(row.containsKey("analysis") ? row.analysis : false)
			)
		}

	// meta sent to create_yml
	ch_scout_yaml_meta = ch_samplesheet
		.filter { row ->
			row.type == "proband"
		}
		.map { row ->
			tuple(
				row.group,
				row.id,
				row.diagnosis,
				row.assay,
				row.type,
				row.clarity_sample_id,
				(row.containsKey("analysis") ? row.analysis : false)
			)
		}

	// BAM-start
	// Check for .bam files in read1 and start from bam if any found.
	ch_bam_start = ch_samplesheet
		.filter {
			row -> row.read1.endsWith("bam")
		}
		.map {
			row ->
			def group = row.group
			def id = row.id
			def bam = row.read1
			def bai = row.read2
			tuple(group, id, bam, bai)
		}



	rename_bam(ch_bam_start) // See process for more info about why this is needed.
	copy_bam(rename_bam.out.bam_bai)
	bamtoyaml(ch_bam_start)
	ch_output_info = ch_output_info.mix(bamtoyaml.out.bamchoice_INFO)

	ch_bam_start_dedup_dummy = Channel.empty()
	if(params.run_melt) {
		dedupdummy(ch_bam_start)
		ch_bam_start_dedup_dummy = dedupdummy.out.dedup_dummy
	}

	ch_bam_bai = Channel.empty()
	ch_bam_bai = ch_bam_bai.mix(copy_bam.out.bam_bai)

	// PED //
	ch_ped_input = ch_samplesheet
		.filter { row -> row.type == "proband" }
		.map { row ->
			def group = row.group
			def id = row.id
			def type = row.type
			def sex = row.sex
			def father = row.father
			def mother = row.mother
			tuple(group, id, type, sex, mother, father)
		}

	create_ped(ch_ped_input)
	ch_ped_base = create_ped.out.ped_base
	ch_ped_father_affected = Channel.empty()
	ch_ped_mother_affected = Channel.empty()

	ch_ped_trio_affected_permutations = Channel.empty()  // Channel for base ped + father and mother set as affected peds
	ch_ped_trio_affected_permutations = ch_ped_trio_affected_permutations.mix(ch_ped_base)
	if(params.mode == "family" && params.assay == "wgs") {

		ch_ped_father_affected = create_ped.out.ped_fa
		ch_ped_mother_affected = create_ped.out.ped_ma

		ch_ped_trio_affected_permutations = ch_ped_trio_affected_permutations.mix(ch_ped_father_affected).mix(ch_ped_mother_affected)
		madeline(ch_ped_trio_affected_permutations)

		ch_versions = ch_versions.mix(madeline.out.versions.first())
		ch_output_info = ch_output_info.mix(madeline.out.madde_INFO)

	}

	// FASTQ //
	if (params.umi) {
		fastp(ch_fastq)
		ch_fastq = fastp.out.fastq_trimmed_reads
		ch_versions = ch_versions.mix(fastp.out.versions.first())
	}

	// ALIGN //
	//TODO: why do we have a params.align conditional anyway?
	ch_dedup_stats = Channel.empty()
	if (params.align) {
		bwa_align(ch_fastq)
		markdup(bwa_align.out.bam_bai)
		ch_dedup_stats = ch_dedup_stats.mix(markdup.out.dedup_metrics)
		ch_output_info = ch_output_info.mix(markdup.out.dedup_bam_INFO)
		ch_bam_bai = ch_bam_bai.mix(markdup.out.dedup_bam_bai)

		ch_versions = ch_versions.mix(bwa_align.out.versions.first())
		ch_versions = ch_versions.mix(markdup.out.versions.first())
	}

	bqsr(ch_bam_bai)
	ch_versions = ch_versions.mix(bqsr.out.versions.first())

	// POST SEQ QC //
	sentieon_qc(ch_bam_bai)
	ch_versions = ch_versions.mix(sentieon_qc.out.versions.first())

	ch_dedup_stats = ch_dedup_stats.mix(ch_bam_start_dedup_dummy)

	sentieon_qc_postprocess(
		sentieon_qc.out.sentieon_qc_metrics.join(ch_dedup_stats, by: [0,1])
	)

	ch_qc_json = ch_qc_json.mix(sentieon_qc_postprocess.out.qc_json)

	IDSNP_CALL(ch_bam_bai, params.idsnps)
	IDSNP_VCF_TO_JSON(IDSNP_CALL.out.vcf)
	ch_versions = ch_versions.mix(IDSNP_CALL.out.versions.first())

	// COVERAGE //
	d4_coverage(ch_bam_bai)
	ch_versions = ch_versions.mix(d4_coverage.out.versions.first())
	ch_output_info = ch_output_info.mix(d4_coverage.out.d4_INFO)

	if (params.gatkcov) {
		gatkcov(ch_gatkcov_meta.join(ch_bam_bai, by: [0, 1]))
		ch_versions = ch_versions.mix(gatkcov.out.versions.first())
	}

	if (params.assay == "swea") {
		depth_onco(ch_bam_bai)
	}

	// CONTAMINATION //
	if (params.antype == "wgs") {
		verifybamid2(ch_bam_bai)
		ch_versions = ch_versions.mix(verifybamid2.out.versions.first())
	}

	// SNV CALLING //
	dnascope(ch_bam_bai.join(bqsr.out.dnascope_bqsr, by: [0, 1]))
	ch_versions = ch_versions.mix(dnascope.out.versions.first())
	gvcf_combine(dnascope.out.gvcf_tbi.groupTuple())
	ch_versions = ch_versions.mix(gvcf_combine.out.versions.first())

	ch_split_normalize = gvcf_combine.out.combined_vcf
	ch_split_normalize_concat_vcf = Channel.empty()

	// TODO: move antypes and similar to constants?
	if (params.antype == "panel") {
		freebayes(ch_bam_bai)
		ch_versions = ch_versions.mix(freebayes.out.versions.first())
		ch_split_normalize_concat_vcf = ch_split_normalize_concat_vcf.mix(freebayes.out.freebayes_variants)
	}

	// MITO SIDE-QUEST
	if (params.antype == "wgs") { // TODO: if params.mito etc ? will probably mess up split_normalize

		fetch_MTseqs(ch_bam_bai)

		ch_output_info = ch_output_info.mix(fetch_MTseqs.out.mtBAM_INFO)

		// MITO BAM QC
		sentieon_mitochondrial_qc(fetch_MTseqs.out.bam_bai)

		build_mitochondrial_qc_json(sentieon_mitochondrial_qc.out.qc_tsv)

		ch_qc_json = ch_qc_json.mix(build_mitochondrial_qc_json.out.qc_json)

		// SNVs
		ch_mutect2_input = fetch_MTseqs.out.bam_bai.groupTuple()
		run_mutect2(ch_mutect2_input)

		split_normalize_mito(
			run_mutect2.out.vcf,
			ch_split_normalize_meta
			)

		run_hmtnote(split_normalize_mito.out.vcf)
		ch_split_normalize_concat_vcf = ch_split_normalize_concat_vcf.mix(run_hmtnote.out.vcf)

		run_haplogrep(run_mutect2.out.vcf)
		ch_output_info = ch_output_info.mix(run_haplogrep.out.haplogrep_INFO)

		// SVs
		run_eklipse(ch_meta.join(fetch_MTseqs.out.bam_bai, by : [0, 1]))
		ch_output_info = ch_output_info.mix(run_eklipse.out.eklipse_INFO)

		// MITO VERSIONS
		ch_versions = ch_versions.mix(run_hmtnote.out.versions.first())
		ch_versions = ch_versions.mix(split_normalize_mito.out.versions.first())
		ch_versions = ch_versions.mix(run_mutect2.out.versions.first())
		ch_versions = ch_versions.mix(sentieon_mitochondrial_qc.out.versions.first())
		ch_versions = ch_versions.mix(fetch_MTseqs.out.versions.first())
		ch_versions = ch_versions.mix(run_eklipse.out.versions.first())
		ch_versions = ch_versions.mix(run_haplogrep.out.versions.first())
	}

	// SNV ANNOTATION
	if (params.annotate) {

		split_normalize(ch_split_normalize, ch_split_normalize_concat_vcf)

		ch_snv_indel_vcf = split_normalize.out.intersected_vcf.mix(ch_vcf_annotation_only)

		SNV_ANNOTATE(ch_snv_indel_vcf, ch_ped_trio_affected_permutations)
		ch_versions = ch_versions.mix(SNV_ANNOTATE.out.versions)
		ch_output_info = ch_output_info.mix(SNV_ANNOTATE.out.output_info)

		// SNPs
		ch_peddy_input_vcf = SNV_ANNOTATE.out.annotated_snv_vcf
			.filter { it ->
				def type = it[1]
				type == "proband"
			}

		// TODO: Move this guy to QC:
		peddy(ch_peddy_input_vcf.join(ch_ped_base, by: [0,1]))
		ch_output_info = ch_output_info.mix(peddy.out.peddy_INFO)


		ch_versions = ch_versions.mix(peddy.out.versions.first())

		if (params.antype == "wgs") {
			// fastgnomad
			fastgnomad(split_normalize.out.norm_uniq_dpaf_vcf)


			// upd
			ch_upd_meta = ch_samplesheet
				.filter { row ->
					row.type == "proband"
				}
				.map { row ->
					tuple(row.group, row.id, row.mother, row.father)
				}

			// upd
			upd(fastgnomad.out.vcf, ch_upd_meta)
			upd_table(upd.out.upd_sites)

			// roh
			roh(fastgnomad.out.vcf)
			overview_plot(
				gatkcov.out.cov_plot.filter{ it ->
					it[2] == "proband"
				},
				roh.out.roh_plot,
				upd.out.upd_bed
			)

			generate_gens_data(dnascope.out.gvcf_tbi.join(gatkcov.out.cov_gens, by: [0,1]))

			ch_gens_v4_meta = gatkcov.out.cov_plot
				.combine(roh.out.roh_plot.map { it -> it[1] })
				.combine(upd.out.upd_bed)
				.combine(upd.out.upd_sites.map { it -> it[1] })

			generate_gens_v4_meta(ch_gens_v4_meta)

			ch_cron_meta = generate_gens_v4_meta.out.meta
				.join(generate_gens_data.out.is_done, by: [0,1])

			gens_v4_cron(ch_cron_meta)

			ch_output_info = ch_output_info.mix(overview_plot.out.oplot_INFO)

			ch_versions = ch_versions.mix(upd.out.versions.first())
			ch_versions = ch_versions.mix(roh.out.versions.first())
		}
	}

	ch_loqusdb_sv = Channel.empty()

	if (params.sv) {
		ch_smn_tsv = Channel.empty()
		ch_panel_svs_present = Channel.empty()
		ch_panel_svs_absent = Channel.empty()
		ch_postprocessed_merged_sv_vcf = Channel.empty()
		if(params.antype  == "wgs") {
			// SMN CALLING //
			SMNCopyNumberCaller(ch_bam_bai)
			ch_output_info = ch_output_info.mix(SMNCopyNumberCaller.out.smn_INFO)

			// Collects each individual's SMNCNC-tsv and creates one tsv-file
			ch_smn_tsv = ch_smn_tsv.mix(
				SMNCopyNumberCaller.out
					.smn_tsv.collectFile(
						keepHeader: true,
						storeDir: "${params.results_output_dir}/smn/")
			)

			// CALL REPEATS //

			expansionhunter(ch_bam_bai.join(ch_expansionhunter_meta, by: [0,1]))
			stranger(expansionhunter.out.expansionhunter_vcf)
			vcfbreakmulti_expansionhunter(
				stranger.out.vcf_annotated,
				ch_stranger_meta
			)
			ch_output_info = ch_output_info.mix(vcfbreakmulti_expansionhunter.out.str_INFO)
			reviewer(expansionhunter.out.bam_vcf)

			ch_versions = ch_versions.mix(SMNCopyNumberCaller.out.versions.first())
			ch_versions = ch_versions.mix(reviewer.out.versions.first())
			ch_versions = ch_versions.mix(expansionhunter.out.versions.first())
			ch_versions = ch_versions.mix(vcfbreakmulti_expansionhunter.out.versions.first())
			ch_versions = ch_versions.mix(stranger.out.versions.first())
		}


		// BIG SV //

		gatk_coverage(ch_bam_bai)
		ch_gatk_coverage = gatk_coverage.out.coverage_tsv
		gatk_call_ploidy(ch_gatk_coverage)
		ch_gatk_ploidy = gatk_call_ploidy.out.call_ploidy

		// TODO: do the joining and combining outside
		gatk_call_cnv(ch_gatk_coverage.join(ch_gatk_ploidy, by: [0,1]).combine(ch_gatk_ref))

		ch_gatk_postprocess_input = gatk_call_cnv.out.gatk_calls
			.groupTuple(by : [0,1])
			.join(ch_gatk_ploidy, by: [0, 1])
			.combine(ch_gatk_ref.groupTuple(by : [3])) // TODO: Explanation here.

		postprocessgatk(ch_gatk_postprocess_input)
		filter_merge_gatk(postprocessgatk.out.called_gatk)
		ch_filtered_merged_gatk_calls = filter_merge_gatk.out.merged_filtered_vcf


		ch_versions = ch_versions.mix(gatk_call_cnv.out.versions.first())
		ch_versions = ch_versions.mix(gatk_coverage.out.versions.first())
		ch_versions = ch_versions.mix(gatk_call_ploidy.out.versions.first())
		ch_versions = ch_versions.mix(postprocessgatk.out.versions.first())


		ch_manta_out = Channel.empty()
		if (params.antype == "wgs") {
			manta(ch_bam_bai)
			ch_manta_out = ch_manta_out.mix(manta.out.vcf)
			tiddit(ch_bam_bai)
			svdb_merge(
				ch_manta_out.groupTuple(),
				tiddit.out.vcf.groupTuple(),
				ch_filtered_merged_gatk_calls.groupTuple()
			)
			ch_postprocessed_merged_sv_vcf = ch_postprocessed_merged_sv_vcf.mix(svdb_merge.out.merged_bndless_vcf)
			ch_panel_svs_present = ch_postprocessed_merged_sv_vcf
			ch_loqusdb_sv = ch_loqusdb_sv.mix(svdb_merge.out.merged_vcf)

			ch_versions = ch_versions.mix(manta.out.versions.first())
			ch_versions = ch_versions.mix(tiddit.out.versions.first())
			ch_versions = ch_versions.mix(svdb_merge.out.versions.first())
		}

		// MELT //
		// TODO: The panel SV-calling code presumes melt is called so just move the process code there:
		if (params.run_melt) {

		    ch_melt_qc_vals = sentieon_qc_postprocess.out.qc_json.map { item ->
				def group = item[0]
				def id = item[1]
				def qc_json = item[2]

				println "qc_json from tuple: ${qc_json}"

				def ins_dev = "NA"   // Default value for stub runs
				def coverage = "NA"
				def ins_size = "NA"

				def INS_SIZE
				def MEAN_DEPTH
				def COV_DEV

				if (qc_json.exists() && qc_json.size() > 0) {
					println "Reading qc_json file: ${qc_json}"

					qc_json.readLines().each { line ->
						if (line =~ /\"(ins_size_dev)\" : \"(\S+)\"/) {
							ins_dev = (line =~ /\"(ins_size_dev)\" : \"(\S+)\"/)[0][2]  // Extract matched value
						}
						if (line =~ /\"(mean_coverage)\" : \"(\S+)\"/) {
							coverage = (line =~ /\"(mean_coverage)\" : \"(\S+)\"/)[0][2]
						}
						if (line =~ /\"(ins_size)\" : \"(\S+)\"/) {
							ins_size = (line =~ /\"(ins_size)\" : \"(\S+)\"/)[0][2]
						}
					}

					INS_SIZE = ins_size ?: "NA"  // Default to "NA" if not present
					MEAN_DEPTH = coverage ?: "NA"
					COV_DEV = ins_dev ?: "NA"

				} else {
					println "Warning: Empty or missing qc_json file for ${id}. Using default stub values."
					INS_SIZE = "NA"
					MEAN_DEPTH = "NA"
					COV_DEV = "NA"
				}
				tuple(group, id, INS_SIZE, MEAN_DEPTH, COV_DEV)
			}

			melt(ch_bam_bai, ch_melt_qc_vals)
			intersect_melt(melt.out.melt_vcf_nonfiltered)

			ch_versions = ch_versions.mix(melt.out.versions.first())
			ch_versions = ch_versions.mix(intersect_melt.out.versions.first())
		}

		if (params.antype == "panel") {
			ch_panel_merge = Channel.empty()
			manta_panel(ch_bam_bai)
			ch_manta_out = ch_manta_out.mix(manta_panel.out.vcf)


			cnvkit_panel(ch_bam_bai, split_normalize.out.intersected_vcf, ch_melt_qc_vals)
			ch_cnvkit_out = cnvkit_panel.out.cnvkit_calls
			ch_output_info = ch_output_info.mix(cnvkit_panel.out.cnvkit_INFO)

			ch_panel_merge = ch_panel_merge.mix(
				ch_cnvkit_out,
				ch_manta_out,
				ch_filtered_merged_gatk_calls
			).groupTuple()

			svdb_merge_panel(ch_panel_merge)

			ch_loqusdb_sv = ch_loqusdb_sv.mix(svdb_merge_panel.out.loqusdb_vcf)

			postprocess_merged_panel_sv_vcf(svdb_merge_panel.out.merged_vcf, intersect_melt.out.merged_intersected_vcf)
			ch_postprocessed_merged_sv_vcf = ch_postprocessed_merged_sv_vcf.mix(
				postprocess_merged_panel_sv_vcf.out.merged_postprocessed_vcf
			)

			// COUNT NUMBER OF SVs, if no called don't do annotation ///
			ch_panel_svs_check = postprocess_merged_panel_sv_vcf.out.merged_postprocessed_vcf.map { item ->
				def group = item[0]
				def id = item[1]
				def merged_vcf = item[2]

				def has_sv = true // to allow for stub

				if (merged_vcf.exists() && merged_vcf.size() > 0) {
					has_sv = false
					merged_vcf.withReader { reader ->
						String line
						while ((line = reader.readLine()) != null) {
							if (!line.startsWith('#')) {
								has_sv = true
								break
							}
						}
					}
				}

				tuple(group, id, merged_vcf, has_sv)
			}

			// Split into two channels: one with SVs, one without
			ch_panel_svs_present = ch_panel_svs_check
				.filter { group, id, merged_vcf, has_sv -> has_sv }
				.map { group, id, merged_vcf, has_sv -> tuple(group, id, merged_vcf) }

			ch_panel_svs_absent = ch_panel_svs_check
				.filter { group, id, merged_vcf, has_sv -> !has_sv }
				.map { group, id, merged_vcf, has_sv -> tuple(group, "proband", merged_vcf) } //this assumes no trio on panel ever. 

			ch_versions = ch_versions.mix(manta_panel.out.versions.first())
			ch_versions = ch_versions.mix(cnvkit_panel.out.versions.first())
			ch_versions = ch_versions.mix(svdb_merge_panel.out.versions.first())
			ch_versions = ch_versions.mix(postprocess_merged_panel_sv_vcf.out.versions.first())
		}

		// ANNOTATE SVs //
		ch_proband_meta = ch_split_normalize_meta
		filter_proband_null_calls(ch_panel_svs_present,ch_proband_meta)
		tdup_to_dup(filter_proband_null_calls.out.filtered_vcf)
		annotsv(tdup_to_dup.out.renamed_vcf)
		vep_sv(tdup_to_dup.out.renamed_vcf)
		postprocess_vep_sv(vep_sv.out.vep_sv_vcf)
		artefact(postprocess_vep_sv.out.merged_processed_vcf)
		bcftools_annotate_dbvar(artefact.out.vcf)

		ch_versions = ch_versions.mix(annotsv.out.versions.first())
		ch_versions = ch_versions.mix(vep_sv.out.versions.first())
		ch_versions = ch_versions.mix(postprocess_vep_sv.out.versions.first())
		ch_versions = ch_versions.mix(artefact.out.versions.first())
		
		ch_ped_prescore = ch_ped_trio_affected_permutations
		ch_add_annotsv_input = bcftools_annotate_dbvar.out.vcf.join(annotsv.out.annotsv_tsv) // ch: group,  path(vcf), path(tbi), path(annotsv_tsv)
		add_annotsv_to_svvcf(ch_add_annotsv_input)
		add_omim_morbid_to_svvcf(add_annotsv_to_svvcf.out.vcf)
		add_callerpenalties_to_svvcf(add_omim_morbid_to_svvcf.out.vcf)

		ch_add_geneticmodels_to_svvcf_input = add_callerpenalties_to_svvcf.out.vcf.cross(ch_ped_prescore)
			.map{
				item ->
				def group = item[0][0]
				def penalty_vcf = item[0][1]
				def type = item[1][1]
				def ped = item[1][2]
				tuple(group, type, ped, penalty_vcf)
			}
		add_geneticmodels_to_svvcf(ch_add_geneticmodels_to_svvcf_input)
		score_sv(add_geneticmodels_to_svvcf.out.annotated_sv_vcf)
		bgzip_scored_genmod(score_sv.out.scored_vcf.mix(ch_panel_svs_absent))
		ch_output_info = ch_output_info.mix(bgzip_scored_genmod.out.sv_INFO)

		svvcf_to_bed(bgzip_scored_genmod.out.sv_rescore_vcf, ch_svvcf_to_bed_meta)

		ch_compound_finder_input = bgzip_scored_genmod.out.sv_rescore  // Take final scored SV VCF
			.join(ch_ped_trio_affected_permutations, by: [0,1])        // join with correct ped on group, type
			.join(SNV_ANNOTATE.out.annotated_snv_vcf, by: [0,1])               // join with final SNV VCF + index on group, type

		compound_finder(ch_compound_finder_input)
		ch_output_info = ch_output_info.mix(compound_finder.out.svcompound_INFO)

		ch_versions = ch_versions.mix(score_sv.out.versions.first())
		ch_versions = ch_versions.mix(bgzip_scored_genmod.out.versions.first())
		ch_versions = ch_versions.mix(compound_finder.out.versions.first())

		// TODO: streamline if-conditions:
		if(params.antype == "wgs" && params.trio && params.mode == "family") {
			plot_pod(
				fastgnomad.out.vcf,
				bgzip_scored_genmod.out.sv_rescore_vcf.join(ch_ped_base, by: 0),
				ch_meta.filter { row ->
					def type = row[3]
					type == "proband"
				}
			)
		}

	} else {
		// TODO: move the entire else-block to top w/ if-not at the beginning

		// add_to_loqusb won't run if no svvcf is generated
		// the code below creates dummy svvcf for no-SV runs
		def dummy_file = file("NA")
		ch_loqusdb_no_sv_dummy = ch_samplesheet
			.map { row ->
				def group = row.group
				def sv_vcf_dummy = dummy_file
				tuple(group, sv_vcf_dummy)
			}
			.first()

		ch_loqusdb_sv = ch_loqusdb_sv.mix(ch_loqusdb_no_sv_dummy)

	}

	// LOQUSDB //
	add_to_loqusdb(
		ch_ped_base.join(SNV_ANNOTATE.out.annotated_snv_vcf, by: [0,1]),
		ch_loqusdb_sv
	)

	// MERGE QC JSONs AND OUTPUT TO CDM //
	merge_qc_json(ch_qc_json.groupTuple(by: [0,1]))
	ch_qc_to_cdm = merge_qc_json.out.qc_cdm_merged.join(ch_qc_extra, by: [0,1])
	qc_to_cdm(ch_qc_to_cdm)

	// OUTPUT INFO
	output_files(ch_output_info.groupTuple())
	// SCOUT YAML
	create_yaml(ch_scout_yaml_meta, ch_ped_base, output_files.out.yaml_INFO)

	emit:
		versions = ch_versions
}

	// TODO: re-implement annotation-only runs and sort out remainder of this block:
	// Input channels for alignment, variant calling and annotation //

	// annotate_only = Channel.create()

	// If input-files has bam files bypass alignment, otherwise go for fastq-channels => three options for fastq, sharded bwa, normal bwa or umi trimming
	// input_files.view().choice(bam_choice, fastq, fastq_sharded, fastq_umi, annotate_only ) { it[2] =~ /\.bam/ ? 0 : ( it[2] =~ /\.vcf.gz/ ? 4 : (params.shardbwa ? 2 : (params.umi ? 3 : 1) )) }

	// annotate_only.into{
	// 	annotate_only_vep;
	// 	annotate_only_cadd
	// }




process fastp {
	cpus 10
	tag "$id"
	time '1h'
	memory '20 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	container "${params.container_fastp}"

	input:
		tuple val(group), val(id), path(fq_r1), path(fq_r2)

	output:
		tuple val(group), val(id), path("${id}_R1_a_q_u_trimmed.fq.gz"), path("${id}_R2_a_q_u_trimmed.fq.gz"), emit: fastq_trimmed_reads
		path("*versions.yml"), emit: versions

	when:
		params.umi

	script:
		"""
		fastp -i $fq_r1 -I $fq_r2 --stdout \\
			-U --umi_loc=per_read --umi_len=3 \\
			-w ${task.cpus} \\
		| fastp --stdin --interleaved_in -f 2 -F 2 \\
			-o ${id}_R1_a_q_u_trimmed.fq.gz \\
			-O ${id}_R2_a_q_u_trimmed.fq.gz \\
			-l 30 \\
			-w ${task.cpus}
		
		${fastp_version(task)}
		"""

	stub:
		"""
		touch "${id}_R1_a_q_u_trimmed.fq.gz"
		touch "${id}_R2_a_q_u_trimmed.fq.gz"

		${fastp_version(task)}
		"""
}
def fastp_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    fastp: \$(echo \$(fastp -v 2>&1) | cut -f 2 -d " ")
	END_VERSIONS
	"""
}


process bwa_align {
	cpus 50
	memory '100 GB' 	// 64 GB peak giab //
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	tag "$id"
	container  "${params.container_sentieon}"

	input:
		tuple val(group), val(id), path(fastq_r1), path(fastq_r2)

	output:
		tuple val(group), val(id), path("${id}_merged.bam"), path("${id}_merged.bam.bai"), emit: bam_bai
		path "*versions.yml", emit: versions

	when:
		params.align

	script:
		"""
		sentieon bwa mem \\
			-M \\
			-K ${params.bwa_K_size} \\
			-R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' \\
			-t ${task.cpus} \\
			${params.genome_file} $fastq_r1 $fastq_r2 \\
			| sentieon util sort \\
			-r ${params.genome_file} \\
			-o ${id}_merged.bam \\
			-t ${task.cpus} --sam2bam -i -

		${bwa_align_versions(task)}
		"""

	stub:
		"""
		touch "${id}_merged.bam"
		touch "${id}_merged.bam.bai"

		${bwa_align_versions(task)}
		"""

}
def bwa_align_versions(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon util --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	    bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
	END_VERSIONS
	"""
}

process markdup {
	cpus 40
	errorStrategy 'retry'
	maxErrors 5
	tag "$id"
	memory '50 GB' // 12GB peak GIAB
	time '3h'
	container  "${params.container_sentieon}"
	publishDir "${params.results_output_dir}/bam", mode: 'copy' , overwrite: 'true', pattern: '*_dedup.bam*'

	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${id}_dedup.bam"), path("${id}_dedup.bam.bai"), emit: dedup_bam_bai
		tuple val(group), val(id), path("dedup_metrics.txt"), emit: dedup_metrics
		tuple val(group), path("${group}_bam.INFO"), emit: dedup_bam_INFO
		path "*versions.yml", emit: versions

	script:
		"""
		sentieon driver \\
			--temp_dir /local/scratch/ \\
			-t ${task.cpus} \\
			-i $bam --shard 1:1-248956422 --shard 2:1-242193529 --shard 3:1-198295559 --shard 4:1-190214555 --shard 5:1-120339935 --shard 5:120339936-181538259 --shard 6:1-170805979 --shard 7:1-159345973 --shard 8:1-145138636 --shard 9:1-138394717 --shard 10:1-133797422 --shard 11:1-135086622 --shard 12:1-56232327 --shard 12:56232328-133275309 --shard 13:1-114364328 --shard 14:1-107043718 --shard 15:1-101991189 --shard 16:1-90338345 --shard 17:1-83257441 --shard 18:1-80373285 --shard 19:1-58617616 --shard 20:1-64444167 --shard 21:1-46709983 --shard 22:1-50818468 --shard X:1-124998478 --shard X:124998479-156040895 --shard Y:1-57227415 --shard M:1-16569 \\
			--algo LocusCollector \\
			--fun score_info ${id}.score

		sentieon driver \\
			--temp_dir /local/scratch/ \\
			-t ${task.cpus} \\
			-i $bam \\
			--algo Dedup --score_info ${id}.score \\
			--metrics dedup_metrics.txt \\
			--rmdup ${id}_dedup.bam

		# TODO: Build this and other INFO outputs in separate workflow/process.
		echo "BAM	$id	/access/${params.subdir}/bam/${id}_dedup.bam" > ${group}_bam.INFO

		${markdup_versions(task)}
		"""

	stub:
		"""
		touch "${id}_dedup.bam"
		touch "${id}_dedup.bam.bai"
		touch "dedup_metrics.txt"
		touch "${group}_bam.INFO"

 		${markdup_versions(task)}
 		"""
}
def markdup_versions(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}

// Rename of bams prior to copy_bam required so that copy_bam can output the
// bam with unaltered filenames. This needs to be fixed in a better way.
// Please see: https://github.com/SMD-Bioinformatics-Lund/nextflow_wgs/issues/306
process rename_bam {
	tag "$id"
	cpus 1
	memory '2GB'
	time '1h'

	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${bam}.tmp"), path("${bai}.tmp"), emit: bam_bai

	script:
		"""
		ln -s ${bam} "${bam}.tmp"
		ln -s ${bai} "${bai}.tmp"
		"""

	stub:
		"""
		ln -s ${bam} "${bam}.tmp"
		ln -s ${bai} "${bai}.tmp"
		"""
}

process copy_bam {

	tag "$id"
	cpus 1
	memory '2GB'
	time '1h'

	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${id}_dedup.bam"), path("${id}_dedup.bam.bai"), emit: bam_bai
	script:
		"""
		ionice -c 2 -n 7 cp ${bam} "${id}_dedup.bam"
		ionice -c 2 -n 7 cp ${bai} "${id}_dedup.bam.bai"
		"""

	stub:
		"""
		touch "${id}_dedup.bam"
		touch "${id}_dedup.bam.bai"
		"""
}


// For melt to work if started
process dedupdummy {
	input:
		tuple val(group), val(id), path(bam), path(bai)
	output:
		tuple val(group), val(id), path("dummy"), emit: dedup_dummy

	when:
		params.run_melt

	script:
	"""
	echo test > dummy
	"""
}


process bqsr {
	cpus 40
	errorStrategy 'retry'
	maxErrors 5
	tag "$id"
	memory '30 GB'
	// 12gb peak giab //
	time '5h'
	container  "${params.container_sentieon}"
	publishDir "${params.results_output_dir}/bqsr", mode: 'copy' , overwrite: 'true', pattern: '*.table'

	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${id}.bqsr.table"), emit: dnascope_bqsr
		path "*versions.yml", emit: versions

	script:
		"""
		sentieon driver -t ${task.cpus} \\
			-r ${params.genome_file} -i $bam \\
			--algo QualCal ${id}.bqsr.table \\
			-k $params.KNOWN_SITES

		${bqsr_version(task)}
		"""

	stub:
		"""
		touch "${id}.bqsr.table"
		${bqsr_version(task)}
		"""
}
def bqsr_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}

// When rerunning sample
process dnascope {
	cpus 54
	memory '100 GB'
	// 12 GB peak giab //
	time '4h'
	tag "$id"
	container  "${params.container_sentieon}"

	input:
		tuple val(group), val(id), val(bam), val(bai), val(bqsr)

	output:
		tuple val(group), val(id), path("${id}.dnascope.gvcf.gz"), path("${id}.dnascope.gvcf.gz.tbi"), emit: gvcf_tbi
		path "*versions.yml", emit: versions

	when:
		params.varcall

	// TODO: Move shards out to some config file:
	script:
		"""
		sentieon driver \\
			-t ${task.cpus} \\
			-r ${params.genome_file} \\
			-q $bqsr \\
			-i $bam \\
			--shard 1:1-248956422  \\
			--shard 2:1-242193529  \\
			--shard 3:1-198295559  \\
			--shard 4:1-190214555  \\
			--shard 5:1-120339935  \\
			--shard 5:120339936-181538259  \\
			--shard 6:1-170805979  \\
			--shard 7:1-159345973  \\
			--shard 8:1-145138636  \\
			--shard 9:1-138394717  \\
			--shard 10:1-133797422  \\
			--shard 11:1-135086622  \\
			--shard 12:1-56232327  \\
			--shard 12:56232328-133275309  \\
			--shard 13:1-114364328  \\
			--shard 14:1-107043718  \\
			--shard 15:1-101991189  \\
			--shard 16:1-90338345  \\
			--shard 17:1-83257441  \\
			--shard 18:1-80373285  \\
			--shard 19:1-58617616  \\
			--shard 20:1-64444167  \\
			--shard 21:1-46709983  \\
			--shard 22:1-50818468  \\
			--shard X:1-124998478  \\
			--shard X:124998479-156040895  \\
			--shard Y:1-57227415  \\
			--shard M:1-16569 \\
			--algo DNAscope --emit_mode GVCF ${id}.dnascope.gvcf.gz

		${dnascope_version(task)}
		"""

	stub:
		"""
		touch "${id}.dnascope.gvcf.gz"
		touch "${id}.dnascope.gvcf.gz.tbi"

		${dnascope_version(task)}
		"""
}
def dnascope_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}


//Collect various QC data:
process sentieon_qc {
	cpus 52
	memory '30 GB'
	tag "$id"
	time '2h'
	container  "${params.container_sentieon}"

	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple (
			val(group),
			val(id),
			path("mq_metrics.txt"),
			path("qd_metrics.txt"),
			path("gc_summary.txt"),
			path("gc_metrics.txt"),
			path("aln_metrics.txt"),
			path("is_metrics.txt"),
			path("assay_metrics.txt"),
			path("cov_metrics.txt"),
			path("cov_metrics.txt.sample_summary"),
			emit: sentieon_qc_metrics
		)

		path "*versions.yml", emit: versions

	script:
		target = ""
		// A bit of cheating here - these are really optional arguments
		panel_command = "touch cov_metrics.txt cov_metrics.txt.sample_summary"
		cov = "WgsMetricsAlgo assay_metrics.txt"

		if (params.onco || params.exome) {
			target = "--interval $params.intervals"
			cov = "CoverageMetrics --cov_thresh 1 --cov_thresh 10 --cov_thresh 30 --cov_thresh 100 --cov_thresh 250 --cov_thresh 500 cov_metrics.txt"
			panel_command = "sentieon driver -r ${params.genome_file} -t ${task.cpus} -i ${bam} --algo HsMetricAlgo --targets_list ${params.intervals} --baits_list ${params.intervals} assay_metrics.txt"
		}

		"""
		sentieon driver \\
			-r ${params.genome_file} $target \\
			-t ${task.cpus} \\
			-i $bam \\
			--algo MeanQualityByCycle mq_metrics.txt \\
			--algo QualDistribution qd_metrics.txt \\
			--algo GCBias --summary gc_summary.txt gc_metrics.txt \\
			--algo AlignmentStat aln_metrics.txt \\
			--algo InsertSizeMetricAlgo is_metrics.txt \\
			--algo $cov
		$panel_command

		${sentieon_qc_version(task)}
		"""

	stub:
		"""
		touch "assay_metrics.txt"
		touch "mq_metrics.txt"
		touch "qd_metrics.txt"
		touch "gc_summary.txt"
		touch "gc_metrics.txt"
		touch "aln_metrics.txt"
		touch "is_metrics.txt"
		touch "cov_metrics.txt"
		touch "cov_metrics.txt.sample_summary"
		${sentieon_qc_version(task)}
		"""
}
def sentieon_qc_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}

process sentieon_qc_postprocess {
	cpus 2
	memory '1 GB'
	tag "$id"
	time '2h'

	input:
		tuple(
			val(group),
			val(id),
			path(mq_metrics),
			path(qd_metrics),
			path(gc_summary),
			path(gc_metrics),
			path(aln_metrics),
			path(is_metrics),
			path(assay_metrics),
			path(cov_metrics),
			path(cov_metrics_sample_summary),
			path(dedup_metrics)
		)
	output:
		tuple val(group), val(id), path("${id}_qc.json"), emit: qc_json

	script:
		assay = (params.onco || params.exome) ? "panel" : "wgs"
		"""
		qc_sentieon.pl \\
			--SID ${id} \\
			--type ${assay} \\
			--align_metrics_file ${aln_metrics} \\
			--insert_file ${is_metrics} \\
			--dedup_metrics_file ${dedup_metrics} \\
			--metrics_file ${assay_metrics} \\
			--gcsummary_file ${gc_summary} \\
			--coverage_file ${cov_metrics} \\
			--coverage_file_summary ${cov_metrics_sample_summary} \\
			> ${id}_qc.json

		"""
}

process d4_coverage {
	cpus 16
	memory '10 GB'
	publishDir "${params.results_output_dir}/cov", mode: 'copy', overwrite: 'true', pattern: '*.d4'
	tag "$id"
	container  "${params.container_d4tools}"

	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${id}_coverage.d4"), emit: ch_final_d4
		tuple val(group), path("${group}_d4.INFO"), emit: d4_INFO
		path "*versions.yml", emit: versions

	when:
		params.run_chanjo2

	script:
	"""
	d4tools create \\
		--threads ${task.cpus} \\
		"${bam}" \\
		"${id}_coverage.d4"

	echo "D4	$id	/access/${params.subdir}/cov/${id}_coverage.d4" > ${group}_d4.INFO

	${d4_coverage_version(task)}
	"""

	stub:
	"""
	touch "${id}_coverage.d4"
	touch "${group}_d4.INFO"

	${d4_coverage_version(task)}
	"""
}
def d4_coverage_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    d4tools: \$(echo \$( d4tools 2>&1 | head -1 ) | sed "s/.*version: //" | sed "s/)//" )
	END_VERSIONS
	"""
}

process verifybamid2 {
	cpus 16
	memory '10 GB'
	publishDir "${params.results_output_dir}/contamination", mode: 'copy', overwrite: 'true', pattern: '*.selfSM'
	tag "$id"
	container  "${params.container_verifybamid2}"

	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		path("${id}.result.selfSM")
		path("${id}.result.Ancestry")
		path "*versions.yml", emit: versions

	script:
		"""
		verifybamid2 \
			--SVDPrefix ${params.verifybamid2_svdprefix} \
			--Reference ${params.genome_file} \
			--BamFile ${bam}

			mv result.selfSM ${id}.result.selfSM
			mv result.Ancestry ${id}.result.Ancestry

		${verifybamid2_version(task)}
		"""

	stub:
		"""
		touch "${id}.result.selfSM"
		touch "${id}.result.Ancestry"

		${verifybamid2_version(task)}
		"""
}
def verifybamid2_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    VerifyBamID2: \$( echo \$( verifybamid2 --help 2>&1 | grep Version ) | sed "s/^.*Version://" )
	END_VERSIONS
	"""
}

// Calculate coverage for paneldepth
process depth_onco {
	cpus 2
	time '1h'
	memory '10 GB'
	publishDir "${params.results_output_dir}/cov", mode: 'copy', overwrite: 'true'
	tag "$id"
	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		path("${id}.lowcov.overlapping.bed"), emit: cov_onco

	script:
		"""
		panel_depth.pl $bam $params.scoutbed > ${id}.lowcov.bed
		overlapping_genes.pl ${id}.lowcov.bed $params.gene_regions > ${id}.lowcov.overlapping.bed
		"""

	stub:
		"""
		touch "${id}.lowcov.overlapping.bed"
		"""
}

process SMNCopyNumberCaller {
	cpus 10
	memory '25GB'
	time '2h'
	publishDir "${params.results_output_dir}/plots/SMNcnc", mode: 'copy' , overwrite: 'true', pattern: '*.pdf*'
	tag "$id"

	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		path("*.tsv"), emit: smn_tsv
		tuple path("*.pdf"), path("*.json")
		tuple val(group), path("${group}_smn.INFO"), emit: smn_INFO
		path "*versions.yml", emit: versions

	// TODO: split off the plotting into own process?
	script:
		"""
		samtools view -H $bam | \\
			sed -e 's/SN:1/SN:chr1/' | sed -e 's/SN:2/SN:chr2/' |  \\
			sed -e 's/SN:3/SN:chr3/' | sed -e 's/SN:4/SN:chr4/' |  \\
			sed -e 's/SN:5/SN:chr5/' | sed -e 's/SN:6/SN:chr6/' |  \\
			sed -e 's/SN:7/SN:chr7/' | sed -e 's/SN:8/SN:chr8/' |  \\
			sed -e 's/SN:9/SN:chr9/' | sed -e 's/SN:10/SN:chr10/' | \\
			sed -e 's/SN:11/SN:chr11/' | sed -e 's/SN:12/SN:chr12/' |  \\
			sed -e 's/SN:13/SN:chr13/' | sed -e 's/SN:14/SN:chr14/' |  \\
			sed -e 's/SN:15/SN:chr15/' | sed -e 's/SN:16/SN:chr16/' |  \\
			sed -e 's/SN:17/SN:chr17/' | sed -e 's/SN:18/SN:chr18/' |  \\
			sed -e 's/SN:19/SN:chr19/' | sed -e 's/SN:20/SN:chr20/' |  \\
			sed -e 's/SN:21/SN:chr21/' | sed -e 's/SN:22/SN:chr22/' |  \\
			sed -e 's/SN:X/SN:chrX/' | sed -e 's/SN:Y/SN:chrY/' |   \\
			sed -e 's/SN:MT/SN:chrM/' | \\
			samtools reheader - $bam > ${id}.bam
		samtools index -b ${id}.bam -@ ${task.cpus}
		echo ${id}.bam > manifest.txt
		smn_caller.py --manifest manifest.txt --genome 38 --prefix ${id} --outDir . --threads ${task.cpus}
		rm ${id}.bam
		source activate py3-env
		python /SMNCopyNumberCaller/smn_charts.py -s ${id}.json -o .
		mv ${id}.tsv ${group}_SMN.tsv
		echo "SMN ${params.accessdir}/smn/${group}_SMN.tsv" > ${group}_smn.INFO

		${smn_copy_number_caller_version(task)}
		"""

	stub:
		"""
		touch "${id}.bam"
		touch "${id}.tsv"
		touch "${id}.pdf"
		touch "${id}.json"
		touch "${group}_SMN.tsv"
		touch "${group}_smn.INFO"

		${smn_copy_number_caller_version(task)}
		"""
}
def smn_copy_number_caller_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
	    smn-copy-number-caller: 1.1.2
	END_VERSIONS
	"""
}


////////////////////////////////////////////////////////////////////////
////////////////////////// EXPANSION HUNTER ////////////////////////////
////////////////////////////////////////////////////////////////////////

// call STRs using ExpansionHunter, and plot alignments with GraphAlignmentViewer
process expansionhunter {
	tag "$group"
	cpus 2
	time '10h'
	memory '40 GB'

	input:
		tuple val(group), val(id), path(bam), path(bai), val(sex), val(type) // TODO: sex + type are not needed here


	output:
		tuple val(group), val(id), path("${group}.eh.vcf"), emit: expansionhunter_vcf
		tuple val(group), val(id), path("${group}.eh_realigned.sort.bam"), path("${group}.eh_realigned.sort.bam.bai"), path("${group}.eh.vcf"), emit: bam_vcf
		path "*versions.yml", emit: versions

	when:
		params.str

	script:
		"""
		source activate htslib10
		ExpansionHunter \
			--reads $bam \
			--reference ${params.genome_file} \
			--variant-catalog $params.expansionhunter_catalog \
			--output-prefix ${group}.eh
		samtools sort ${group}.eh_realigned.bam -o ${group}.eh_realigned.sort.bam
		samtools index ${group}.eh_realigned.sort.bam

		${expansionhunter_version(task)}
		"""

	stub:
		"""
		source activate htslib10
		touch "${group}.eh.vcf"
		touch "${group}.eh_realigned.sort.bam"
		touch "${group}.eh_realigned.sort.bam.bai"

		${expansionhunter_version(task)}
		"""
}
def expansionhunter_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    expansionhunter: \$(echo \$(ExpansionHunter --version 2>&1) | sed 's/.*ExpansionHunter v// ; s/]//')
	    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
	END_VERSIONS
	"""
}

// annotate expansionhunter vcf
process stranger {
	tag "$group"
	memory '1 GB'
	time '10m'
	cpus 2
	container  "${params.container_stranger}"

	input:
		tuple val(group), val(id), path(eh_vcf)

	output:
		tuple val(group), val(id), path("${group}.fixinfo.eh.stranger.vcf"), emit: vcf_annotated
		path "*versions.yml", emit: versions

	script:
		"""
		stranger ${eh_vcf} -f $params.expansionhunter_catalog > ${group}.eh.stranger.vcf
		grep ^# ${group}.eh.stranger.vcf > ${group}.fixinfo.eh.stranger.vcf
		grep -v ^# ${group}.eh.stranger.vcf | sed 's/ /_/g' >> ${group}.fixinfo.eh.stranger.vcf

		${stranger_version(task)}
		"""

	stub:
		"""
		touch "${group}.fixinfo.eh.stranger.vcf"
		${stranger_version(task)}
		"""
}
def stranger_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    stranger: \$( stranger --version )
	END_VERSIONS
	"""
}


process reviewer {
	tag "$group"
	cpus 2
	time '1h'
	memory '1 GB'
	errorStrategy 'ignore'
	container  "${params.container_reviewer}"
	publishDir "${params.results_output_dir}/plots/reviewer/${group}", mode: 'copy' , overwrite: 'true', pattern: '*.svg'

	input:
		tuple val(group), val(id), path(bam), path(bai), path(vcf)

	output:
		path("*svg")
		path "*versions.yml", emit: versions

	script:
		version_str = reviewer_version(task)
		"""
		grep LocusId ${params.expansionhunter_catalog} | sed 's/[",^ ]//g' | cut -d':' -f2 | perl -na -e 'chomp; \
		system("REViewer --reads ${bam} \
			--vcf ${vcf} \
			--reference ${params.genome_file} \
			--catalog ${params.expansionhunter_catalog} \
			--locus \$_ \
			--output-prefix ${id}");'

		echo "${version_str}" > "${task.process}_versions.yml"
		"""

	stub:
		version_str = reviewer_version(task)
		"""
		touch "${id}.svg"
		echo "${version_str}" > "${task.process}_versions.yml"
		"""
}
def reviewer_version(task) {
	// TODO: Reconcile this version stub with others.
	"""${task.process}:
	    reviewer: \$(echo \$(REViewer --version 2>&1) | sed 's/^.*REViewer v//')"""
}

// split multiallelic sites in expansionhunter vcf
// FIXME: Use env variable for picard path...
process vcfbreakmulti_expansionhunter {
	publishDir "${params.results_output_dir}/vcf", mode: 'copy' , overwrite: 'true', pattern: '*.vcf.gz'
	tag "$group"
	time '1h'
	memory '50 GB'

	input:
		tuple val(group), val(id), path(eh_vcf_anno)
		tuple val(group2), val(id2), val(sex), val(mother), val(father), val(phenotype), val(diagnosis), val(type), val(assay), val(clarity_sample_id), val(ffpe), val(analysis)

	output:
		path("${group}.expansionhunter.vcf.gz"), emit: expansionhunter_scout
		tuple val(group), path("${group}_str.INFO"), emit: str_INFO
		path "*versions.yml", emit: versions

	script:
		if (father == "") { father = "null" }
		if (mother == "") { mother = "null" }
		if (params.mode == "family") {
			"""
			java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf INPUT=${eh_vcf_anno} OUTPUT=${eh_vcf_anno}.rename.vcf NEW_SAMPLE_NAME=${id}
			vcfbreakmulti ${eh_vcf_anno}.rename.vcf > ${group}.expansionhunter.vcf.tmp
			familyfy_str.pl --vcf ${group}.expansionhunter.vcf.tmp --mother $mother --father $father --out ${group}.expansionhunter.vcf
			bgzip ${group}.expansionhunter.vcf
			tabix ${group}.expansionhunter.vcf.gz
			echo "STR	${params.accessdir}/vcf/${group}.expansionhunter.vcf.gz" > ${group}_str.INFO

			${vcfbreakmulti_expansionhunter_version(task)}
			"""
		}
		else {
			"""
			java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf INPUT=${eh_vcf_anno} OUTPUT=${eh_vcf_anno}.rename.vcf NEW_SAMPLE_NAME=${id}
			vcfbreakmulti ${eh_vcf_anno}.rename.vcf > ${group}.expansionhunter.vcf
			bgzip ${group}.expansionhunter.vcf
			tabix ${group}.expansionhunter.vcf.gz
			echo "STR	${params.accessdir}/vcf/${group}.expansionhunter.vcf.gz" > ${group}_str.INFO

			${vcfbreakmulti_expansionhunter_version(task)}
			"""

		}

	stub:
		"""
		touch "${group}.expansionhunter.vcf.gz"
		touch "${group}_str.INFO"

		${vcfbreakmulti_expansionhunter_version(task)}
		"""
}
def vcfbreakmulti_expansionhunter_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcflib: 1.0.9
	    rename-sample-in-vcf: \$(echo \$(java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf --version 2>&1) | sed 's/-SNAPSHOT//')
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	END_VERSIONS
	"""
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


// MELT always give VCFs for each type of element defined in mei_list
// If none found -> 0 byte vcf. merge_melt.pl merges the three, if all empty
// it creates a vcf with only header
// merge_melt.pl gives output ${id}.melt.merged.vcf
process melt {
	cpus 3
	errorStrategy 'retry'
	container  "${params.container_melt}"
	tag "$id"
	// memory seems to scale with less number of reads?
	memory '70 GB'
	time '3h'
	publishDir "${params.results_output_dir}/vcf", mode: 'copy' , overwrite: 'true', pattern: '*.vcf'

	input:
		tuple val(group), val(id), path(bam), path(bai)
		tuple val(group1), val(id1), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV)

	output:
		tuple val(group), val(id), path("${id}.melt.merged.vcf"), emit: melt_vcf_nonfiltered
		path "*versions.yml", emit: versions

	script:
		"""
		java -jar /opt/MELTv2.2.2/MELT.jar Single \\
			-bamfile $bam \\
			-r 150 \\
			-h ${params.genome_file} \\
			-n /opt/MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed \\
			-z 500000 \\
			-d 50 \\
			-t $params.mei_list \\
			-w . \\
			-c $MEAN_DEPTH \\
			-e $INS_SIZE \\
			-exome
		merge_melt.pl $params.meltheader $id

		${melt_version(task)}
		"""

	stub:
		"""
		touch "${id}.melt.merged.vcf"
		${melt_version(task)}
		"""
}
def melt_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    melt: \$(echo \$(java -jar /opt/MELTv2.2.2/MELT.jar -h | grep "^MELTv" | cut -f1 -d" " | sed "s/MELTv//" ) )
	END_VERSIONS
	"""
}

process intersect_melt {
	cpus 2
	tag "$id"
	memory '2 GB'
	time '1h'
	publishDir "${params.results_output_dir}/vcf", mode: 'copy' , overwrite: 'true', pattern: '*.vcf'

	input:
		tuple val(group), val(id), path(vcf)

	output:
		tuple val(group), val(id), path("${id}.melt.merged.intersected.vcf"), emit: merged_intersected_vcf
		path "*versions.yml", emit: versions

	when:
		params.run_melt

	script:
		"""
		bedtools intersect -a $vcf -b $params.intersect_bed -header > ${id}.melt.merged.intersected.vcf
		${intersect_melt_version(task)}
		"""

	stub:
		"""
		touch "${id}.melt.merged.intersected.vcf"
		${intersect_melt_version(task)}
		"""
}
def intersect_melt_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bedtools: \$(echo \$(bedtools --version 2>&1) | sed -e "s/^.*bedtools v//" )
	END_VERSIONS
	"""
}


process bamtoyaml {
	cpus 1
	time "5m"
	memory "2MB"

	input:
	tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), path("${group}_bamstart.INFO"), emit: bamchoice_INFO

	script:
		"""
		echo "BAM	$id	/access/${params.subdir}/bam/${bam.getName()}" > ${group}_bamstart.INFO
		"""

	stub:
		"""
		touch "${group}_bamstart.INFO"
		"""
}


process gvcf_combine {
	cpus 16
	tag "$group"
	memory '5 GB'
	time '5h'
	container  "${params.container_sentieon}"

	input:
		tuple val(group), val(id), path(gvcfs), path(gvcf_idxs)

	output: // Off to split_normalize, together with other stuff
		tuple val(group), val(id), path("${group}.combined.vcf.gz"), path("${group}.combined.vcf.gz.tbi"), emit: combined_vcf
		path "*versions.yml", emit: versions

	script:
		all_gvcfs = gvcfs.collect { it.toString() }.sort().join(' -v ')
		"""
		sentieon driver \\
			-t ${task.cpus} \\
			-r ${params.genome_file} \\
			--algo GVCFtyper \\
			-v $all_gvcfs ${group}.combined.vcf.gz

		${gvcf_combine_version(task)}
		"""

	stub:
		all_gvcfs = gvcfs.collect { it.toString() }.sort().join(' -v ')
		"""
		touch "${group}.combined.vcf.gz"
		touch "${group}.combined.vcf.gz.tbi"

		${gvcf_combine_version(task)}
		"""
}
def gvcf_combine_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}

// Create ped
process create_ped {
	tag "$group"
	time '20m'
	publishDir "${params.results_output_dir}/ped", mode: 'copy' , overwrite: 'true'
	memory '1 GB'

	input:
		tuple val(group), val(id), val(type), val(sex), val(mother), val(father)

	output:
		tuple val(group), val(type), path("${group}_base.ped"), emit: ped_base
		tuple val(group), val(type_ma), path("${group}_ma.ped"), emit: ped_ma, optional: true
		tuple val(group), val(type_fa), path("${group}_fa.ped"), emit: ped_fa, optional: true

	script:
		if ( father == "" ) {
			father = "0"
		}
		if ( mother == "" ) {
			mother = "0"
		}
		type_fa = "fa"
		type_ma = "ma"
		"""
		create_ped.pl --mother $mother --father $father --group $group --id $id --sex $sex
		"""

	stub:
		// TODO: _ma and _fa.ped stub files should only be created for trios.
		//       otherwise risk for messing up wgs-single/panel stub runs.
		type_fa = "fa"
		type_ma = "ma"
		"""
		touch "${group}_base.ped"
		touch "${group}_ma.ped"
		touch "${group}_fa.ped"

        echo $type_fa $type_ma > type.val
		"""
}

//madeline ped, run if family mode
process madeline {
	publishDir "${params.results_output_dir}/ped", mode: 'copy' , overwrite: 'true', pattern: '*.xml'
	memory '1 GB'
	time '1h'
	cpus 2
	container  "${params.container_madeline}"

	input:
		tuple val(group), val(type), path(ped)

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
		echo "MADDE	$type ${params.accessdir}/ped/${ped}.madeline.xml" > ${group}_madde.INFO

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

process freebayes {
	cpus 1
	time '2h'
	memory '10 GB'
	container  "${params.container_twist_myeloid}"

	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), path("${id}.pathfreebayes.vcf_no_header.tsv.gz"), emit: freebayes_variants
		path "*versions.yml", emit: versions


	when:
		params.antype == "panel"


	script:
		if (params.onco) {
			"""
			freebayes -f ${params.genome_file} --pooled-continuous --pooled-discrete -t $params.intersect_bed --min-repeat-entropy 1 -F 0.03 $bam > ${id}.freebayes.vcf
			vcfbreakmulti ${id}.freebayes.vcf > ${id}.freebayes.multibreak.vcf
			bcftools norm -m-both -c w -O v -f ${params.genome_file} -o ${id}.freebayes.multibreak.norm.vcf ${id}.freebayes.multibreak.vcf
			vcfanno_linux64 -lua $params.VCFANNO_LUA $params.vcfanno ${id}.freebayes.multibreak.norm.vcf > ${id}.freebayes.multibreak.norm.anno.vcf
			grep ^# ${id}.freebayes.multibreak.norm.anno.vcf > ${id}.freebayes.multibreak.norm.anno.path.vcf
			grep -v ^# ${id}.freebayes.multibreak.norm.anno.vcf | grep -i pathogenic > ${id}.freebayes.multibreak.norm.anno.path.vcf2
			cat ${id}.freebayes.multibreak.norm.anno.path.vcf ${id}.freebayes.multibreak.norm.anno.path.vcf2 > ${id}.freebayes.multibreak.norm.anno.path.vcf3
			filter_freebayes.pl ${id}.freebayes.multibreak.norm.anno.path.vcf3 | bgzip -c > "${id}.pathfreebayes.vcf_no_header.tsv.gz"

			${freebayes_version(task)}
			"""
		}
		else {
			"""
			touch "${id}.pathfreebayes.vcf_no_header.tsv.gz"

			${freebayes_version(task)}
			"""
		}

	stub:
		"""
		touch "${id}.pathfreebayes.vcf_no_header.tsv.gz"

		${freebayes_version(task)}
		"""
}
def freebayes_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g')
	    vcflib: 1.0.9
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	    vcfanno: \$(echo \$(vcfanno_linux64 2>&1 | grep version | cut -f3 -d' ')  )
	END_VERSIONS
	"""
}

/////////////// MITOCHONDRIA SNV CALLING ///////////////
///////////////                          ///////////////

// create an MT BAM file
process fetch_MTseqs {
	cpus 2
	memory '10GB'
	time '1h'
	tag "$id"
	publishDir "${params.results_output_dir}/bam", mode: 'copy', overwrite: 'true', pattern: '*.bam*'

	input:
		tuple val(group), val(id), path(bam), path(bai)

    output:
        tuple val(group), val(id), file ("${id}_mito.bam"), path("${id}_mito.bam.bai"), emit: bam_bai
		tuple val(group), path("${group}_mtbam.INFO"), emit: mtBAM_INFO
		path "*versions.yml", emit: versions

	script:
		"""
		sambamba view -f bam $bam M > ${id}_mito.bam
		samtools index -b ${id}_mito.bam
		echo "mtBAM	$id	/access/${params.subdir}/bam/${id}_mito.bam" > ${group}_mtbam.INFO

		${fetch_MTseqs_version(task)}
		"""

	stub:
		"""
		touch "${id}_mito.bam"
		touch "${id}_mito.bam.bai"
		touch "${group}_mtbam.INFO"

		${fetch_MTseqs_version(task)}
		"""
}
def fetch_MTseqs_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sambamba: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
	    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
	END_VERSIONS
	"""
}


process sentieon_mitochondrial_qc {

    // Fetch mitochondrial coverage statistics
    // Calculate mean_coverage and pct_above_500x

    cpus 30
    memory '20 GB'
	tag "$id"
	time '2h'
	container  "${params.container_sentieon}"

	input:
        tuple val(group), val(id), path(bam), path(bai)

	output:
    	tuple val(group), val(id), path("${id}_mito_coverage.tsv"), emit: qc_tsv
		path "*versions.yml", emit: versions

	when:
	    params.antype == "wgs"

	script:
		"""
		sentieon driver \\
			-r ${params.genome_file} \\
			-t ${task.cpus} \\
			-i $bam \\
			--algo CoverageMetrics \\
			--omit_base_output  \\
			--omit_locus_stat \\
			--omit_sample_stat \\
			--cov_thresh 500 \\
			mt_cov_metrics.txt

		head -1 mt_cov_metrics.txt.sample_interval_summary > "${id}_mito_coverage.tsv"
		grep "^M" mt_cov_metrics.txt.sample_interval_summary >> "${id}_mito_coverage.tsv"
		${sentieon_mitochondrial_qc_version(task)}
		"""

	stub:
		"""
		touch "${id}_mito_coverage.tsv"
		${sentieon_mitochondrial_qc_version(task)}
		"""
}
def sentieon_mitochondrial_qc_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}

process build_mitochondrial_qc_json {
    memory '1 GB'
    cpus 2
    tag "$id"
    time "1h"

    input:
        tuple val(group), val(id), path(mito_qc_file)
    output:
        tuple val(group), val(id), path("${id}_mito_qc.json"), emit: qc_json

	script:
		"""
		mito_tsv_to_json.py ${mito_qc_file} > "${id}_mito_qc.json"
		"""
	stub:
		"""
		touch "${id}_mito_qc.json"
		"""
}


// gatk FilterMutectCalls in future if FPs overwhelms tord/sofie/carro
process run_mutect2 {
	cpus 4
	memory '50 GB'
	time '1h'
	tag "$group"
	publishDir "${params.results_output_dir}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf'

	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${group}.mutect2.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	when:
		!params.onco

	script:
		bams = bam.join(' -I ')

		"""
		source activate gatk4-env
		gatk Mutect2 \
		--mitochondria-mode \
		-R $params.genome_file \
		-L M \
		-I $bams \
		-O ${group}.mutect2.vcf

		${run_mutect2_version(task)}
		"""

	stub:
		bams = bam.join(' -I ')
		"""
		source activate gatk4-env
		touch "${group}.mutect2.vcf"

		${run_mutect2_version(task)}
		"""
}
def run_mutect2_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
	END_VERSIONS
	"""
}

// split and left-align variants
process split_normalize_mito {
	cpus 2
	memory '1GB'
	time '1h'

	input:
		tuple val(group), val(id), path(mito_snv_vcf)
		tuple val(group2), val(proband_id), val(sex), val(type)

	output:
		tuple val(group), path("${group}.mutect2.breakmulti.filtered5p.0genotyped.proband.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		# Old workaround to remove false-positive that crashes bcftools norm:
		# TODO: deal w/ this in some better way
		grep -vP "^M\\s+955" ${mito_snv_vcf} > ${mito_snv_vcf}.fix

		bcftools norm -m-both -o ${mito_snv_vcf}.breakmulti ${mito_snv_vcf}.fix
		bcftools sort ${mito_snv_vcf}.breakmulti | bgzip > ${mito_snv_vcf}.breakmulti.fix
		tabix -p vcf ${mito_snv_vcf}.breakmulti.fix
		bcftools norm -f $params.rCRS_fasta -o ${mito_snv_vcf.baseName}.adjusted.vcf ${mito_snv_vcf}.breakmulti.fix
		bcftools view -i 'FMT/AF[*]>0.05' ${mito_snv_vcf.baseName}.adjusted.vcf -o ${group}.mutect2.breakmulti.filtered5p.vcf
		bcftools filter -S 0 --exclude 'FMT/AF[*]<0.05' ${group}.mutect2.breakmulti.filtered5p.vcf -o ${group}.mutect2.breakmulti.filtered5p.0genotyped.vcf
		filter_mutect2_mito.pl ${group}.mutect2.breakmulti.filtered5p.0genotyped.vcf ${proband_id} > ${group}.mutect2.breakmulti.filtered5p.0genotyped.proband.vcf

		${split_normalize_mito_version(task)}
		"""

	stub:
		"""
		echo "${proband_id}" > proband.id
		touch "${group}.mutect2.breakmulti.filtered5p.0genotyped.proband.vcf"
		${split_normalize_mito_version(task)}
		"""
}
def split_normalize_mito_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcflib: 1.0.9
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	END_VERSIONS
	"""
}

// use python tool HmtNote for annotating vcf
// future merging with diploid genome does not approve spaces in info-string
// TODO: what is this future merging issue and does it still apply?
process run_hmtnote {
	cpus 2
	memory '5GB'
	time '1h'

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.fixinfo.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		source activate tools
		hmtnote annotate ${vcf} ${group}.hmtnote --offline
		grep ^# ${group}.hmtnote > ${group}.fixinfo.vcf
		grep -v ^# ${group}.hmtnote | sed 's/ /_/g' >> ${group}.fixinfo.vcf

		${run_hmtnote_version(task)}
		"""

	stub:
		"""
		source activate tools
		touch "${group}.fixinfo.vcf"

		${run_hmtnote_version(task)}
		"""
}
def run_hmtnote_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    hmtnote: \$(echo \$(hmtnote --version 2>&1) | sed 's/^.*hmtnote, version //; s/Using.*\$//' )
	END_VERSIONS
	"""
}

// run haplogrep 2 on resulting vcf
process run_haplogrep {
	time '1h'
	memory '50 GB'
	cpus 2
	publishDir "${params.results_output_dir}/plots/mito", mode: 'copy', overwrite: 'true', pattern: '*.png'

	input:
		tuple val(group), val(id), path(mito_snv_vcf)

	output:
		path "${group}.haplogrep.png"
		tuple val(group), path("${group}_haplo.INFO"), emit: haplogrep_INFO
		path "*versions.yml", emit: versions

	script:
		version_str = run_haplogrep_version(task)
		"""
		for sample in \$(bcftools query -l "${mito_snv_vcf}"); do

			bcftools view -c1 -Oz -s "\$sample" -o "\${sample}.vcf.gz" "${mito_snv_vcf}"
			java  -Xmx16G -Xms16G -jar /opt/bin/haplogrep.jar classify \
			--in "\${sample}.vcf.gz" \\
			--out "\${sample}.hg2.vcf" \\
			--format vcf \\
			--lineage 1

			dot "\${sample}.hg2.vcf.dot" -Tps2 > "\${sample}.hg2.vcf.ps"

			gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dGraphicsAlphaBits=4 -r1200 -dDownScaleFactor=3 -sOutputFile=\${sample}.hg2.vcf.png \${sample}.hg2.vcf.ps

		done
		montage -mode concatenate -tile 3x1 *.png ${group}.haplogrep.png
		echo "IMG haplogrep ${params.accessdir}/plots/mito/${group}.haplogrep.png" > "${group}_haplo.INFO"

		echo "${version_str}" > "${task.process}_versions.yml"
		"""

	stub:
		version_str = run_haplogrep_version(task)
		"""
		touch "${group}.haplogrep.png"
		touch "${group}_haplo.INFO"

		echo "${version_str}" > "${task.process}_versions.yml"
		"""
}
def run_haplogrep_version(task) {
	// TODO: Reconcile this version stub with others.
	"""${task.process}:
	    haplogrep: \$(echo \$(java -jar /opt/bin/haplogrep.jar classify 2>&1) | sed "s/htt.*Classify v// ; s/ .*//")
	    montage: \$(echo \$(gm -version 2>&1) | head -1 | sed -e "s/GraphicsMagick //" | cut -d" " -f1 )"""
}

// use eKLIPse for detecting mitochondrial deletions
process run_eklipse {

	tag "$id"
	cpus 2
	// in rare cases with samples above 50 000x this can peak at 500+ GB of VMEM. Add downsampling!
	memory '100GB'
	time '60m'
	publishDir "${params.results_output_dir}/plots/mito", mode: 'copy', overwrite: 'true', pattern: '*.txt'
	publishDir "${params.results_output_dir}/plots/mito", mode: 'copy', overwrite: 'true', pattern: '*.png'

	input:
		tuple val(group), val(id), val(sex), val(type), path(bam), path(bai)

	output:
		tuple path("*.png"), path("${id}.hetplasmid_frequency.txt")
		tuple val(group), path("${id}_eklipse.INFO"), emit: eklipse_INFO, optional: true
		path "*versions.yml", emit: versions

	script:
		yml_info_command = ""
		if (type == "proband") {
			yml_info_command = "echo 'IMG eklipse ${params.accessdir}/plots/mito/${id}_eklipse.png' > ${id}_eklipse.INFO"
		}
		"""
		source activate htslib10
		echo "${bam}\tsample" > infile.txt
		python /eKLIPse/eKLIPse.py \
		-in infile.txt \
		-ref /eKLIPse/data/NC_012920.1.gb
		mv eKLIPse_*/eKLIPse_deletions.csv ./${id}_deletions.csv
		mv eKLIPse_*/eKLIPse_genes.csv ./${id}_genes.csv
		mv eKLIPse_*/eKLIPse_sample.png ./${id}_eklipse.png
		hetplasmid_frequency_eKLIPse.pl --bam ${bam} --in ${id}_deletions.csv
		mv hetplasmid_frequency.txt ${id}.hetplasmid_frequency.txt
		$yml_info_command

		${run_eklipse_version(task)}
		"""

	stub:
		yml_info_command = ""
		if (type == "proband") {
			yml_info_command = "echo 'IMG eklipse ${params.accessdir}/plots/mito/${id}_eklipse.png' > ${id}_eklipse.INFO"
		}
		"""
		source activate htslib10
		touch "${id}.hetplasmid_frequency.txt"
		touch "${id}.png"
		$yml_info_command

		${run_eklipse_version(task)}
		"""
}
def run_eklipse_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    eklipse: 1.8
	END_VERSIONS
	"""
}

//eklipseM_INFO.collectPath(name: "eklipse.INFO").set{ eklipse_INFO }

// Splitting & normalizing variants, merging with Freebayes/Mutect2, intersecting against exome/clinvar introns
process split_normalize {
	cpus 2
	publishDir "${params.results_output_dir}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	tag "$group"
	memory '50 GB'
	time '1h'
	input:
		tuple val(group), val(ids), path(vcf), path(idx) // is ids supposed to be tuple?
		tuple val(group2), path(vcfconcat)

	output:
		tuple val(group), path("${group}.norm.uniq.DPAF.vcf"), emit: norm_uniq_dpaf_vcf
		tuple val(group), val(id), path("${group}.intersected.vcf"), emit: intersected_vcf
		tuple val(group), path("${group}.multibreak.vcf"), emit: multibreak_vcf
		path "*versions.yml", emit: versions

	script:
	id = ids[0]
	// rename M to MT because genmod does not recognize M
	if (params.onco || params.assay == "modycf") {
		"""
		if [[ -s "$vcfconcat" ]]; then
			zcat $vcf $vcfconcat > ${id}.concat.freebayes.vcf
		else
			# Otherwise it crashes for modycf where $vcfconcat is empty
			zcat $vcf > ${id}.concat.freebayes.vcf
		fi
		vcfbreakmulti ${id}.concat.freebayes.vcf > ${group}.multibreak.vcf
		bcftools norm -m-both -c w -O v -f ${params.genome_file} -o ${group}.norm.vcf ${group}.multibreak.vcf
		bcftools sort ${group}.norm.vcf | vcfuniq > ${group}.norm.uniq.vcf
		wgs_DPAF_filter.pl ${group}.norm.uniq.vcf > ${group}.norm.uniq.DPAF.vcf
		bedtools intersect \\
			-a ${group}.norm.uniq.DPAF.vcf \\
			-b ${params.intersect_bed} \\
			-u -header > ${group}.intersected.vcf

		${split_normalize_version(task)}
		"""
	} else {
		"""
		vcfbreakmulti ${vcf} > ${group}.multibreak.vcf
		bcftools norm -m-both -c w -O v -f ${params.genome_file} -o ${group}.norm.vcf ${group}.multibreak.vcf
		bcftools sort ${group}.norm.vcf | vcfuniq > ${group}.norm.uniq.vcf
		wgs_DPAF_filter.pl ${group}.norm.uniq.vcf > ${group}.norm.uniq.DPAF.vcf
		bedtools intersect \\
			-a ${group}.norm.uniq.DPAF.vcf \\
			-b ${params.intersect_bed} \\
			-u -header > ${group}.intersected_diploid.vcf
		java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar MergeVcfs \\
		I=${group}.intersected_diploid.vcf I=${vcfconcat} O=${group}.intersected.vcf
		sed 's/^M/MT/' -i ${group}.intersected.vcf
		sed 's/ID=M,length/ID=MT,length/' -i ${group}.intersected.vcf

		${split_normalize_version(task)}
		"""
	}

	stub:
		id = ids[0]
		"""
		echo ${id} > id.val
		touch "${group}.norm.uniq.DPAF.vcf"
		touch "${group}.intersected.vcf"
		touch "${group}.multibreak.vcf"

		${split_normalize_version(task)}
		"""
}
def split_normalize_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcflib: 1.0.9
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	    bedtools: \$(echo \$(bedtools --version 2>&1) | sed -e "s/^.*bedtools v//" )
	    merge-vcfs: \$(echo \$(java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar MergeVcfs --version 2>&1 | sed 's/-SNAPSHOT//'))
	END_VERSIONS
	"""
}

/////////////// Collect QC, emit: single file ///////////////

process merge_qc_json {
    cpus 2
    errorStrategy 'retry'
    maxErrors 5
    publishDir "${params.results_output_dir}/qc", mode: 'copy' , overwrite: 'true', pattern: '*.QC'
    tag "$id"
    time '1h'
	memory '1 GB'

    input:
        tuple val(group), val(id), path(qc)

    output:
    	tuple val(group), val(id), path("${id}.QC"), emit: qc_cdm_merged

    script:
        qc_json_files = qc.join(' ')
		"""
		merge_json_files.py ${qc_json_files} > ${id}.QC
		"""

	stub:
		"""
		touch "${id}.QC"
		"""
}

// Load QC data, emit: CDM (via middleman)
process qc_to_cdm {
	cpus 2
	errorStrategy 'retry'
	maxErrors 5
	publishDir "${params.crondir}/qc", mode: 'copy' , overwrite: 'true'
	tag "$id"
	time '1h'

	input:
		tuple val(group), val(id), path(qc_json), val(diagnosis), val(r1), val(r2)

	output:
		path("${id}.cdm"), emit: cdm_done


	when:
		!params.noupload


	script:
		parts = r1.split('/')
		idx =  parts.findIndexOf {it ==~ /......_......_...._........../}
		rundir = parts[0..idx].join("/")
		"""
	    echo "--run-folder ${rundir} --sample-id ${id} --subassay ${diagnosis} --assay ${params.cdm_assay} --qc ${params.results_output_dir}/qc/${id}.QC" > ${id}.cdm
		"""
}


process peddy {

	publishDir "${params.results_output_dir}/ped", mode: 'copy' , overwrite: 'true', pattern: '*.ped'
	publishDir "${params.results_output_dir}/ped", mode: 'copy' , overwrite: 'true', pattern: '*.csv'

	cpus 4
	tag "$group"
	time '1h'
	memory '20GB'

	input:
		tuple val(group), val(type), path(vcf), path(idx), path(ped)

	output:
		tuple path("${group}.ped_check.csv"),path("${group}.peddy.ped"), path("${group}.sex_check.csv"), emit: peddy_files
		tuple val(group), path("${group}_peddy.INFO"), emit: peddy_INFO
		path "*versions.yml", emit: versions

	when:
		!params.annotate_only && params.run_peddy

	script:
		"""
		source activate py3-env
		python -m peddy --sites hg38 -p ${task.cpus} $vcf $ped --prefix $group
		echo "PEDDY	${params.accessdir}/ped/${group}.ped_check.csv,${params.accessdir}/ped/${group}.peddy.ped,${params.accessdir}/ped/${group}.sex_check.csv" > ${group}_peddy.INFO

		${peddy_version(task)}
		"""

	stub:
		"""
		source activate py3-env
		touch "${group}.ped_check.csv"
		touch "${group}.peddy.ped"
		touch "${group}.sex_check.csv"
		touch "${group}_peddy.INFO"

		${peddy_version(task)}
		"""
}
def peddy_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    peddy: \$(echo \$(python -m peddy --version 2>&1) | sed 's/^.*peddy, version //')
	END_VERSIONS
	"""
}

// Extract all variants (
process fastgnomad {
	cpus 2
	memory '40 GB'
	tag "$group"
	publishDir "${params.results_output_dir}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	time '2h'

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.SNPs.vcf"), emit: vcf

	when:
		params.antype == "wgs"

	script:
		"""
		gzip -c $vcf > ${vcf}.gz
		annotate -g $params.FASTGNOMAD_REF -i ${vcf}.gz > ${group}.SNPs.vcf
		"""

	stub:
		"""
		touch "${group}.SNPs.vcf"
	"""
}


// Call UPD regions
process upd {
	tag "$group"
	time '1h'
	memory '1 GB'
	cpus 2

	input:
		tuple val(group), path(vcf)
		tuple val(group2), val(id), val(mother), val(father)

	output:
		path("upd.bed"), emit: upd_bed
		tuple val(group), path("upd.sites.bed"), emit: upd_sites
		path "*versions.yml", emit: versions

	script:
		if( params.mode == "family" && params.trio ) {
			"""
			upd --vcf $vcf --proband $id --mother $mother --father $father --af-tag GNOMADAF regions > upd.bed
			upd --vcf $vcf --proband $id --mother $mother --father $father --af-tag GNOMADAF sites > upd.sites.bed

			${upd_version(task)}
			"""
		}
		else {
			"""
			touch "upd.bed"
			touch "upd.sites.bed"

			${upd_version(task)}
			"""
		}

	stub:
		"""
		touch "upd.bed"
		touch "upd.sites.bed"

		${upd_version(task)}
		"""
}
def upd_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    upd: \$(echo \$(upd --version 2>&1))
	END_VERSIONS
	"""
}


process upd_table {
	publishDir "${params.results_output_dir}/plots", mode: 'copy' , overwrite: 'true'
	tag "$group"
	time '1h'
	memory '1 GB'
	cpus 2

	input:
		tuple val(group), path(upd_sites)

	output:
		path("${group}.UPDtable.xls")

	when:
		params.mode == "family" && params.trio

	script:
		"""
		upd_table.pl $upd_sites > ${group}.UPDtable.xls
		"""

	stub:
		"""
		touch "${group}.UPDtable.xls"
		"""
}


// Call ROH regions
process roh {
	tag "$group"
	time '1h'
	memory '1 GB'
	cpus 2

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("roh.txt"), emit: roh_plot
		path "*versions.yml", emit: versions

	script:
		"""
		bcftools roh --rec-rate 1e-9 --AF-tag GNOMADAF ${vcf} -o roh.txt
		${roh_version(task)}
		"""

	stub:
		"""
		touch "roh.txt"
		${roh_version(task)}
		"""
}
def roh_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS
	"""
}

// // Create coverage profile using GATK
process gatkcov {
	publishDir "${params.results_output_dir}/cov", mode: 'copy' , overwrite: 'true', pattern: '*.tsv'
	tag "$group"
	cpus 2
	memory '80 GB'
	time '5h'

	input:
		tuple val(group), val(id), val(sex), val(type), path(bam), path(bai)

	output:
		tuple val(group), val(id), val(type), val(sex), path("${id}.standardizedCR.tsv"), path("${id}.denoisedCR.tsv"), emit: cov_plot
		tuple val(group), val(id), path("${id}.standardizedCR.tsv"), path("${id}.denoisedCR.tsv"), emit: cov_gens
		path "*versions.yml", emit: versions

	script:

		def PON = [F: params.GATK_PON_FEMALE, M: params.GATK_PON_MALE]

		"""
		source activate gatk4-env

		gatk CollectReadCounts \\
			-I $bam -L $params.COV_INTERVAL_LIST \\
			--interval-merging-rule OVERLAPPING_ONLY -O ${bam}.hdf5

		gatk --java-options "-Xmx30g" DenoiseReadCounts \\
			-I ${bam}.hdf5 --count-panel-of-normals ${PON[sex]} \\
			--standardized-copy-ratios ${id}.standardizedCR.tsv \\
			--denoised-copy-ratios ${id}.denoisedCR.tsv

		gatk PlotDenoisedCopyRatios \\
			--standardized-copy-ratios ${id}.standardizedCR.tsv \\
			--denoised-copy-ratios ${id}.denoisedCR.tsv \\
			--sequence-dictionary $params.GENOMEDICT \\
			--minimum-contig-length 46709983 --output . --output-prefix $id

		${gatkcov_version(task)}
		"""

	stub:
		"""
		source activate gatk4-env
		touch "${id}.standardizedCR.tsv"
		touch "${id}.denoisedCR.tsv"

		${gatkcov_version(task)}
		"""
}
def gatkcov_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
	END_VERSIONS
	"""
}


// Plot ROH, UPD and coverage in a genomic overview plot
process overview_plot {

	cpus 2
	tag "$group"
	time '1h'
	memory '5 GB'
	publishDir "${params.results_output_dir}/plots", mode: 'copy' , overwrite: 'true', pattern: "*.png"

	input:
		tuple val(group), val(id), val(type), val(sex), path(cov_stand), path(cov_denoised)
		tuple val(group2), path(roh)
		path(upd_bed)

	output:
		path("${group}.genomic_overview.png")
		tuple val(group), path("${group}_oplot.INFO"), emit: oplot_INFO

	script:
		"""
		genome_plotter.pl --dict ${params.GENOMEDICT} \\
			--sample ${id} \\
			--upd ${upd_bed} \\
			--roh ${roh} \\
			--sex ${sex} \\
			--cov ${cov_denoised} \\
			--out ${group}.genomic_overview.png
		echo "IMG overviewplot	${params.accessdir}/plots/${group}.genomic_overview.png" > ${group}_oplot.INFO
		"""

	stub:
		"""
		touch "${group}.genomic_overview.png"
		touch "${group}_oplot.INFO"
		"""
}

process generate_gens_data {
	publishDir "${params.results_output_dir}/plot_data", mode: 'copy' , overwrite: 'true', pattern: "*.gz*"
	publishDir "${params.crondir}/gens", mode: 'copy', overwrite: 'true', pattern: "*.gens"
	tag "$group"
	cpus 1
	time '3h'
	memory '5 GB'

	input:
		tuple val(group), val(id), path(gvcf), path(gvcf_index), path(cov_stand), path(cov_denoise)

	output:
		tuple path("${id}.cov.bed.gz"), path("${id}.baf.bed.gz"), path("${id}.cov.bed.gz.tbi"), path("${id}.baf.bed.gz.tbi"), path("${id}.overview.json.gz")
		tuple val(group), val(id), val(true), 	emit: is_done
		path("${id}.gens"), 					emit: gens_middleman

	when:
		params.prepare_gens_data

	script:
		"""
		generate_gens_data.pl $cov_stand $gvcf $id $params.GENS_GNOMAD
		echo "gens load sample --sample-id $id --case-id $group --genome-build 38 --baf ${params.gens_accessdir}/${id}.baf.bed.gz --coverage ${params.gens_accessdir}/${id}.cov.bed.gz --overview-json ${params.gens_accessdir}/${id}.overview.json.gz" > ${id}.gens
		"""

	stub:
		"""
		touch "${id}.cov.bed.gz"
		touch "${id}.baf.bed.gz"
		touch "${id}.cov.bed.gz.tbi"
		touch "${id}.baf.bed.gz.tbi"
		touch "${id}.overview.json.gz"
		touch "${id}.gens"
		"""
}

process generate_gens_v4_meta {
	publishDir "${params.results_output_dir}/plot_data", mode: 'copy', overwrite: 'true', pattern: "*.tsv"
	publishDir "${params.results_output_dir}/plot_data", mode: 'copy', overwrite: 'true', pattern: "*.bed"
	tag "$id"
	cpus 1
	time '1h'
	memory '10 GB'

	input:
		tuple val(group), val(id), val(type), val(sex), path(cov_stand), path(cov_denoise), path(roh), path(upd_bed), path(upd_sites)

	output:
		tuple val(group), val(id), val(type), val(sex), path("${id}.gens_track.roh.bed"), path("${id}.gens_track.upd.bed"), path("${id}.meta.tsv"), path("${id}.chrom_meta.tsv"), emit: meta
	
	when:
		params.prepare_gens_data
	
	script:
		"""
		if [ "$type" == "proband" ]; then
			prepare_gens_v4_input.py \\
				--roh ${roh} \\
				--upd_regions ${upd_bed} \\
				--upd_sites ${upd_sites} \\
				--cov ${cov_denoise} \\
				--chrom_lengths ${params.GENOMEDICT} \\
				--sample ${id} \\
				--sex ${sex} \\
				--out_gens_track_roh "${id}.gens_track.roh.bed" \\
				--out_gens_track_upd "${id}.gens_track.upd.bed" \\
				--out_meta "${id}.meta.tsv" \\
				--out_chrom_meta "${id}.chrom_meta.tsv" \\
				--analysis_mode "${params.mode}"
		else
			touch "${id}.gens_track.roh.bed"
			touch "${id}.gens_track.upd.bed"
			touch "${id}.meta.tsv"
			touch "${id}.chrom_meta.tsv"
		fi
		"""

	stub:
		"""
		touch "${id}.gens_track.roh.bed"
		touch "${id}.gens_track.upd.bed"
		touch "${id}.meta.tsv"
		touch "${id}.chrom_meta.tsv"

		touch "${id}.gens"
		"""
}

process gens_v4_cron {
	publishDir "${params.crondir}/gens", mode: 'copy', overwrite: 'true', pattern: "*.gens_v4_const"
	tag "$id"
	cpus 1
	time '10m'
	memory '1 GB'

	input:
		tuple val(group), val(id), val(type), val(sex), path(track_roh), path(track_upd), path(meta_tsv), path(chrom_meta_tsv), val(_gens_input_data_is_done)
	
	output:
		path("${id}.gens_v4_const"), emit: gens_v4_middleman
	
	when:
		params.prepare_gens_data

	script:
		def meta_opts = type == "proband" ? 
			"--meta ${params.gens_accessdir}/${meta_tsv.getName()} --meta ${params.gens_accessdir}/${chrom_meta_tsv.getName()}":
			""
		"""
		echo "gens load sample \\
			--sample-id $id \\
			--case-id $group \\
			--genome-build 38 \\
			--sex $sex \\
			--sample-type $type \\
			--baf ${params.gens_accessdir}/${id}.baf.bed.gz \\
			--coverage ${params.gens_accessdir}/${id}.cov.bed.gz \\
			${meta_opts}" > ${id}.gens_v4_const

		if [[ "$type" == "proband" ]]; then
			echo "gens load sample-annotation \\
				--sample-id $id \\
				--case-id $group \\
				--genome-build 38 \\
				--file ${params.gens_accessdir}/${track_roh.getName()} \\
				--name \\\"LOH\\\"" >> ${id}.gens_v4_const

			# Only load UPD track for proband with family
			if [[ "${params.mode}" == "family" ]]; then
				echo "gens load sample-annotation \\
					--sample-id $id \\
					--case-id $group \\
					--genome-build 38 \\
					--file ${params.gens_accessdir}/${track_upd.getName()} \\
					--name \\\"UPS\\\"" >> ${id}.gens_v4_const
			fi
		fi
		"""
	
	stub:
		"""
		touch "${id}.gens_v4_const"
		"""
}

// SV-calling //

// GATK panel+wgs //

process gatk_coverage {
	cpus 2
	memory '50GB'
	time '2h'
	container  "${params.container_gatk}"
	tag "$id"
	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${id}.tsv"), emit: coverage_tsv
		path "*versions.yml", emit: versions


	when:
		params.gatkcnv

	script:
		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		gatk --java-options "-Xmx20g" CollectReadCounts \\
			-L $params.gatk_intervals \\
			-R $params.genome_file \\
			-imr OVERLAPPING_ONLY \\
			-I $bam \\
			--format TSV -O ${id}.tsv

		${gatk_coverage_version(task)}
		"""

	stub:
		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		touch "${id}.tsv"

		${gatk_coverage_version(task)}
		"""
}
def gatk_coverage_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
	END_VERSIONS
	"""
}

process gatk_call_ploidy {
	cpus 10
	memory '50GB'
	time '2h'
	container  "${params.container_gatk}"
	tag "$id"

	input:
		tuple val(group), val(id), path(coverage_tsv)

	output:
		tuple val(group), val(id), path("ploidy.tar"), emit: call_ploidy
		path "*versions.yml", emit: versions

	script:
		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		gatk --java-options "-Xmx20g" DetermineGermlineContigPloidy \\
			--model ${params.ploidymodel} \\
			-I ${coverage_tsv} \\
			-O ploidy/ \\
			--output-prefix ${group}
		tar -cvf ploidy.tar ploidy/

		${gatk_call_ploidy_version(task)}
		"""

	stub:
		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		touch "ploidy.tar"

		${gatk_call_ploidy_version(task)}
		"""
}
def gatk_call_ploidy_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
	END_VERSIONS
	"""
}

process gatk_call_cnv {
	cpus 8
	memory '50GB'
	time '3h'
	container  "${params.container_gatk}"
	tag "$id"

	input:
		tuple val(group), val(id), path(tsv), path(ploidy), val(i), val(refpart) //TODO: is reffart a path or val


	output:
	//TODO: wtf is i
		tuple val(group), val(id), val(i), path("${group}_${i}.tar"), emit: gatk_calls
		path "*versions.yml", emit: versions

	script:
		"""
		export THEANO_FLAGS="base_compiledir=."
		set +u
		source activate gatk
		export HOME=/local/scratch
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		tar -xvf ploidy.tar
		mkdir ${group}_${i}
		gatk --java-options "-Xmx25g" GermlineCNVCaller \\
			--run-mode CASE \\
			-I $tsv \\
			--contig-ploidy-calls ploidy/${group}-calls/ \\
			--model ${refpart} \\
			--output ${group}_${i}/ \\
			--output-prefix ${group}_${i}
		tar -cvf ${group}_${i}.tar ${group}_${i}/

		${gatk_call_cnv_version(task)}
		"""

	stub:
		"""
		export THEANO_FLAGS="base_compiledir=."
		set +u
		source activate gatk
		export HOME=/local/scratch
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		source activate gatk
		touch "${group}_${i}.tar"

		${gatk_call_cnv_version(task)}
		"""
}
def gatk_call_cnv_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
	END_VERSIONS
	"""
}

process postprocessgatk {
	cpus 5
	memory '50GB'
	time '3h'
	container  "${params.container_gatk}"
	publishDir "${params.results_output_dir}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'
	tag "$id"

	input:

		tuple val(group), val(id), val(i), path(tar), path(ploidy), val(shard_no), val(shard)

	output:
		tuple val(group), val(id), path("genotyped-intervals-${group}-vs-cohort30.vcf.gz"), path("genotyped-segments-${group}-vs-cohort30.vcf.gz"), path("denoised-${group}-vs-cohort30.vcf.gz"), emit: called_gatk
		path "*versions.yml", emit: versions


	script:
		def modelshards = shard.join(' --model-shard-path ') // join each reference shard
		def caseshards = []
		// TODO: put into func
		// // join each shard(n) that's been called
	    i.each { shard_name ->
        	def shard_path = group + '_' + shard_name + '/' + group + '_' + shard_name + '-calls'
        	caseshards << shard_path
    	}
		caseshards = caseshards.join( ' --calls-shard-path ')
		version_str = postprocessgatk_version(task)
		"""
		export THEANO_FLAGS="base_compiledir=."

		for model in ${tar}; do
			tar -xvf \$model
		done

		tar -xvf ${ploidy}

		set +u
		source activate gatk
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}

		gatk --java-options "-Xmx25g" PostprocessGermlineCNVCalls \\
			--allosomal-contig X --allosomal-contig Y \\
			--contig-ploidy-calls ploidy/${group}-calls/ \\
			--sample-index 0 \\
			--output-genotyped-intervals genotyped-intervals-${group}-vs-cohort30.vcf.gz \\
			--output-genotyped-segments genotyped-segments-${group}-vs-cohort30.vcf.gz \\
			--output-denoised-copy-ratios denoised-${group}-vs-cohort30.vcf.gz \\
			--sequence-dictionary ${params.GENOMEDICT} \\
			--calls-shard-path ${caseshards} \\
			--model-shard-path ${modelshards}

		echo "${version_str}" > "${task.process}_versions.yml"
		"""

	stub:
		def modelshards = shard.join(' --model-shard-path ') // join each reference shard
		def caseshards = []
		// TODO: lsp complains about indexing var
		// // join each shard(n) that's been called
	    i.each { shard_name ->
        	def shard_path = group + '_' + shard_name + '/' + group + '_' + shard_name + '-calls'
        	caseshards << shard_path
    	}
		caseshards = caseshards.join( ' --calls-shard-path ')
		version_str = postprocessgatk_version(task)
		"""
		export THEANO_FLAGS="base_compiledir=."
		set +u
		source activate gatk
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		source activate gatk
		touch "genotyped-intervals-${group}-vs-cohort30.vcf.gz"
		touch "genotyped-segments-${group}-vs-cohort30.vcf.gz"
		touch "denoised-${group}-vs-cohort30.vcf.gz"

		echo "${modelshards}"
		echo "${caseshards}"

		echo "${version_str}" > "${task.process}_versions.yml"
		"""
}
def postprocessgatk_version(task) {
	// TODO: Reconcile this version stub with others.
	"""${task.process}:
	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')"""
}


process filter_merge_gatk {
	cpus 2
	tag "$group"
	time '2h'
	memory '1 GB'
	publishDir "${params.results_output_dir}/sv_vcf", mode: 'copy', overwrite: 'true'

	input:
		tuple val(group), val(id), path(gentotyped_intervals), path(genotyped_segments), path(denoised_copy_ration)

	output:
		tuple val(group), val(id), path("${id}.gatk.filtered.merged.vcf"), emit: merged_filtered_vcf

	script:
		"""
		filter_gatk.pl $genotyped_segments > ${id}.gatk.filtered.vcf
		mergeGATK.pl ${id}.gatk.filtered.vcf > ${id}.gatk.filtered.merged.vcf
		"""

	stub:
		"""
		touch "${id}.gatk.filtered.merged.vcf"
		"""
}


process manta {
	cpus  56
	publishDir "${params.results_output_dir}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'
	tag "$id"
	time '15h'
	memory '150 GB'
	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${id}.manta.vcf.gz"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		bams = bam.join('--bam ')

		"""
		configManta.py --bam $bams --reference ${params.genome_file} --runDir .
		python runWorkflow.py -m local -j ${task.cpus}
		mv results/variants/diploidSV.vcf.gz ${id}.manta.vcf.gz
		mv results/variants/diploidSV.vcf.gz.tbi ${id}.manta.vcf.gz.tbi

		${manta_version(task)}
		"""

	stub:
		"""
		touch "${id}.manta.vcf.gz"
		${manta_version(task)}
		"""
}
def manta_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    manta: \$( configManta.py --version )
	END_VERSIONS
	"""
}

process manta_panel {
	cpus  20
	publishDir "${params.results_output_dir}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'
	tag "$id"
	time '1h'
	memory '50 GB'


	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${id}.manta.vcf.gz"), emit: vcf
		path "*versions.yml", emit: versions

	when:
		params.sv && params.antype == "panel"

	script:
		"""
		configManta.py --bam $bam --reference ${params.genome_file} --runDir . --exome --callRegions $params.bedgz --generateEvidenceBam
		python runWorkflow.py -m local -j ${task.cpus}
		mv results/variants/diploidSV.vcf.gz ${id}.manta.vcf.gz
		mv results/variants/diploidSV.vcf.gz.tbi ${id}.manta.vcf.gz.tbi

		${manta_panel_version(task)}
		"""

	stub:
		"""
		touch "${id}.manta.vcf.gz"
		${manta_panel_version(task)}
		"""
}
def manta_panel_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    manta: \$( configManta.py --version )
	END_VERSIONS
	"""
}


process cnvkit_panel {
	cpus  5
	container  "${params.container_twist_myeloid}"
	publishDir "${params.results_output_dir}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	publishDir "${params.results_output_dir}/plots/", mode: 'copy', overwrite: 'true', pattern: '*.png'
	tag "$id"
	time '1h'
	memory '20 GB'
	input:
		tuple val(group), val(id), path(bam), path(bai)
		tuple val(group2), val(id2), path(intersected_vcf)
		tuple val(group3), val(id3), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV)

	output:
		tuple val(group), val(id), path("${id}.cnvkit_filtered.vcf"), emit: cnvkit_calls
		tuple val(group), val(id), path("${id}.call.cns"), emit: unfiltered_cns
		tuple val(group), val(id), path("${group}.genomic_overview.png"), emit: genomic_overview_plot
		tuple val(group), path("${group}_oplot.INFO"), emit: cnvkit_INFO
		path "*versions.yml", emit: versions


	script:
		"""
		cnvkit.py batch $bam -r $params.cnvkit_reference -p 5 -d results/
		cnvkit.py call results/*.cns -v $intersected_vcf -o ${id}.call.cns
		filter_cnvkit.pl ${id}.call.cns $MEAN_DEPTH > ${id}.filtered
		cnvkit.py export vcf ${id}.filtered -i "$id" > ${id}.cnvkit_filtered.vcf
		cnvkit.py scatter -s results/*dedup.cn{s,r} -o ${group}.genomic_overview.png -v $intersected_vcf -i $id
		echo "IMG overviewplot	${params.accessdir}/plots/${group}.genomic_overview.png" > ${group}_oplot.INFO

		${cnvkit_panel_version(task)}
		"""

	stub:
		"""
		touch "${id}.cnvkit_filtered.vcf"
		touch "${id}.call.cns"
		touch "${group}.genomic_overview.png"
		touch "${group}_oplot.INFO"

		${cnvkit_panel_version(task)}
		"""
}
def cnvkit_panel_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
	END_VERSIONS
	"""
}

process svdb_merge_panel {
	container  "${params.container_svdb}"
	cpus 2
	cache 'deep'
	tag "$group"
	publishDir "${params.results_output_dir}/sv_vcf/merged/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	time '1h'
	memory '1 GB'
	input:
		tuple val(group), val(id), path(vcfs)

	output:
		tuple val(group), val(id), path("${group}.merged.vcf"), emit: merged_vcf
		tuple val(group), path("${group}.merged.vcf"), emit: loqusdb_vcf
		path "*versions.yml", emit: versions


	script:
		if (vcfs.size() > 1) {
			// for each sv-caller add idx, find vcf and find priority, add in priority order! //
			// index of vcfs added
			manta_idx = vcfs.findIndexOf{ it =~ 'manta' }
			cnvkit_idx = vcfs.findIndexOf{ it =~ 'cnvkit' }
			gatk_idx = vcfs.findIndexOf{ it =~ 'gatk' }

			// find vcfs //
			manta = manta_idx >= 0 ? vcfs[manta_idx].collect {it + ':manta ' } : null
			cnvkit = cnvkit_idx >= 0 ? vcfs[cnvkit_idx].collect {it + ':cnvkit ' } : null
			gatk = gatk_idx >= 0 ? vcfs[gatk_idx].collect {it + ':gatk ' } : null
			tmp = manta + gatk + cnvkit
			tmp = tmp - null
			vcfs_svdb = tmp.join(' ')

			// find priorities //
			mantap = manta_idx >= 0 ? 'manta' : null
			gatkp = gatk_idx >= 0 ? 'gatk' : null
			cnvkitp = cnvkit_idx >= 0 ? 'cnvkit' : null
			tmpp = [mantap, gatkp, cnvkitp]
			tmpp = tmpp - null
			priority = tmpp.join(',')

			"""
			svdb \\
				--merge \\
				--vcf ${vcfs_svdb} \\
				--no_intra \\
				--pass_only \\
				--bnd_distance 2500 \\
				--overlap 0.7 \\
				--priority ${priority} \\
				--ins_distance 0 > ${group}.merged.tmp


			# copy callers out of INFO.tuple to INFO.SCOUT_CUSTOM
			add_callers_to_scout_custom.py \\
				--callers $priority \\
				--merged_vcf ${group}.merged.tmp > ${group}.merged.callers.tmp

			add_vcf_header_info_records.py \\
				--vcf ${group}.merged.callers.tmp \\
				--info SCOUT_CUSTOM . String "Custom annotations for scout" '' '' \\
				--output ${group}.merged.vcf

			${svdb_merge_panel_version(task)}
			"""
		}
		else {
			"""
			mv $vcfs ${group}.merged.vcf
			${svdb_merge_panel_version(task)}
			"""
		}

	stub:
		"""
		touch "${group}.merged.vcf"
		${svdb_merge_panel_version(task)}
		"""
}
def svdb_merge_panel_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
	END_VERSIONS
	"""
}

process postprocess_merged_panel_sv_vcf {
	cpus 2
	tag "$group"
	publishDir "${params.results_output_dir}/sv_vcf/merged/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	time '1h'
	memory '1 GB'


	input:
		tuple val(group), val(id), path(merged_vcf)
		tuple val(group2), val(id2), path(melt_vcf)

	output:
		tuple val(group), val(id), path("${group}.merged.bndless.genotypefix.melt.vcf"), emit: merged_postprocessed_vcf
		path "*versions.yml", emit: versions


	script:
		"""
		# Remove BNDs
		grep -v "BND" $merged_vcf > ${group}.merged.bndless.vcf

		# Any 0/0 GT -> 0/1, otherwise loqus will reject them.
		modify_cnv_genotypes_for_loqusdb.pl --merged_panel_sv_vcf ${group}.merged.bndless.vcf > ${group}.merged.bndless.genotypefix.vcf

		# Add MELT data to info vars:
		add_vcf_header_info_records.py \\
			--vcf ${group}.merged.bndless.genotypefix.vcf \\
			--info MELT_RANK . String "Evidence level 1-5, 5 - highest" '' '' \\
			--info MELT_QC . String "Quality of call" '' '' \\
			--output ${group}.merged.bndless.genotypefix.headers.vcf

		# Combine with MELT:
		vcf-concat  ${group}.merged.bndless.genotypefix.headers.vcf $melt_vcf | vcf-sort -c > ${group}.merged.bndless.genotypefix.melt.vcf
		${postprocess_merged_panel_sv_version(task)}
		"""

	stub:
		"""
		touch ${group}.merged.bndless.genotypefix.melt.vcf
		${postprocess_merged_panel_sv_version(task)}
		"""

}
def postprocess_merged_panel_sv_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcftools: \$(vcftools --version | cut -f 2 -d " " | tr -d "()")
	END_VERSIONS
	"""
}

process tiddit {
	cpus  2
	publishDir "${params.results_output_dir}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	time '10h'
	tag "$id"
	memory '15 GB'

	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${id}.tiddit.filtered.vcf"), emit: vcf
		path "*versions.yml", emit: versions


	when:
		params.sv && params.antype == "wgs"

	script:
		"""
		TIDDIT.py --sv -o ${id}.tiddit --bam $bam
		grep -E \"#|PASS\" ${id}.tiddit.vcf > ${id}.tiddit.filtered.vcf
		${tiddit_version(task)}
		"""

	stub:
		"""
		touch "${id}.tiddit.filtered.vcf"
		${tiddit_version(task)}
		"""
}
def tiddit_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    tiddit: \$(echo \$(TIDDIT.py 2>&1) | sed 's/^.*TIDDIT-//; s/ .*\$//')
	END_VERSIONS
	"""
}

process svdb_merge {
	cpus 2
	container  "${params.container_svdb}"
	tag "$group"
	publishDir "${params.results_output_dir}/sv_vcf/merged/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	time '2h'
	memory '1 GB'

	// TODO: make interface to this more general?
	//       e.g. tuple group, id, caller_name, calls
	input:
		tuple val(group), val(id), path(mantaV)
		tuple val(group2), val(id2), path(tidditV)
		tuple val(group3), val(id3), path(gatkV)

	output:
		tuple val(group), val(id), path("${group}.merged.bndless.vcf"), emit: merged_bndless_vcf
		tuple val(group), path("${group}.merged.vcf"), emit: merged_vcf
		path "*versions.yml", emit: versions

	script:
		if (params.mode == "family") {
			def vcfs = []
			def manta = []
			def tiddit = []
			def gatk = []

			/*
			 Order in which VCFs are merged matters when the merged SV
			 is annotated with final position/length, which affects
			 artefact matching in loqusdb.

			 A possibly better way to sort here would be to sort the
			 file by familial-relation (e.g. always sort proband-mother-father)
			 this would ensure the same merge-order regardless of sample-id
			 */

			mantaV = mantaV.collect { it.toString() }.sort()
			gatkV = gatkV.collect { it.toString() }.sort()
			tidditV = tidditV.collect { it.toString() }.sort()


			def vcf_idx = 1

			mantaV.each { _manta_vcf ->
				def tmp = mantaV[vcf_idx - 1] + ":manta" + "${vcf_idx}"
				def tmp1 = tidditV[vcf_idx - 1] + ":tiddit" + "${vcf_idx}"
				def tmp2 = gatkV[vcf_idx - 1] + ":gatk" + "${vcf_idx}"
				vcfs = vcfs + tmp + tmp1 + tmp2
				def mt = "manta" + "${vcf_idx}"
				def tt = "tiddit" + "${vcf_idx}"
				def ct = "gatk" + "${vcf_idx}"
				manta << mt
				tiddit << tt
				gatk << ct
				vcf_idx = vcf_idx + 1
			}

			prio = manta + tiddit + gatk
			prio = prio.join(',')
			vcfs = vcfs.join(' ')
			"""
			svdb \\
				--merge \\
				--vcf $vcfs \\
				--no_intra \\
				--pass_only \\
				--bnd_distance 2500 \\
				--overlap 0.7 \\
				--priority $prio \\
				--ins_distance 0 > ${group}.merged.vcf

			grep -v BND ${group}.merged.vcf > ${group}.merged.bndless.vcf

			${svdb_merge_version(task)}
			"""
		}

		else {
			tmp = mantaV.collect {it + ':manta ' } + tidditV.collect {it + ':tiddit ' } + gatkV.collect {it + ':gatk ' }
			vcfs = tmp.join(' ')
			"""
			svdb \\
				--merge \\
				--vcf $vcfs \\
				--no_intra \\
				--pass_only \\
				--bnd_distance 2500 \\
				--overlap 0.7 \\
				--priority manta,tiddit,gatk \\
				--ins_distance 0 > ${group}.merged.vcf

			grep -v BND ${group}.merged.vcf > ${group}.merged.bndless.vcf

			${svdb_merge_version(task)}
			"""
		}

	stub:
		"""
		touch "${group}.merged.vcf"
		touch "${group}.merged.bndless.vcf"

		${svdb_merge_version(task)}
		"""
}
def svdb_merge_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
	END_VERSIONS
	"""
}


process add_to_loqusdb {
	cpus 1
	publishDir "${params.crondir}/loqus", mode: 'copy' , overwrite: 'true'
	tag "$group"
	memory '100 MB'
	time '25m'
	input:
		tuple val(group), val(type), path(ped), path(vcf), path(tbi)
		tuple val(group2), path(svvcf)

	output:
		path("${group}*.loqus"), emit: loqusdb_done

	when:
		!params.noupload && !params.reanalyze

	script:

		def sv_variants_arg=""

		// Handle missing svvcf:
		if(svvcf.baseName != "NA") {
			sv_variants_arg="--sv-variants ${params.accessdir}/sv_vcf/merged/${svvcf}"
		}

		"""
		echo "-db $params.loqusdb load -f ${params.accessdir}/ped/${ped} --variant-file ${params.accessdir}/vcf/${vcf} ${sv_variants_arg}" > "${group}.loqus"
		"""

	stub:
		"""
		touch "${group}.loqus"
		"""
}

process filter_proband_null_calls {
	tag "$group"
	cpus 2
	memory '5 GB'
	time '20m'
	container  "${params.container_bcftools}"

	input:
		tuple val(group), val(id), path(sv_vcf)
		tuple val(group2), val(proband_id), val(sex), val(type)
	
	output:
		tuple val(group), val(id), path("${group}.proband.calls.vcf"), emit: filtered_vcf

	script:
		"""
		SAMPLE_INDEX=\$(bcftools query -l $sv_vcf | grep $proband_id -n | cut -f 1 -d ':' | awk '{print \$1 - 1}')
    	bcftools view -e \"GT[\$SAMPLE_INDEX]='./.'\" -o ${group}.proband.calls.vcf -O v $sv_vcf
		"""

	stub:
		"""
		touch ${group}.proband.calls.vcf
		"""

}

process tdup_to_dup {
	// in order for vep_sv to match on SVTYPE for gnomadSV TDUP needs to go
	tag "$group"
	cpus 2
	memory '5 GB'
	time '20m'

	input:
		tuple val(group), val(id), path(sv_vcf)

	output:
		tuple val(group), val(id), path("${group}.tdups_to_dups.vcf"), emit: renamed_vcf

	script:
		"""
		sed -r 's/SVTYPE=TDUP/SVTYPE=DUP/g' $sv_vcf > ${group}.tdups_to_dups.vcf
		"""

	stub:
		"""
		touch ${group}.tdups_to_dups.vcf
		"""
}

process annotsv {

	container  "${params.container_annotsv}"
	cpus 2
	tag "$group"
	publishDir "${params.results_output_dir}/annotsv/", mode: 'copy', overwrite: 'true', pattern: '*.tsv'
	time '5h'
	memory '20 GB'

	input:
		tuple val(group), val(id), path(sv_vcf)

	output:
		tuple val(group), path("${group}_annotsv.tsv"), emit: annotsv_tsv
		path "*versions.yml", emit: versions

	script:
		version_str = annotsv_version(task)
		"""
		AnnotSV -SvinputFile ${sv_vcf} \\
			-annotationMode full \\
			-annotationsDir $params.ANNOTSV_DB \\
			-outputDir ${group} \\
			-includeCI 0 \\
			-genomeBuild GRCh38
		if [ -f ${group}/*.annotated.tsv ]; then
			mv ${group}/*.annotated.tsv ${group}_annotsv.tsv
		else
		    echo "1\n" > ${group}_annotsv.tsv
		fi
		echo "${version_str}" > "${task.process}_versions.yml"
		"""

	stub:
		version_str = annotsv_version(task)
		"""
		export ANNOTSV="/AnnotSV"
		touch "${group}_annotsv.tsv"

		echo "${version_str}" > "${task.process}_versions.yml"
		"""
}
def annotsv_version(task) {
	"""${task.process}:
    annotsv: \$( echo \$(AnnotSV --version) | sed -e "s/AnnotSV //g ; s/Copyright.*//" )"""
}

process vep_sv {
	cpus 10
	container  "${params.container_vep}"
	tag "$group"
	memory '50 GB'
	time '1h'

	input:
		tuple val(group), val(id), path(vcf)

	output:
		tuple val(group), val(id), path("${group}.vep.vcf"), emit: vep_sv_vcf
		path "*versions.yml", emit: versions

	script:
		"""
		vep \\
			-i ${vcf} \\
			-o ${group}.vep.vcf \\
			--offline \\
			--merged \\
			--sift b --polyphen b --ccds --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --af --af_1kg --max_af --pubmed --uniprot --mane --tsl --appris --variant_class --gene_phenotype --mirna \\
			--synonyms $params.VEP_SYNONYMS \\
			--vcf \\
			--no_stats \\
			--fork ${task.cpus} \\
			--force_overwrite \\
			--plugin LoFtool \\
			--fasta $params.VEP_FASTA \\
			--dir_cache $params.VEP_CACHE \\
			--dir_plugins $params.VEP_PLUGINS \\
			--max_sv_size $params.VEP_MAX_SV_SIZE \\
			--distance $params.VEP_TRANSCRIPT_DISTANCE \\
			--custom file=$params.GNOMAD_SV,short_name=gnomad,fields=AF%FREQ_HOMALT%SVTYPE,format=vcf,reciprocal=1,overlap_cutoff=70,same_type=1 \\
			-cache \\
			--format vcf

		${vep_sv_version(task)}
		"""

	stub:
		"""
		touch "${group}.vep.vcf"
		${vep_sv_version(task)}
		"""
}
def vep_sv_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
	END_VERSIONS
	"""
}

process postprocess_vep_sv {
	cpus  2
	memory '10GB'
	time '1h'
	tag "$group"
	container  "${params.container_svdb}"

	input:
		tuple val(group), val(id), path(vcf)

	output:
		tuple val(group), path("${group}.vep.clean.merge.vcf"), emit: merged_processed_vcf
		path "*versions.yml", emit: versions

	script:
		"""
		# Filter variants with FILTER != . or PASS and variants missing CSQ field. SVDB creates invalid vcf header if input VCF file name contains mixtures of . _ -
		# Output file name needs to be generic
		postprocess_vep_vcf.py $vcf > postprocessed.vcf
		svdb --merge --overlap 0.9 --notag --vcf postprocessed.vcf --ins_distance 0 > ${group}.vep.clean.merge.tmp.vcf

	# --notag above will remove set
		add_vcf_header_info_records.py \\
			--vcf ${group}.vep.clean.merge.tmp.vcf \\
			--info set 1 String "Source VCF for the merged record in SVDB" '' '' \\
			--info VARID 1 String "The variant ID of merged samples" '' '' \\
			--output ${group}.vep.clean.merge.headers.tmp.vcf

		# Prepare annotations for scout:
		normalize_caller_names_in_svdb_fields.py ${group}.vep.clean.merge.headers.tmp.vcf --callers manta gatk tiddit > ${group}.vep.clean.merge.vcf
		${postprocess_vep_sv_version(task)}
		"""
	stub:
		"""
		touch "${group}.vep.clean.merge.vcf"
		${postprocess_vep_sv_version(task)}
		"""
}
def postprocess_vep_sv_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
	END_VERSIONS
	"""
}

// Query artefact db
process artefact {
	cpus 2
	tag "$group"
	time '10h'
	memory '10 GB'
	container  "${params.container_svdb}"


	input:
		tuple val(group), path(sv)
	output:
		tuple val(group), path("${group}.artefact.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		// use loqusdb dump not svdb database //
		if (params.gatkcnv) {
			"""
			svdb \\
			--query --bnd_distance 25000 --overlap 0.7 --in_occ Obs --out_occ ACOUNT --in_frq Frq --out_frq AFRQ  \\
			--db $params.svdb \\
			--ins_distance 0 \\
			--query_vcf $sv > ${group}.artefact.vcf

			${artefact_version(task)}
			"""
		}
		// for oncov1-0 still use svdb database remove in future//
		else {
			"""
			svdb \\
			--sqdb $params.svdb \\
			--query \\
			--query_vcf $sv \\
			--out_occ ACOUNT \\
			--ins_distance 0 \\
			--out_frq AFRQ > ${group}.artefact.vcf

			${artefact_version(task)}
			"""
		}

	stub:
		"""
		touch "${group}.artefact.vcf"
		${artefact_version(task)}
		"""
}
def artefact_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
	END_VERSIONS
	"""
}

process bcftools_annotate_dbvar {
	// this could be made a generic process for annotation of SVs
	tag "$group"
	cpus 2
	memory '5 GB'
	time '20m'
	container  "${params.container_bcftools}"

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.sv.dbvar.annotated.vcf.gz"), path("${group}.sv.dbvar.annotated.vcf.gz.tbi"), emit: vcf


	script:
		"""
		bcftools annotate \\
    	   -a $params.DBVAR_DEL -h $params.DBVAR_HEADERS \\
    	   -c CHROM,FROM,TO,dbvar,dbVar_status,CLNSIG,CLNASS,CLNACC \\
    	   -O v --min-overlap 0.7  $vcf | bcftools view -i 'INFO/SVTYPE="DEL"' -O z -o deletions_annotated.vcf.gz
    	   
		bcftools index --tbi deletions_annotated.vcf.gz
		bcftools annotate \\
    	   -a $params.DBVAR_DUP -h $params.DBVAR_HEADERS \\
    	   -c CHROM,FROM,TO,dbvar,dbVar_status,CLNSIG,CLNASS,CLNACC \\
    	   -O v --min-overlap 0.7 $vcf | bcftools view -i 'INFO/SVTYPE="TDUP" || INFO/SVTYPE="DUP"' -O z -o duplications_annotated.vcf.gz 
    	   
		bcftools index --tbi duplications_annotated.vcf.gz
		bcftools annotate \\
    	   -a $params.DBVAR_INS -h $params.DBVAR_HEADERS \\
    	   -c CHROM,FROM,TO,dbvar,dbVar_status,CLNSIG,CLNASS,CLNACC \\
    	   -O v --min-overlap 0.7 $vcf | bcftools view -i 'INFO/SVTYPE="INS"' -O z -o insertions_annotated.vcf.gz
    	   
		bcftools index --tbi insertions_annotated.vcf.gz
		bcftools view -i 'INFO/SVTYPE!="TDUP" && INFO/SVTYPE!="DUP" && INFO/SVTYPE!="DEL" && INFO/SVTYPE!="INS"' $vcf -O z -o others.vcf.gz
		bcftools index --tbi others.vcf.gz

		bcftools concat deletions_annotated.vcf.gz duplications_annotated.vcf.gz insertions_annotated.vcf.gz others.vcf.gz -O z -o ${group}.sv.dbvar.annotated.vcf.gz -a
		bcftools index --tbi ${group}.sv.dbvar.annotated.vcf.gz
		${bcftools_annotate_dbvar_version(task)}
		"""
	stub:
		"""
		touch ${group}.sv.dbvar.annotated.vcf.gz ${group}.sv.dbvar.annotated.vcf.gz.tbi
		${bcftools_annotate_dbvar_version(task)}
		"""
}
def bcftools_annotate_dbvar_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS
	"""
}

process add_annotsv_to_svvcf {
	cpus 2
	container "${params.container_pysam_cmdvcf}"
	memory "5 GB"
	time "20m"

	input:
		tuple val(group), path(vcf), path(tbi), path(tsv)

	output:
		tuple val(group), path("${group}.sv.annotatedSV.vcf"), emit: vcf
	
	script:
		"""
		add_annotsv.py -i $vcf -t $tsv:ACMG_class -o ${group}.sv.annotatedSV.vcf
		${add_annotsv_to_svvcf_version(task)}
		"""
	stub:
		"""
		touch ${group}.sv.annotatedSV.vcf
		${add_annotsv_to_svvcf_version(task)}
		"""
}
def add_annotsv_to_svvcf_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    python: \$(python --version 2>&1 | sed -e 's/Python //g')
		cmdvcf: 0.1
	END_VERSIONS
	"""
}



process add_omim_morbid_to_svvcf {
	cpus 2
	container "${params.container_pysam_cmdvcf}"
	memory "5 GB"
	time "20m"

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.sv.omim_morbid.vcf"), emit: vcf
	
	script:
		"""
		add_omim_morbid_sv.py -i $vcf -m $params.OMIM_MORBID_GENES -o ${group}.sv.omim_morbid.vcf
		${add_omim_morbid_to_svvcf_version(task)}
		"""
	stub:
		"""
		touch ${group}.sv.omim_morbid.vcf
		${add_omim_morbid_to_svvcf_version(task)}
		"""
}
def add_omim_morbid_to_svvcf_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    python: \$(python --version 2>&1 | sed -e 's/Python //g')
		cmdvcf: 0.1
	END_VERSIONS
	"""
}

process add_callerpenalties_to_svvcf {
	cpus 2
	container "${params.container_pysam_cmdvcf}"
	memory "5 GB"
	time "20m"

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.sv.penalty.vcf"), emit: vcf

	script:
		"""
		sv_varcall_penalties.py -i $vcf -o ${group}.sv.penalty.vcf
		${add_callerpenalties_to_svvcf_version(task)}
		"""

	stub:
		"""
		touch ${group}.sv.penalty.vcf
		${add_callerpenalties_to_svvcf_version(task)}
		"""
}
def add_callerpenalties_to_svvcf_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    python: \$(python --version 2>&1 | sed -e 's/Python //g')
		cmdvcf: 0.1
	END_VERSIONS
	"""
}


process add_geneticmodels_to_svvcf {
	cpus 2
	container "${params.container_pysam_cmdvcf}"
	memory "5 GB"
	time "20m"

	input:
		tuple val(group), val(type), path(ped), path(vcf)

	output:
		tuple val(group), val(type), path("${group}.annotatedSV.vcf"), emit: annotated_sv_vcf

	script:
		"""
		add_geneticmodels_to_svvcf.py -i $vcf -p $ped -o ${group}.annotatedSV.vcf
		${add_geneticmodels_to_svvcf_version(task)}
		"""
	stub:
		"""
		touch ${group}.annotatedSV.vcf
		${add_geneticmodels_to_svvcf_version(task)}
		"""
}
def add_geneticmodels_to_svvcf_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    python: \$(python --version 2>&1 | sed -e 's/Python //g')
		cmdvcf: 0.1
	END_VERSIONS
	"""
}

process score_sv {
	tag "$group $params.mode"
	cpus 2
	memory '10 GB'
	time '2h'
	container  "${params.container_genmod}"

	input:
		tuple val(group), val(type), path(in_vcf)

	output:
		tuple val(group), val(type), path("*.sv.scored.vcf"), emit: scored_vcf
		path "*versions.yml", emit: versions

	script:
		def model = (params.mode == "family" && params.antype == "wgs") ? params.svrank_model : params.svrank_model_s
		def group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group
		"""
		genmod score --family_id ${group_score} --score_config ${model} --rank_results --outfile "${group_score}.sv.scored.vcf" ${in_vcf}

		${score_sv_version(task)}
		"""

	stub:
		group_score = group
		"""
		touch "${group_score}.sv.scored.vcf"

		${score_sv_version(task)}
		"""
}
def score_sv_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    genmod: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
	END_VERSIONS
	"""
}

process bgzip_scored_genmod {
	tag "$group"
	cpus 4
	memory '1 GB'
	time '5m'
	publishDir "${params.results_output_dir}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'
	publishDir "${params.results_output_dir}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz.tbi'
	container  "${params.container_bcftools}"

	input:
		tuple val(group), val(type), path(scored_sv_vcf)

	output:
		tuple val(group), val(type), path("${group_score}.sv.scored.sorted.vcf.gz"), path("${group_score}.sv.scored.sorted.vcf.gz.tbi"), emit: sv_rescore
		tuple val(group), path("${group_score}.sv.scored.sorted.vcf.gz"), emit: sv_rescore_vcf
		tuple val(group), path("${group}_sv.INFO"), emit: sv_INFO
		path "*versions.yml", emit: versions

	script:
		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group
		"""
			bcftools sort -O v -o ${group_score}.sv.scored.sorted.vcf ${scored_sv_vcf}
			bgzip -@ ${task.cpus} ${group_score}.sv.scored.sorted.vcf -f
			tabix ${group_score}.sv.scored.sorted.vcf.gz -f
			echo "SV\t$type\t${params.accessdir}/vcf/${group_score}.sv.scored.sorted.vcf.gz" > ${group}_sv.INFO

			${bgzip_score_sv_version(task)}
		"""
	stub:
		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group
		"""
			touch "${group_score}.sv.scored.sorted.vcf.gz"
			touch "${group_score}.sv.scored.sorted.vcf.gz.tbi"
			touch "${group}_sv.INFO"

			${bgzip_score_sv_version(task)}
		"""
}
def bgzip_score_sv_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS
	"""
}

process compound_finder {
	cpus 2
	tag "$group ${params.mode}"
	publishDir "${params.results_output_dir}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
	memory '10 GB'
	time '2h'

	input:
		tuple val(group), val(type), path(sv_vcf), path(sv_tbi), path(ped), path(snv_vcf), path(snv_tbi)

	output:
		tuple val(group), path("${group_score}.snv.rescored.sorted.vcf.gz"), path("${group_score}.snv.rescored.sorted.vcf.gz.tbi"), emit: vcf
		tuple val(group), path("${group}_svp.INFO"), emit: svcompound_INFO
		path "*versions.yml", emit: versions


	when:
		params.mode == "family" && params.assay == "wgs"


	script:
		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group

		"""
		compound_finder.pl \\
			--sv $sv_vcf --ped $ped --snv $snv_vcf \\
			--osv ${group_score}.sv.rescored.sorted.vcf \\
			--osnv ${group_score}.snv.rescored.sorted.vcf \\
			--skipsv
		bgzip -@ ${task.cpus} ${group_score}.snv.rescored.sorted.vcf -f
		tabix ${group_score}.snv.rescored.sorted.vcf.gz -f
		echo "SVc	$type	${params.accessdir}/vcf/${group_score}.sv.scored.sorted.vcf.gz,${params.accessdir}/vcf/${group_score}.snv.rescored.sorted.vcf.gz" > ${group}_svp.INFO

		${compound_finder_version(task)}
		"""

	stub:
		group_score = group
		"""
		touch "${group_score}.snv.rescored.sorted.vcf.gz"
		touch "${group_score}.snv.rescored.sorted.vcf.gz.tbi"
		touch "${group}_svp.INFO"

		${compound_finder_version(task)}
		"""
}
def compound_finder_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	END_VERSIONS
	"""
}


process output_files {
	cpus 2
	memory '1GB'
	time '1h'

	input:
		tuple val(group), path(files, stageAs: "?/*") // Needed to stage files with identical names

	output:
		tuple val(group), path("${group}.INFO"), emit: yaml_INFO

	script:
		files = files.join( ' ' )

		"""
		cat $files > ${group}.INFO
		"""

	stub:
		"""
		touch "${group}.INFO"
		"""
}


process svvcf_to_bed {
	publishDir "${params.results_output_dir}/bed", mode: 'copy' , overwrite: 'true'
	tag "group"
	memory '1 GB'
	time '1h'
	cpus 2

	input:
		tuple val(group), path(vcf)
		tuple val(group2), val(proband_id), val(sex), val(type)

	output:
		path("${group}.sv.bed")

	when:
		params.antype != "panel"


	script:
		"""
		cnv2bed.pl --cnv ${vcf} --pb ${proband_id} > ${group}.sv.bed
		"""

	stub:
		"""
		touch "${group}.sv.bed"
		"""
}

process plot_pod {
	container  "${params.container_pod}"
	publishDir "${params.results_output_dir}/pod", mode: 'copy' , overwrite: 'true'
	tag "$group"
	time '1h'
	memory '1 GB'
	cpus 2

	input:
		tuple val(group), path(snv)
		tuple val(group2), path(cnv), val(type), path(ped)
		tuple val(group3), val(id), val(sex), val(type3)

	output:
		tuple path("${id}_POD_karyotype.pdf"), path("${id}_POD_results.html")


	script:
		"""
		parental_origin_of_duplication.pl --snv $snv --cnv $cnv --proband $id --ped $ped
		"""

	stub:
		"""
		touch "${id}_POD_karyotype.pdf"
		touch "${id}_POD_results.html"
		"""
}

process create_yaml {
	publishDir "${params.results_output_dir}/yaml", mode: 'copy' , overwrite: 'true', pattern: '*.yaml'
	publishDir "${params.results_output_dir}/yaml/alt_affect", mode: 'copy' , overwrite: 'true', pattern: '*.yaml.*a'
	publishDir "${params.crondir}/scout", mode: 'copy' , overwrite: 'true', pattern: '*.yaml'
	errorStrategy 'retry'
	maxErrors 5
	tag "$group"
	time '5m'
	memory '1 GB'

	input:
		tuple val(group), val(id), val(diagnosis), val(assay), val(type), val(clarity_sample_id), val(analysis)
		tuple val(group2), val(type2), path(ped)
		tuple val(group3), path(INFO)

	output:
		tuple val(group), path("${group}.yaml*"), emit: scout_yaml

	script:
		assay = params.dev ? "dev,${analysis}" : "${assay},${analysis}"
		"""
		create_yml.pl \\
			--g "${group},${clarity_sample_id}" \\
			--d "$diagnosis" \\
			--panelsdef "$params.panelsdef" \\
			--out "${group}.yaml" \\
			--ped "$ped" \\
			--files "$INFO" \\
			--assay "$assay" \\
			--antype "$params.antype" \\
			--extra_panels "$params.extra_panels"
		"""

	stub:
		"""
		touch "${group}.yaml"
		"""
}

process combine_versions {
	publishDir "${params.results_output_dir}/versions", mode: 'copy', overwrite: 'true', pattern: '*.versions.yml'

	input:
		val(group)
		val(versions)

	output:
		path("${group}.versions.yml")

	script:
		versions_joined = versions.sort{ my_it -> my_it.name }.join(" ")
		"""
		cat $versions_joined > ${group}.versions.yml
		"""

	stub:
		versions_joined = versions.sort{ my_it -> my_it.name }.join(" ")
		"""
		cat $versions_joined > ${group}.versions.yml
		"""
}
