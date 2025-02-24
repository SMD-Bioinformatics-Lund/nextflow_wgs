# Outputs

When running the pipeline, two separate results folders are produced. One contains the results and is located in the path specified by the `params.outdir` param. The second is the CRON-dir (`params.crondir`) and by default located in `outdir/cron`. In the setup run at CMD, Lund, this folder is transferred from the computational cluster (Hopper) to a Region Skåne server (Lennart).

## Cron-folder

| Subfolder | Mode | Process              | Name            | Description                       |
| --------- | ---- | -------------------- | --------------- | --------------------------------- |
| gens      | wgs  | `generate_gens_data` | `<group>.gens`  | Input data to Gens                |
| loqus     |      | `add_to_loqusdb`     | `<group>.loqus` | Artefact frequencies to LoqusDB   |
| qc        |      | `qc_to_cdm`          | `<group>.cdm`   | Sample quality information to CDM |
| scout     |      | `scout`              | `<group>.yaml`  | Scout input YAML file             |

## Results-folder

| Subfolder      | Function    | Mode  | Process                         | Name                                             | Description                                                                       |
| -------------- | ----------- | ----- | ------------------------------- | ------------------------------------------------ | --------------------------------------------------------------------------------- |
| annotsv        | Annotation  |       | `annotsv`                       | `<group>_annotsv.tsv`                            | Annotation from AnnotSV                                                           |
| bam            | Mapping     |       | `markdup`                       | `<group>_dedup.bam[.bai]`                        | Deduplicated BAM-files                                                            |
| bam            | Mapping     | wgs   | `fetch_MTseqs`                  | `<group>_mito.bam[.bai]`                         | Mapped mitochondrial reads                                                        |
| bed            | SV calling  |       | `svvcf_to_bed`                  | `<group>.sv.bed`                                 | SV ranges                                                                         |
| bqsr           | SNV calling |       | `bqsr`                          | `<group>.bqsr.table`                             | Base-quality score recalibrations for Sentieon DNAScope                           |
| cov            | Coverage    |       | `chanjo_sambamba`               | `<group>.bwa.chanjo.cov`                         | Coverage calculations                                                             |
| cov            | Coverage    | onco  | `depth_onco`                    | `<group>.lowcov.overlapping.bed`                 | Coverage calculations                                                             |
| cov            | Coverage    | wgs   | `gatkcov`                       | `<group>.denoisedCR.tsv`                         | Coverage calculations                                                             |
| cov            | Coverage    | wgs   | `gatkcov`                       | `<group>.standardizedCR.tsv`                     | Coverage calculations                                                             |
| cov            | Coverage    | wgs   | `d4_coverage`                   | `<group>_coverage.d4`                            | Coverage calculations                                                             |
| ped            | Pedigree    | trio  | `create_ped`                    | `<group>_[ma\|fa\|base].ped`                     | Pedigree information                                                              |
| ped            | Pedigree    | trio  | `peddy`                         | `<group>_[ma\|fa\|base].peddy.ped`               | Pedigree information                                                              |
| ped            | Pedigree    | trio  | `peddy`                         | `<group>.sex_check.csv`                          | Control of assigned sex                                                           |
| ped            | Pedigree    | trio  | `madeline`                      | `<group>_[ma\|fa\|base].ped.madeline.xml`        | Madeline pedigree data                                                            |
| plot_data      | Visualize   | wgs   | `generate_gens_data`            | Various files                                    | Files needed for visualization in Gens                                            |
| plots/SMNcnc   | Visualize   | wgs   | `SMNCopyNumberCaller`           | `smn_<group>.pdf`                                | Copy number reports for SMN1 and SMN2                                             |
| plots/reviewer | Visualize   | wgs   | `reviewer`                      | `<group>/<protein_ids>.svg`                      | Visualizations of repeat expansions                                               |
| plots/mito     | Visualize   | wgs   | `run_haplogrep`                 | `<group>.haplogrep.png`                          | Mitochondria haplogroup overview                                                  |
| plots/mito     | Visualize   | wgs   | `run_eklipse`                   | `<group>_eklipse.png`                            | Mitochondrial deletions visualized                                                |
| plots          | Visualize   | trio  | `upd_table`                     | `<group>.UPDTable.xls`                           | Uniparental disomy overview                                                       |
| overview_plot  | Visualize   | trio  | `overview_plot`                 | `<group>.genomic_overview.png`                   | High-level illustration of uniparental disomy and runs of homozygosity            |
| qc             | QC          |       | `merge_qc_json`                 | `<group>.QC`                                     | Sample QC information                                                             |
| pod            | Visualize   | trio  | `plot_pod`                      | `<group>_POD_karyotype.pdf`                      | CNV origin illustration                                                           |
| pod            | Visualize   | trio  | `plot_pod`                      | `<group>_POD_results.html`                       | CNV origin illustration                                                           |
| sv_vcf         | SV calling  | wgs   | `postprocessgatk`               | `denoised-<group>-vs-cohort30.vcf.gz`            | GATK CNVCaller output                                                             |
| sv_vcf         | SV calling  | wgs   | `postprocessgatk`               | `genotyped-intervals-<group>-vs-cohort30.vcf.gz` | GATK CNVCaller output, call result for each binned interval                       |
| sv_vcf         | SV calling  | wgs   | `postprocessgatk`               | `genotyped-segments-<group>-vs-cohort30.vcf.gz`  | GATK CNVCaller output, intervals combined in reference and variant segments       |
| sv_vcf         | SV calling  | wgs   | `filter_merge_gatk`             | `<group>.gatk.filtered.merged.vcf`               | GATK CNVCaller output, filtered non-variant intervals, combined adjacent segments |
| sv_vcf         | SV calling  | wgs   | `manta`                         | `<group>.manta.vcf.gz`                           | Manta caller output                                                               |
| sv_vcf         | SV calling  | panel | `manta_panel`                   | `<group>.manta.vcf.gz`                           | Manta caller output (adapted for panels)                                          |
| sv_vcf         | SV calling  | panel | `delly_panel`                   | `<group>.delly.vcf.gz`                           | Delly caller output                                                               |
| sv_vcf         | SV calling  | panel | `cnvkit_panel`                  | `<group>.cnvkit_filtered.vcf`                    | CNVKit caller output                                                              |
| sv_vcf         | SV calling  | panel | `svdb_merge_panel`              | `merged/<group>.merged.filtered.melt.vcf`        | Merged panel results                                                              |
| sv_vcf         | SV calling  | wgs   | `tiddit`                        | `<group>.tiddit.filtered.vcf`                    | TIDDIT caller outputs                                                             |
| sv_vcf         | SV calling  | wgs   | `svdb_merge`                    | `merged/<group>.merged.vcf`                      | Merged SV results                                                                 |
| sv_vcf         | SV calling  | wgs   | `svdb_merge`                    | `merged/<group>.merged.bndless.vcf`              | Merged SV results, only CNV calls                                                 |
| vcf            | SNV calling | wgs   | `vcfbreakmulti_expansionhunter` | `<group>.expansionhunter.vcf.gz`                 | Expansion-hunter variants                                                         |
| vcf            | SNV calling | onco  | `melt`                          | `<group>.melt.merged.vcf`                        | Mobile elements variants                                                          |
| vcf            | SNV calling | onco  | `intersect_melt`                | `<group>.melt.merged.intersected.vcf`            | Mobile elements within the ranges of interest                                     |
| vcf            | SNV calling |       | `split_normalize`               | `<group>.multibreak.vcf`                         | Combined variants                                                                 |
| vcf            | SNV calling |       | `split_normalize`               | `<group>.norm.uniq.DPAF.vcf`                     | Combined variants filtered for number of reads per allele                         |
| vcf            | SNV calling |       | `run_mutect2`                   | `<group>.mutect2.vcf`                            | Mitochondrial variants                                                            |
| vcf            | SNV calling |       | `vcf_completion`                | `<group>.scored.vcf.gz[.tbi]`                    | BGZipped and tabix indexed final VCF                                              |
| vcf            | SNV calling |       | `fastgnomad`                    | `<group>.SNPs.vcf`                               | Frequency-filtered variants to investigate inheritance patterns (UPD and ROH)     |
| vcf            | SNV calling |       | `score_sv`                      | `<group>.sv.scored.sorted.vcf.gz[.tbi]`          | Combined and scored SV variants                                                   |
| vcf            | SNV calling | wgs   | `compound_finder`               | `<group>.snv.rescored.sorted.vcf.gz`             | For trio, rescored based on inheritance patterns                                  |
| versions       | Versions    |       | `combine_versions`              | `<group>.versions.yaml`                          | Versions for all executed processes                                               |
| yaml           | Scout       |       | `create_yaml`                   | `<group>.yaml`                                   | YAML file for input into Scout                                                    |
