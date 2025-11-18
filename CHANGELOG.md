# CHANGELOG

### TBD
* Use generate_gens_data.PY for Gens v4 input

### 3.19.0
* constitutional profile
  - bed and interval files
  - SV references
* SV logic changes for panels
 - added a function to count SV calls, if 0 no annotation

### 3.18.5
* Update bed file:
  - New user-added genes
  - Clinvar updated from 2024-12-30 to 2025-08-17
  
### 3.18.4
* Reinstate publishDir for verifybam. Will only have effect on wgs.

### 3.18.3
* Update deploy script to also copy modules folder

### 3.18.2
* Fix father/mother mismatch flip in new Gens UPD tracks
* Split UPD and ROH into two separate tracks for new Gens tracks

### 3.18.1
* Changed annotation keys for dbVar. CLNSIG CLNACC now available for scout for SVs

### 3.18.0
* Add ID-SNPs

### 3.17.2
* Fix so that modycf does not crash in split_normalize
* Sort QC output JSON on keys to allow easy diffing in PipeEval
* Add GitHub action to run Pytest

### 3.17.1
* Fix median calculation bug in qc_sentieon.pl

### 3.17.0
* New process "prepare_gens_v4_input.py" with new script to generate inputs needed for Gens
* Add POD scripts to the bin/pod folder (still executed from within a container, but should be run from bin)

### 3.16.1
* Increase manta max run time

### 3.16.0
* AnnotSV bump 3.4.4
* bcftools bump 1.21
* prescore_sv.pl retirement, replace by four python scripts
  * add_annotsv.py: adds annotsv score to variants
  * add_genetic_models_to_svvcf.py: adds genetic models, and homo- hemizygous deletion annotations to vcf
  * sv_varcall_penalies.py: adds caller specific penalties to variants with only one caller
  * add_omim_morbid_sv.py: adds omim morbid status to variants overlapping genes in OMIM-AUTO gene list
* dbvar/clinvar for SVs as a separate bcftools annotate process
* gnomadSV with VEP, and bumped to latest version

### 3.15.14
* added a oncov1 profile, using the same logic as oncov2 but with relevant normal references

### 3.15.13
* remove unnecessary sentieon mito qc output

### 3.15.12
* correct output file name for myeloid const for freebayes

### 3.15.11
* Compress gvcf_combine gVCF

### 3.15.10
* theano base_dir flag was not exported correctly for postprocess GATK. Should stop random crashes from happening

### 3.15.9
Update sentieon from 202308.01 to 202308.03

### 3.15.8
Add ./workflow dirpath to deploy script

### 3.15.7
* Annotation-only runs restored for single SNV/indel VCFs
* Add entry-point for VCFs (expecting equivalent of `split_normalize.out.vcf_intersected`). 

### 3.15.6
* Add workaround so copybam actually emits copied BAMs

### 3.15.5
* Update GNOMAD genomes and exomes to v4.1 in `annotate_vep` and `indel_vep`
* Update popmax attribute to grpmax (same attribute renamed in v4: https://gnomad.broadinstitute.org/news/2023-11-genetic-ancestry/) 
  - The INFO-level keys created by `modify_vcf_scout.pl` still use "popmax" instead of grpmax.
* `myeloid_const` profile now uses same rank models as onco 
  - The desired behavior of disabled artefact scoring already existed in original onco rank models. Thus separate rank model files for myeloid were not actually needed.

### 3.15.4
* Update dbSNFP to 4.9
* Add BayesDel vep annotations to SNVs
  - Adds BayesDel addAF/noAF predictions and scores

### 3.15.3 [HOTFIX]
* Match `*q.gz*` files as fastq.gz

### 3.15.2
* Re-add GERP++ snv/indel annotations

### 3.15.1
* Update VEP to 113.3
* Update VEP fasta/cache to 113.0
* Remove `SVTYPE` bug workaround from `vep_sv`

### 3.15.0

* Rewrite main.nf to DSL2
* Refactor everything into one workflow `NEXTFLOW_WGS`.
* Remove contamination check for non-wgs profiles
* Remove `melt_qc_val` process (now exists in main workflow)
* Remove `dummy_svvcf_for_loqusdb` process
* Temp disable annotation-only runs (will be re-added later)

### 3.14.5
* Include --case-id flag with group-ID in Gens load command

### 3.14.4
* Routine update of bed intersect file


### 3.14.3
* Fix rankscore parsing in `cnv2bed.pl`

### 3.14.2
* dev-assay for create-yaml, triggered when running with `--dev` flag

### 3.14.1

* Adds basic flake8-based linting 
* Removes unused scripts from `bin/`
* Fixes wrong var name in `bin/normalize_caller_names_in_svdb_fields.py`
* Fix wrong var assignment in `bin/normalize_caller_names_in_svdb_fields.py` that led to caller names not being normalized for wgs trios.

### 3.14.0

#### Updated:

* SVDB to 2.8.2

#### Added
* Util script `add_vcf_header_info_records.py`
* Util script `normalize_caller_names_in_svdb_fields.py`
* Util script `add_callers_to_scout_custom.py`
* Util script `modify_cnv_genotypes_for_loqusdb.pl` (formerly `filter_panel_cnv.pl`)
* `main.nf`: New process `postprocess_merged_panel_sv_vcf`
* `main.nf`: New process `add_omim`

#### Changed:

* Rename `cleanVCF.py` -> `postprocess_vep_vcf.py` 
* Move most container path specs into config
* Split out `add_omim.pl` to own process and move hardcoded db path into config.
* `main.nf`: Rename `postprocess_vep` to `postprocess_vep_sv`

#### Deleted:

* Remove `merge_callsets.pl` from `svdb_merge` process
* Remove delly process and associated params / script code
* Remove onco v1 profile
* Remove `filter_panel_cnv.pl` 
* Remove `omim_genes.txt` (converted to external db)

### 3.13.1
* Add support for including specific panels from other institutes in Scout yaml.
* Scout yaml now get panels and other info added for "pipeeval" profile
* Some refactor of create_yaml.pl

### 3.12.3
* Replace RNU2-4 gene coordinates with RNU4-2 coordinates in wgs intersect bed. 

### 3.12.2
* Sort processes in versions yaml and images in Scout yaml

### 3.12.1
* Update Genmod version allowing control of penalty
* Assign non-scored components for single GIAB
* Fix `compound_finder.pl` such that it converts floats to int (i.e. 5.0 -> 5)
* Remap all start-from-BAM channels that flips the ID <-> group
* Copy in bam file to work dir when running from bam rather than accessing it directly in its original location
* Sort the order of vcfs to `gvcf_combine` for stable SNV-calls

### 3.11.2
* * Ensure that input VCFs are always supplied in the same alphanumeric order to `svdb_merge` when running trio analysis (see [#172](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/issues/172))

### 3.11.1
* Add a process to get contamination values from verifybamid2 software.
* Update configs/nextflow.hopper.config with a specific verifybamid2 container.
* Update configs/nextflow.hopper.config with specific SVDPrefix files for panel and wgs.

### 3.10.4
* Added --format vcf to `vep_sv` to fix for cases where vcf file carries no variants. 
 
### 3.10.3
* Add workaround to enable loqusdb export for runs where SV calling is disabled
* Rename myeloid_const loqusdb to `loqusdb_myeloid_const`
* Disable artefact scoring in `myeloid_const` rank models
	
### 3.10.2
* Fix mito QC stats JSON conversion for samples started from old bams with updated sample ids. 

### 3.10.1
* Update config for bed intersect
* Some fixes to the logging of the bed intersect script

### 3.9.10
* Use reduced gene_panel JSON to avoid adding dead/archived panels to new scout cases 
* Add lennart-side script/worker CRON job to generate new gene panel JSON

### 3.9.9
* Extend the update_bed.pl script to handle multiple input files
* Rewrite to Python and add tests

### 3.9.8
* Reverted removed code in gene panel matches, caused missing gene panels for onco samples

### 3.9.7
* Solved trio eklipse image being wrongly added to yaml
* Removed outdated regex matches for genepanel, would remove important gene panels
* General clean-up of create_yml.pl

### 3.9.6
* Fix bug where wrong tuple value unpacked as group and sample id in `bqsr` when starting run from bam


### 3.9.5
* Fixed faulty if-condition for annotsv, would result in empty annotsv tsv everytime

### 3.9.4
* Use -K flag in bwa-mem for consistent results

### 3.9.3

* Re-optimized profiles wgs and onco. More memory allocations
* Added flag for reanalyze for bjorn to hook into

### 3.9.2

* Add updated and more communicative deploy script
* Remove or rename other deploy scripts

### 3.9.1

* Update MODY-cf configs to use the same as onco
* Clean up in MODY-cf config post merge

### 3.9.0
	
* Give the Sentieon container path by a parameter in the config file
* Update the Sentieon container to 202308 version
* Split out the `sentieon_qc` post-processing into its own process `sentieon_qc_postprocess`
* Update the Perl script used in `sentieon_qc_postprocess` to take input parameters as explicit arguments
* Update intersect file to latest used version of ClinVar (20231230)
* Update fastp to 0.23.4 and move to own container to fix reproducibility issue ([#143](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/issues/143))
* Update CADD to v1.7
* Increase `inher_models` processing time
* Updated VEP from 103.0 to 111.0
* Updated VEP fasta from 98.0 to 111.0
* Updated VEP cache from 103.0 to 111.0
* Moved VEP parameters from processes to config
* Disabled vep `--everything` to disable VEP annotation w/ GNOMAD
* Removed deprecated `--af_esp` from `--everything`
* Tentative update of scout ranking.  
* cleanVCF.py now removes records missing CSQ-field.
* Add `SVTYPE` VEP 111 bug workaround in `vep_sv` process. (See  [Ensembl/ensembl-vep#1631](https://github.com/Ensembl/ensembl-vep/issues/1631#issuecomment-1985973568))
* Add VEP105 - 111 annotations to all rank models in use
* Fix onco model filename version (v5 rank model was misnamed as v4 in production)

### 3.8.2
* Re-enable D4 file generation (for Chanjo2)

### 3.8.1
* Disable Chanjo2

### 3.8.0
* Add mody-cf profile

### 3.7.14
* Run D4 coverage for full file

### 3.7.13
* Further simplifications of the checklist template

### 3.7.12
* Trim down size of checklist template, and add check for entering used test samples

### 3.7.11
* Add d4 file path directly to Scout YAML

### 3.7.10
* Tag Mitochondrial variants with GQ, loqusdb enabling

### 3.7.9
* Add CRON file to load Chanjo2

### 3.7.8
* added csv-file to onComplete function to accomodate CCCP

### 3.7.7
#### create_yml.pl
* small name change for myeloid constitutional to match clarity
* removed custom_images header for samples without images as pydantic would crash in scout load

### 3.7.6
* Add d4 coverage calculations to the workflow

### 3.6.6
* Fix genmod caller-penalty bug for GATK GQC vals ([#170](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/issues/170))

### 3.6.5
* Remove bgzip and gunzip from versions
* Some cleanup in version documentation and code

### 3.6.4
* Use new docs as main entry point in repo
* Start removing old docs

### 3.6.3
* Update software responsible list in docs

### 3.6.2
* Added changelog reminder to github workflows

### 3.6.1
* Adding a new variant catalogue for expansionhunter/stranger/reviewer

### 3.5.11
* Add `documentation` to change type category in PR template.

### 3.6.0
* Changed melt configs, added flags: exome, removed flags: cov (was being used improperly)
* Added priors to `mei_list`, and changed `mei_list` to a new location in config
* Changes has been verified, report can be found internally

### 3.5.10
* Changed path to normal-pool-refs for gens. Uses masked hg38 references

### 3.5.9
* Add first iteration of updated documentation

### 3.5.8
* Move out resource files from `main.nf` to `nextflow.config`
* Move the selected fields for PHYLOP and PHASTCONS in vep to be specified in the process, similarly to the other plugins/custom fields

### 3.5.7
* Clean out unused files in repo root directory

### 3.5.6
* Add Github PR template/test documentation

### 3.5.5
* Update the cron log directory to use the `params.crondir` folder as base

### 3.5.4
* Add version outputs from all processes that use external software.
* Add stubs to processes to allow performing stub runs.

### 3.5.3
- Hotfix, increase melt sensivity by increasing amount of reads melt are alowed to use in RAM. 

### 3.5.2
- MELT is no longer filtered on location based upon regex names INTRONIC/null/PROMOTER, instead added a intersect towards bedfile. This will show splice site variants

### 3.5.1
* Add REVEL (Rare Exome Variant Ensemble Learner) Scores to VEP annotations (VEP `REVEL_rankscore` and `REVEL_score`)
	
### 3.5.0

#### Added

* Two processes for computing mitochondrial seq QC data from mt bam files and saving to JSON:
* Script `bin/merge_json_files.py` to merge 1 or more JSON files into one JSON.  Used to generate the final `{id}.QC` from the json output of the processes `sentieon_qc` and `sentieon_mitochondrial_qc`.
* Script `bin/mito_tsv_to_json.py` to extract and convert mtQC data from `sentieon_mitochondrial_qc` process output to json

#### Changed

* process `sentieon_qc` outputs to intermediate `{id}_qc.json` file instead of the final `{id}.QC`

### 3.4.3
#### patch genes
- added two more genes to expansionhunter variant catalogue.

### 3.4.2
#### hotfix
- dont print Mitochondrion, we handle the mitochondrion seperatly in the pipeline, caused loqusdb errors

### 3.4.1
#### minor additions/edits
- fixed filepaths for access-dir for myeloid profile in nextflow.config
- fixed assay name for create_yml.pl so yaml-file gets correct institute owner for myeloid samples

### 3.4.0
#### new features
- added a script to update wgs-bed file with current clinvar intron + intergenic regions. Also produces a log file of what's been added and removed
- added support to dry run vcf for testing scoring

### 3.3.3
#### minor improvements
- merged cnv2bed branch, small updates to color scheme for Alamut import files for CNVs
- increase time limit of create_pedigree

### 3.3.2
#### performance improvements
- added retries to vcfanno, file-system caching bug out?
- removed deep caching from freebayes(onco only) weird bug?

### 3.3.1
#### bug fixes
- alt affect type was lost for SNV<->SV compound, would get mixed up
  - added type and joined upon the value
- compounds for only SNVs for alt affect duos was wrongly renamed, added a sed-command

### 3.3.0
#### new features
- oncov2-0 and wgs profiles now both use loqusdb dumps for SV artefact annotations
- create pedigree has completely changed, now it is a separate perl-skript
  - creates one pedigree per affections status of parent, i.e in a trio, three ped-files with mother affect/father affected/no parent affected(default loaded into scout)
  - will calculate all states per genomod score, per vcf
  - optionally load these cases into scout, located in a subcategory in yaml-output-folder
- oncov2-0 now implemented, uses the old onco profile and oncov1-0 uses oncov1-0 profile (to be discontinued)
  - No longer use delly SV-caller, instead use GATK + CNVkit + manta
  - new version of MELT that catch much more important variation
  - indicator of onco-version in rankmodel-name of yaml-file, visable in scout case page

### 3.2.6
#### enhancement
- Added regex to support wgs-hg38-XXXX. suffix to run wgs-profile with different flags. ie --noupload true, no cdm/loqusdb upload for reruns

### 3.2.5
#### Bugfix
- Fixed a serious bug in prescore_sv.pl, would randomly chose proband-id for duos


### 3.2.4
#### Feature improvement
- Added SVs to loqusdb load. Using scored snv-vcf for correct MT->M notation

### 3.2.3
#### Bugfixes
- genmod patch not taking effect in singularity, switched to a smaller genmod container with patch
- updated processes in main.nf for above container

### 3.2.2
#### Bugfixes and improvement
- Changes to custom_images in yaml, case/str
- Added support for reviewer, activate when scout is updated

### 3.2.1
#### Bugfixes and improvement
- Image sizes for mitochondrial plots in yaml
- resource management in processes
- gatk-ref moved to cached directory

### 3.2.0
### Minor release
- Added support for loading images into scout, each process generating a plot can now be added as a path to scout-yaml
- Some support for Grace (new cluster)

### 3.1.2
#### Bugfixes
- CDM load-file only created for one individual of family, fixed (join function corrected)
- increased memory allocation for onco_depth
- removed shards from dedup, caused malformed output for dedup_metrics. Works as intended still

### 3.1.1
#### Bugfixes
- params.assay for onco has historical name, depth_onco process used wrong value
- using shard-specification for depup caused faulty dedupmetrics file

### 3.1.0
#### Performance Improvements
- Removed all distribution of sentieon except first alignment step option
- Removed bam.toRealPath() from all processes. Bai files are now given along with bam files if alignment is to be skipped. More down below
#### New/Updated Features
- VCF start removed temporarily
- BAM start now work better. Add headers bam + bai to csv with corresponding files
- BATCH start now available for onco-samples. Thorough channel joining and removal of distributed sentieon made it possible (does not work for wgs profile!)

### 3.0.9
#### bug fixes
- fixed grep for multi-allelic FP loci

### 3.0.8
#### feature
- added hash element rstreshold (rankscore threshold), that if defined overwrites defualt -1 to createyml.pl

### 3.0.7
#### feature
- added panel depth as alternative to chanjo for panel-data

### 3.0.6
#### bug fixes
- correctly assigned theanoflag for gatk coverage and ploidy, would in rare cases cause crashes

### 3.0.5
#### new funcion
- GENS middleman command added to `generate_gens_data`. Needed for loading of data into GENS thorugh cron and middleman
#### bug fixees
- REViewer now loops through a perl shell script instead of bash. Low covered loci error no longer crash all other svg-image generation
- fixed a typo which named all svgs as 7156, a validation and verification sample

### 3.0.4
- recurring multi-allelic variant @MT:955 keeps vcfmultibreak in a never ending loop
- grep -v ^MT 955

### 3.0.3
- ignore errors of REViewer
- loqusdb faulty input caused wrong imports to loqusdb, now fixed

### 3.0.2
- GAV replaced with REViewer
- Stranger 0.8 with updated variant catalogue

### 3.0.1
- source activate gatk to all gatk processes that use 4.1.9 and set +eu. unbound errors
- increased memory allocation for several gatk and mito processes
- sharded merged bam and non-sharded bam now produces output for dedup too. Locuscollector no longer redirects bam. This saves upto 70% of temporary files!

### 3.0.0
#### Summary of changes
new functions
### main.nf
- added mito-calling
  - mutect2
  - hmtnote
  - haplogrep
  - eklipse
  - modifications to filter_indels (new VEP fields)
  - modifications to `modify_vcf_scout.pl`, ignore maxentscan for M
- added SMNCopyNumberCalling
- New VEP container and version (103)
- gatk cnv calling
  - adjustments to all affected scripts
### container
- new container specifically for madeline2
- main container now includes all software except madeline2 and VEP
  - new conda environments
  - updates to Expansionhunter
  - updates to Stranger
  - updated GATK version
  - updated Sentieon
  - added: haplogrep, hmtnote, eklipse, melt, graphalignmentviewer, SMNcopynumbercaller, CNVkit and imagemagick

### minor improvements
- group and sample IDs of outputs re-thought
- contig synonyms for VEP
- BAM start working better


### 2.1.12
#### params.panelsdef
- added path to symlinked latest weekly definition.

### 2.1.11
#### qc-values added
- added pf missmatch and error rates to qc-json

### 2.1.10
#### scout presentation
- rankmodels now separate VEP-consequence from AnnotSVrank and dbvar in Consequence and Clinical_significance respectively


### 2.1.9
#### bugfix
- rescore.nf had wrongly named variable in output for bamfiles 

### 2.1.8
- create\_yml.pl now recieved gene_panel content from hopper-json. no longer require scout-vm connectivity

### 2.1.7
- clincalwes now has correct loqusdb not piggybacking of onco

### 2.1.6
- timelimit increases, scratch and stage in/out for processes

### 2.1.5

- create\_yml.pl added ahus analysis for wgs_hg38 assay. Stinking mess initiated, please correct

### 2.1.4

- create_yml.pl added hemato analysis for clinicalwesv1-0 assay. Correct institute for myeloid normals

### 2.1.3

#### Features
- create_yml.pl now has a hash with all definable scout import-fields per assay, allowing easier additions and modifications to/of assays.

### 2.1.2

#### Fixes
- Fix a bug that generated corrupt Gens json files...

### 2.1.1

#### Features
- Generate a json with data for the Gens overview plot, to allow quicker loading in Gens.


### 2.1.0

#### Features
- added rescoring function through rescore.nf
#### misc
- fixed naming of expansionhunter vcfs
- new rankmodels for wgs profile. 5.1 (loqusdb cutoffs and VEP-csq scoring)


### 2.0.2

#### Features
- added specific delly filtering script
  - now correctly filters breakpoints outside panel
- `filter_panel_cnv.pl` now only annotates for scout
  - delly precise/imprecise annotation
- new artifact database for SV-calling for WGS and oncogenetics

### 2.0.1

#### Features
- added container and git-hash to logging

#### Fixes
- pathing through freebayes cause nextflow to not recieve completion status, now wgs is run through freebayes with touch-command only

### 2.0

#### Features
- Optional input files, fastq/bam/vcf
- hg38 alignment and annotations
- profiles, wgs/onco/exome
- several new sv-variant callers, melt, cnvkit, delly
- POD-tool for duplication events in trios
- Freebayes calling for difficult homopylomers in onco-samples
- Yaml-creation for scout import overhauled
- new container with needed software
- new rank-models for onco (both snv and sv) and wgs (sv-rank not live yet)
- Various small improvements of code

#### Fixes
- Optimization of cpu/memory/time for each process
- Numerous small improvements of several scripts


### 1.5.4

#### Features
- Last hg19 version

#### Fixes
- Minor file-paths issues resolved

### 1.5.3

#### Fixes
- Fix incorrect filtering of UPD calls in genomeplotter

### 1.5.2

#### Fixes
- Only plot UPDs in overview plot if > 100 informative sites
- Don't run UPD process on duos

### 1.5

#### Features

#### Fixes
- Use the correct scout server for create_yaml and loqusdb annotation
- Removed some hardcoded assumptions in create_yml.pl
- Always create gvcfs, and publish
- Publish chanjo coverage file
- Fix ID mixup in gatkcov process
- Added retry strategies to cdm, locuscollector and bqsr processes
- create_yml: Only add each panel once and sort panels alphabetically

### 0.1.4

#### Features

#### Fixes
- Increase allocated  memory for Clinvar SnpSift process due to occasional crashes
- Further fixes to output folders
- Retry create_yaml and loqus process up to 5 times

### 0.1.3

#### Features
- Allow for non-distributed BWA (new default). For urgent cases use --shardbwa
- Add diagnosis field to CDM

#### Fixes
- Change back to adding intersected vcf for loqusdb instead of full genomic vcf
- Change name of 1000G INFO field to make it show up in Scout
- Removed some hardcoded paths to scripts and added them the the bin/ folder

### 0.1.2

#### Features

#### Fixes
- Add 1000G in a special field for positions missing gnomAD (typically non-exonic) and add it to rankmodel

### 0.1.1

#### Features

#### Fixes
- Properly add selected gene panels to YAML
- Add STR-vcf to YAML
- Rename sample in STR-vcf to agree with sample name instead of bam filename
