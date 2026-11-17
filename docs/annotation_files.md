# Annotation files

To run the pipeline, you need to setup a number of annotation files. Default values are specified in the `configs/nextflow.hopper.config`. A summary of the required annotation files is shown below.

| Parameter                 | Function             | Format          | Description                                                                                                                                                                                                                                             |
| ------------------------- | -------------------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `genome_file`             | Reference genome     | `fasta`         | The main reference                                                                                                                                                                                                                                      |
| `GENOMEDICT`              | Reference genome     | `dict`          | Dict file for the FASTA file built using [Picard Tools](https://broadinstitute.github.io/picard/)                                                                                                                                                       |
| `rCRS_fasta`              | Reference genome     | `fasta`         | Mitochondrial FASTA reference sequence found [here](https://www.ncbi.nlm.nih.gov/nuccore/251831106)                                                                                                                                                     |
| `bwa_shards`              | Alignment            | integer         | Number of shards to split the reads into prior to alignment                                                                                                                                                                                             |
| `shardbwa`                | Alignment            | boolean         | Boolean specifying whether to do alignment in sharded mode                                                                                                                                                                                              |
| `bqsr_known_polymorphic_sites_vcf` | Alignment          | `vcf`           | Gold-standard indels used in Sentieon's base quality score recalibration (BQSR)                                                                                                                                                                         |
| `VEP_CACHE`               | Annotation (VEP)     | folder          | VEP files for offline run. Instructions on how to setup a cache can be found [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache).                                                                                           |
| `VEP_FASTA`               | Annotation (VEP)     | `fasta`         | Reference sequence used for optimization within the VEP cache                                                                                                                                                                                           |
| `SYNONYMS`                | Annotation (VEP)     | `tsv`           | Chromosome synonyms used by VEP (for instance to recognize `M` as `MT`)                                                                                                                                                                                 |
| `CADD`                    | Annotation (VEP)     | `tsv`           | Precalculated CADD indices ([download page](https://cadd.gs.washington.edu/download)). [CADD](https://cadd.gs.washington.edu) is used to score deleteriousness of SNVs and indels in human.                                                             |
| `GNOMAD_GENOMES`          | Annotation (VEP)     | `vcf`           | Allele frequencies for variants ([download page](https://gnomad.broadinstitute.org/downloads))                                                                                                                                                          |
| `GNOMAD_EXOMES`           | Annotation (VEP)     | `vcf`           |                                                                                                                                                                                                                                                         |
| `GNOMAD_MT`               | Annotation (VEP)     | `vcf`           |                                                                                                                                                                                                                                                         |
| `GNOMAD_SV`               | Annotation (VEP)     | `vcf.gz`        | Original gnomAD SV reference; VEP continues to use this file. |
| `GNOMAD_SV_DEL`           | Deletion annotation  | `bed.gz`        | Deletion intervals prepared from `GNOMAD_SV` for reciprocal overlap matching with bcftools. |
| `MAXENTSCAN`              | Annotation (VEP)     | folder          | A direct path to the MaxEntScan scripts                                                                                                                                                                                                                 |
| `dbNSFP`                  | Annotation (VEP)     | `tsv`           | Path to gzipped tab delimited file containing the dbNSFP annotations.                                                                                                                                                                                   |
| `PHYLOP`                  | Annotation (VEP)     | `tsv`           | Measures evolutionary conservation at individual alignment sites. This parameter expects a tab-delimited file with pre-calculated phyloP scores.                                                                                                        |
| `PHASTCONS`               | Annotation (VEP)     | `tsv`           | Conservation scores. This parameter expects a tab-delimited file with pre-calculated phyloP scores.                                                                                                                                                     |
| `vcfanno_config`         | Annotation (VCFAnno) | `toml`          | Config file pointing to which files to retrieve annotation from, which fields that should be extracted and how these should be inserted in the target VCF. This additional annotations are different depending on which profile is used for processing. |
| `vcfanno_lua`            | Annotation (VCFAnno) | `lua`           | Custom script for VCFAnno annotations                                                                                                                                                                                                                   |
| `COV_INTERNAL_LIST`       | CNV                  | `interval_list` | List of pre-calculated intervals to use in `gatk CollectReadCount` ([more info](https://gatk.broadinstitute.org/hc/en-us/articles/360037427371-PreprocessIntervals-BETA-))                                                                              |
| `GATK_PON_MALE`           | CNV                  | `hdf5`          | Gender specific panel of normals used in `DenoiseReadCounts` GATK step                                                                                                                                                                                  |
| `GATK_PON_FEMALE`         | CNV                  | `hdf5`          | Gender specific panel of normals used in `DenoiseReadCounts` GATK step                                                                                                                                                                                  |
| `expansionhunter_catalog` | STR                  | `json`          | Catalog of variants to look for and their genomic locations. A default catalogue can be downloaded [here](https://github.com/Illumina/ExpansionHunter/tree/master/variant_catalog/grch38)                                                               |
| `rank_model`              | Rank models          | `ini`           | Rank model for SNVs, trios                                                                                                                                                                                                                              |
| `rank_model_s`            | Rank models          | `ini`           | Rank model for SNVs, singles                                                                                                                                                                                                                            |
| `svrank_model`            | Rank models          | `ini`           | Rank model for SVs, trios                                                                                                                                                                                                                               |
| `svrank_model_s`          | Rank models          | `ini`           | Rank model for SVs, singles                                                                                                                                                                                                                             |
| `gene_regions`            | Onco                 | `bed`           | Ranges with genes of interest specific for onco mode                                                                                                                                                                                                    |
| `FASTGNOMAD_REF`          | UPD & ROH            | `dat`           | Compact file with highly frequent SNV-locations. Used for filtering.                                                                                                                                                                                    |
| `GENS_GNOMAD`             | Gens                 | `tsv`           | Gnomad frequencies for Gens                                                                                                                                                                                                                             |
| `loqusdb`                 | LoqusDB              | string          | Path-name to the `loqusdb` from which to retrieve local artifact frequencies                                                                                                                                                                            |
| `panelsdef`               | Scout                | `json`          | Location of panels used in Scout                                                                                                                                                                                                                        |
| `scoutbed`                | Scout                | `bed`           | Bedfile used in coverage calculations for Scout                                                                                                                                                                                                         |
| `svdb`                    | LoqusDB              | `vcf`           | LoqusDB SV 
GATK
`gatk_ref_csv`              | GENS/gatk cnvs       | `csv`           | GENS pons per sex and input ploidy, call models for GATK coverage and interval files |

artefacts                                                                                                                                                                                                                                    |
## Preparing annotation files

This is a non-comprehensive documentation on how some of the reference files can be prepared. If you miss instructions for some annotation file, please open an issue in the [GitHub repository](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/issues).

### Preparing gnomAD deletion annotations

VEP continues to annotate SVs using the original `GNOMAD_SV` VCF. A subsequent
bcftools step replaces only the annotations for VEP's ordinary `deletion` class,
requiring at least 70% overlap in each direction. Insertions, duplications,
mobile-element deletions and other classes keep their existing VEP annotations.

Prepare the deletion reference once, using the same gnomAD release as `GNOMAD_SV`:

```bash
python3 bin/gnomad_del.py reference \
    --vcf /path/to/gnomad.v4.1.sv.sites.slim.vcf.gz \
    --output /path/to/gnomad.v4.1.sv.sites.slim.del.bed.gz
```

Set `params.GNOMAD_SV_DEL` to the output file. Keep `params.GNOMAD_SV` pointing to
the original VCF. Rebuild the deletion reference whenever the original changes.
The generated file is a gzipped BED with a payload column; bcftools' region
annotation mode reads it directly, so no tabix index is needed.

The pipeline uses three small stages with the existing Python and bcftools
containers: prepare temporary deletion intervals, annotate them with bcftools,
and copy the results into the original VCF's gnomAD CSQ fields. Removing the VCF
anchor base in the temporary copy makes the overlap calculation independent of
whether the deletion has symbolic or sequence alleles. The original coordinates,
alleles, other annotations and genotypes are retained. Match order and VEP's
50-record limit are preserved.

bcftools retrieves candidates with a small numerical tolerance; the helper then
checks exactly 70% in both directions using integer arithmetic. This prevents
floating-point rounding at the threshold from changing long-variant matches.

For local tests without VEP, `gnomad_del.py query` also accepts a raw benchmark
VCF, and `restore` writes deletion annotations as INFO fields. Complex replacement
alleles that cannot be represented as a pure deletion are left unchanged with a
message; use VEP classification for those records. These matching tests do not
replace a complete pipeline evaluation or recompute rank scores.

### Preparing reference and indexes

The pipeline is running with the hg38 reference genome, but with no `chr` prefix. 

The `chr` prefix from each FASTA entry.

```
sed 's/^>chr/>/' <fasta>.fna > <fasta>_nochr.fna
```

Generate a bwa index. This requires that you have [bwa](https://github.com/lh3/bwa) installed.

```
bwa index <fasta>_nochr.fna
```

Create the sequence dictionary. This requires [Picard Tools](https://broadinstitute.github.io/picard/).

```
java -jar picard.jar CreateSequenceDictionary \
    REFERENCE=<fasta>_nochr.fna \
    OUTPUT=<fasta>_nochr.dict
```

Generate the samtools fasta index using [samtools](http://www.htslib.org/).

```
samtools faidx <fasta>_nochr.fna
```

### VCF anno

Example `vcf_anno` as used by CMD. Replace `<base_dir>` with the base folder for each file. Here the following information is retrieved:

* Frequency information from `loqusdb`
* Clinical significance information and reviewer quality from `ClinVar`
* Ensemble gene annotations

```
[[annotation]]
file="<base_dir>/loqusdb_latest.vcf.gz"
fields=["Frq"]
names=["loqusdb_freq"]
ops=["concat"]

[[annotation]]
file="<base_dir>/clinvar38_latest.vcf.gz"
fields=["CLNSIG","CLNREVSTAT","CLNACC"]
names=["CLNSIG","CLNREVSTAT","CLNACC"]
ops=["first","first","first"]

[[annotation]]
file="<base_dir>/ensembl_genes_38.sort.padded200bp.noNeg.bed.gz"
columns=[4]
names=["Annotation"]
ops=["concat"]
```

### GATK meta-csv

Example of meta-csv for gatk- references/pons. Headers (description). PON can be left empty for panels that dont use GENS. if platform/sex combination does not exist in this meta-csv pipeline will exit with error message to indicate this

* platform (illumina/illuminax)
* sex (M/F)
* intervals (path_to_interval_file)
* ploidy_model (path_to_ploidy_model)
* ref_folders (path_to_ref_folders)
* pon (path_to_sex_specific_gens_pon)
