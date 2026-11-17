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
| `GNOMAD_SV_INTERVALS`     | SV annotation        | directory       | Three matching tables prepared from the gnomAD SV VCF for bcftools annotation. |
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

### Preparing gnomAD SV annotations

VEP provides gene consequences and variant classes. gnomAD SV matching is handled
entirely by bcftools and the coordinate helper, with results written once per
variant in INFO. Ordinary deletions require at least 70% overlap in each direction,
including sequence alleles. Insertions retain master/VEP 114.2's two matching
paths, described below. Duplications, inversions, copy-number variants, tandem
and mobile-element classes also use same-class matching; there is no gnomAD custom
annotation in VEP.

Prepare the reference directory once from the original gnomAD SV VCF:

```bash
python3 bin/gnomad_sv.py reference \
    --vcf /path/to/gnomad.v4.1.sv.sites.slim.vcf.gz \
    --output /path/to/gnomad.v4.1.sv.sites.slim.intervals
```

Set `params.GNOMAD_SV_INTERVALS` to that directory. Rebuild it whenever the source
VCF changes. It contains `overlap.bed.gz`, `ins_site.bed.gz` and
`ins_length.bed.gz`. Rebuild all three tables when moving from the earlier
DEL/INS-only prototype: the payload now includes the reference SV class.
bcftools reads these gzipped BED tables directly in region-annotation mode;
no tabix indexes are needed. The helper rejects candidates of a different class.
The original gnomAD VCF is only needed for this preparation step.

The pipeline uses three small stages with the existing Python and bcftools
containers: prepare temporary matching VCFs, annotate them with bcftools, and copy
the results into the original VCF's gnomAD INFO fields. The temporary records use
input row numbers as IDs, allowing duplicate or missing variant IDs. Original
coordinates, alleles, CSQ and genotypes are retained.

`gnomad` and `gnomad_SVTYPE` are String lists; `gnomad_AF` and
`gnomad_FREQ_HOMALT` are Float lists. Values use standard comma separators, with
`.` for missing fields to preserve alignment between matches. Records without
matches have no gnomAD INFO values. Both live SV rank models read
`INFO/gnomad_AF` and retain their existing maximum-AF aggregation and thresholds.
The Hopper profiles use these models from `rank_models/live/sv` in the pipeline
checkout. Explicit model overrides must also read the new INFO field.

#### Exact parsing and matching rules

The INS compatibility rules follow VEP 114.2's
[input parser](https://github.com/Ensembl/ensembl-vep/blob/release/114.2/modules/Bio/EnsEMBL/VEP/Parser/VCF.pm),
[VCF matcher](https://github.com/Ensembl/ensembl-vep/blob/release/114.2/modules/Bio/EnsEMBL/VEP/AnnotationSource/File/VCF.pm)
and [reference scan](https://github.com/Ensembl/ensembl-vep/blob/release/114.2/modules/Bio/EnsEMBL/VEP/AnnotationSource/File.pm).

First read the CSQ header and select supported SV classes using `VARIANT_CLASS`.
`INS:ME:ALU`, for example, belongs to a different class and must not be collapsed
into ordinary INS. Skip records VEP left without CSQ. Normalise chromosome names
using the synonyms and the `chr`/`M` aliases. Do not split multiallelic records or
guess a class for complex replacements. VEP's `indel` class and unsupported
reference types such as CPX/CTX are not coerced into a different class.

The rules below use **zero-based, half-open intervals**. `P` is the original
one-based VCF POS. For sequence alleles, remove their common prefix and then
their common suffix, without realigning against a FASTA. If the prefix had `k`
bases, define `B = P - 1 + k`.

| Record | Parse | Matching representation and rule |
| --- | --- | --- |
| Sequence DEL | After trimming, REF remains and ALT is empty. Let `L` be the remaining REF length. | `[B, B + L)`, reciprocal 70%. |
| Symbolic DEL, or DEL with unknown bases and `SVTYPE=DEL` | Use the first `abs(SVLEN)`; otherwise use END. | `[P, P + L)` or `[P, END)`, reciprocal 70%. Missing or non-positive spans cannot be annotated. |
| Sequence INS handled as a regular VEP variation feature | After trimming, REF is empty and ALT remains. | Query site `[B - 1, B)`. Compare with the reference VCF's raw REF span, without an inserted-length check. With gnomAD's one-base REF this requires the same insertion position. |
| Symbolic INS, or INS routed through VEP's structural path | A non-ACGT allele plus `SVTYPE=INS` takes this path, including Manta's mixed records with N-padded REF. Set `S = P + 1`; set `T = P + abs(first SVLEN)`, otherwise END, otherwise `P + 1`. | `[min(S,T) - 1, max(S,T))`, reciprocal 70%. This is VEP's matching convention, not a claim that the insertion occupies that reference span. |
| Other structural classes | Classify symbolic ALT, retaining tandem and mobile-element subtypes. Use SVTYPE for non-standard ALT spellings, as VEP does (`<TDUP>` with `SVTYPE=DUP` is ordinary duplication). Use SVLEN before END. | Positive metadata span `[P,T)`, reciprocal 70%, same class only. Structural mobile-element INS uses the INS-length path above. Query breakends ignore a mate END when SVLEN is absent. |

In particular, `A -> ATTT` is a sequence insertion at the site after the A. Its
three inserted bases do **not** extend the legacy site-matching window. Its
SVLEN is ignored for matching on this path. `ATTT -> A`, by contrast, becomes a
three-base deleted interval after removing the unchanged A.

For structural INS, SVLEN takes precedence over END. END often differs from
`P + SVLEN` in this gnomAD file. When SVLEN is absent, preserve VEP's END boundary
rather than inventing an insertion length. For example, `END=P` gives the
two-base boundary `[P-1,P+1)`; if END is also absent, the boundary is `[P,P+1)`.

For classes other than ordinary deletion, a hit must have been reachable by
VEP's original indexed scan. Preserve the reference's original index end
(`END`, or `POS - 1 + len(REF)` when absent). Require that end to be greater than
the query's original VEP seek start, and require the reference's adjusted start
to be no later than the query's original VEP end plus one. For site queries,
the seek start and VEP end are both `B`; for structural queries, they are `P`
and `T`. This prevents adding matches VEP could not have fetched even when its
SVLEN-derived intervals overlap. Reference raw REF spans go in `ins_site`;
positive VEP metadata spans go in `ins_length`.

For interval matching, bcftools retrieves candidates with a small numerical
tolerance; the helper then checks exactly 70% in both directions using integer
arithmetic. For sequence INS it uses unrestricted site overlap instead. No
inserted-sequence identity comparison is performed. In all modes, sort accepted
hits by original reference row and retain the first 50. If there are more hits,
append `...` to the ID list. Ranking continues to select the maximum reported AF.

For local tests without VEP, `gnomad_sv.py query` also accepts a raw benchmark
VCF. The same INFO output is used with and without CSQ. Complex replacement
alleles that cannot be represented as a pure deletion or insertion are left
unchanged with a message; use VEP classification for those records. These commands
annotate the VCF; run scoring with the updated INFO-based SV rank model afterwards.

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
