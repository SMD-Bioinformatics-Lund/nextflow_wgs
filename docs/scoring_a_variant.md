# Score a single SNV

Run `score_variant.nf` with the annotation resources, containers, and rank models used by the WGS pipeline. The standalone workflow takes one GRCh38 SNV and requires no sample sheet or sequencing reads.

```bash
printf '11\t5227002\t.\tT\tA\t.\t.\t.\n' > snv.vcf
nextflow run score_variant.nf -c configs/nextflow.hopper.config -profile wgs \
    --variant snv.vcf --score_outdir snv-score
```

The Hopper configuration requires access to the cluster's annotation resources and Singularity containers. `genome_file` must point to the uncompressed GRCh38 FASTA and its adjacent `.fai` index.

Provide a plain or gzip-compressed VCF containing exactly one record. A headerless row with the eight tab-separated VCF columns also works.

- Give a single A/C/G/T REF base and a different single A/C/G/T ALT base.
- The Python step rebuilds the record from its coordinates and alleles, generates a singleton pedigree, and discards input sample columns and existing annotations. The default genotype is `0/1` with unknown sex. Set `--query_genotype 1/1` for a homozygous alternate call and `--query_sex 1` or `2` for male or female, respectively.

The output directory contains `scored.vcf` with annotations and the Genmod `RankScore`, plus `score.tsv` with one row containing chromosome, position, end, alleles, and the numeric score. The default output directory is `variant-score`. Specify a different directory for each query whose results you want to keep.

The workflow evaluates the SNV with the configured singleton rank model (`rank_model_s`) using a synthetic sample. The score omits family context and read-based adjustments. Its value can differ from the full WGS pipeline's score for a patient call.

To generate the annotation input VCF and pedigree on their own, run:

```bash
python3 bin/prepare_variant_query.py --variant snv.vcf \
    --reference /path/to/GRCh38.fa --output-dir prepared
```
