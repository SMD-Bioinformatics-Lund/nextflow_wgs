# Production health audit: P0 clinical-risk findings

Audit date: 2026-09-01  
Audited revision: `4d79255bc5d74c312e2978ccd449bcfcdb4e2190`  
Branch at audit: `alkc/migrate-cftr-workflow-input`

This document records only the audit findings classified as P0: defects on an
active clinical path that can directly suppress clinically relevant variants or
materially misrepresent evidence used for interpretation. It intentionally
excludes reliability-only findings such as possible workflow crashes.

The code behavior below is confirmed. Whether any historical patient was
affected requires retrospective case analysis; that must not be inferred from
the code review alone.

## P0-1: DPAF deletes a family-wide variant when any sample fails depth or allele-fraction thresholds

### Evidence

`bin/wgs_DPAF_filter.pl:20-38` iterates over every sample column in a VCF
record. It drops the entire record if any sample has missing/zero DP, or if any
non-reference genotype has fewer than three alternate reads or VAF below 0.2.
The filter also assumes the FORMAT fields are ordered `GT:AD:DP` rather than
resolving their positions from the FORMAT column.

`workflows/split_normalize_snvs.nf:10-15` applies this filter to the combined
family VCF before target intersection, annotation, inheritance modeling, and
ranking. Thus, a good proband call can disappear solely because an unaffected
parent has zero depth or another family member has a weak call.

This was reproduced during the audit with a trio record containing:

```text
proband  0/1:10,10:20
mother   0/0:0,0:0
father   0/0:20,0:20
```

The output retained the VCF header but removed the variant record.

### Clinical mechanism

This is a false-negative mechanism. A potentially diagnostic proband variant
can be absent from all downstream annotated and ranked results because of a
relative's sequencing state, rather than because the proband evidence failed.

Simply exempting `0/0` relatives from the record filter is insufficient. A
zero-depth parent left as `0/0` can falsely support a de novo interpretation.
Record retention, per-sample genotype quality/no-calling, and inheritance
classification must be handled separately.

### Verify on a minimal fixture

Create a small, valid trio VCF with the genotype combination above and run:

```bash
perl bin/wgs_DPAF_filter.pl trio.vcf > trio.dpaf.vcf
grep -v '^#' trio.dpaf.vcf
```

The current implementation produces no data record. Repeat with a parent whose
genotype is missing and with FORMAT fields reordered; these exercise the other
unsafe assumptions.

### Verify on real cases

For a completed family analysis, locate the exact pre- and post-filter files in
the Nextflow work directory:

- pre-filter: the `vcflib_vcfuniq` output, normally `<group>.uniq.vcf.gz`
- post-filter: the `wgs_dpaf_filter` output, normally
  `<group>.norm.uniq.DPAF.vcf.gz`

Compare them without changing the production run:

```bash
bcftools isec -p isec-output \
  '<group>.uniq.vcf.gz' \
  '<group>.norm.uniq.DPAF.vcf.gz'

bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT:%AD:%DP]\n' \
  isec-output/0000.vcf
```

`isec-output/0000.vcf` contains records present before DPAF and absent after it.
For each removed record, independently calculate whether the proband meets the
intended support threshold and identify which relative caused record deletion.
Then confirm whether the record is absent from the final scored VCF and any
case-delivery artifact.

For historical impact assessment, perform this comparison across prior trio or
family runs, constrained to the clinically interpreted target regions and using
the annotation sources appropriate to the original analysis date. Prioritize
removed variants where the proband passes support thresholds and the deletion
is caused only by a relative's no-call, zero depth, or weak alternate evidence.

### Required approach

1. Freeze or add an explicit review gate for affected family workflows until
   the behavior is clinically dispositioned.
2. Define the intended policy with laboratory and clinical stakeholders:
   normally, record retention should depend on a supported non-reference call
   in a clinically relevant family member, not acceptable coverage in every
   sample.
3. Convert insufficient per-sample evidence to an explicit missing/uncertain
   genotype before inheritance modeling; do not silently convert it to
   reference.
4. Parse FORMAT by field name and explicitly support phased, haploid,
   multiallelic, and missing genotypes.
5. Add regression cases covering a strong proband call with zero-depth,
   missing, and weak-call relatives; reordered FORMAT; phased/haploid calls;
   multiallelic sites; and a weak singleton call that should be removed.
6. Validate the replacement on synthetic families, reference truth sets, and a
   representative set of real families. Review changes in inheritance labels,
   not only variant counts.
7. Use the retrospective comparison to decide whether historical cases require
   reanalysis or amendment.

Release acceptance criterion: no supported clinically relevant family-member
call may be removed solely because another family member lacks adequate data,
and inadequate relative evidence must not be interpreted as reference evidence.

## P0-2: Panel SV post-processing changes observed `0/0` genotypes to `0/1`

### Evidence

`bin/modify_cnv_genotypes_for_loqusdb.pl:31-34` rewrites the first sample field
from a leading `0/0` to `0/1`.

Despite the script name and comment, the altered VCF is not limited to LoqusDB:

- `main.nf:693-699` sends the original `svdb_merge_panel.out.loqusdb_vcf` to
  the LoqusDB channel, while the other merged VCF enters post-processing.
- `main.nf:2963-2996` removes BNDs, performs the genotype rewrite, adds MELT
  data, and emits the rewritten VCF.
- `main.nf:732-755` sends that result through proband filtering, AnnotSV, VEP,
  artefact annotation, dbVar annotation, and later scoring.
- `main.nf:3216-3219` only removes a proband `./.` call; it does not restore or
  distinguish a fabricated heterozygous call.

The `0/0` to `0/1` transformation was reproduced during the audit with a
minimal panel SV record.

### Clinical mechanism

The final clinical SV path can represent a reference genotype as a
heterozygous event. This can create a false-positive candidate, alter zygosity
and inheritance assessment, and make downstream evidence appear stronger than
the caller output supports.

### Required approach

1. Place an immediate review gate on panel SV results produced through this
   path, particularly calls whose upstream genotype was `0/0`.
2. Trace several current cases from the caller VCF through the merged,
   genotype-fixed, annotated, scored, and delivered artifacts; compare genotype
   fields at each stage.
3. Remove clinical-path genotype fabrication. If a database loader requires a
   special representation, create a separate loader-only artifact with explicit
   provenance that can never feed annotation, scoring, or reporting.
4. Retrospectively identify delivered variants whose clinical-path genotype
   differs from the upstream caller genotype and clinically review them.
5. Add an invariant test: clinical output genotype must equal the evidence-based
   genotype unless a documented genotype-refinement algorithm records why it
   changed.

Release acceptance criterion: a loader compatibility transformation cannot
enter any clinical interpretation artifact, and every reported genotype must be
traceable to caller evidence or a validated refinement step.

## P0-3: Mitochondrial calls bypass Mutect2 filtering and receive fabricated genotype quality

### Evidence

- `main.nf:1628-1656` runs Mutect2 in mitochondrial mode but does not run
  `FilterMutectCalls`; the source comment explicitly defers it.
- `main.nf:1694-1704` unconditionally removes position `M:955`, filters mainly
  on sample AF above 0.05, and passes the result to
  `filter_mutect2_mito.pl`.
- `bin/filter_mutect2_mito.pl:30-40` retains proband calls with DP greater than
  50, while `bin/filter_mutect2_mito.pl:46-59` assigns placeholder GQ 99 to
  every individual and changes genotype to `1/1` when AF is at least 0.99.
- `main.nf:377-395` annotates these calls and merges them into the main SNV VCF
  path.

### Clinical mechanism

Raw/unfiltered mitochondrial candidates can enter the clinical variant set with
an invented high-confidence genotype value. This can conceal the lack of
validated filtering from downstream systems and reviewers. The unconditional
position deletion is also a site-specific false-negative mechanism independent
of the current evidence in a sample.

### Required approach

1. Flag mitochondrial results from this path for explicit manual review; do not
   use placeholder GQ as evidence.
2. Establish and validate a mitochondrial calling/filtering workflow, including
   Mutect2-recommended filtering inputs where applicable, heteroplasmy limits,
   strand/read-position artifacts, NUMT handling, depth, and assay-specific
   performance.
3. Preserve caller-native confidence/filter fields. If downstream software
   requires GQ, represent unavailable quality as missing rather than 99, or use
   a validated conversion with provenance.
4. Replace the unconditional `M:955` deletion with a validated, documented
   artifact rule that is auditable and cannot discard a true call solely by
   coordinate.
5. Review historical mitochondrial candidates near reporting thresholds and at
   the excluded coordinate, using orthogonal evidence where needed.

Release acceptance criterion: only calls passing a validated mitochondrial QC
model enter clinical interpretation, no confidence value is fabricated, and
coordinate-specific suppression has documented validation and traceability.

## P0-4: Oncology FreeBayes rescue calls are merged with fabricated read support

### Evidence

- `workflows/call_snvs.nf:44-66` enables the FreeBayes branch when configured
  and concatenates its output into the main combined VCF.
- `workflows/call_snvs.nf:265-273` calls pooled FreeBayes with minimum alternate
  fraction 0.03, selects records by a case-insensitive `pathogenic` text match,
  and runs `filter_freebayes.pl`.
- `bin/filter_freebayes.pl:17-35` accepts selected pathogenic/likely pathogenic
  records and emits `GT:AD:DP:GQ:PGT:PID:PL`, but hard-codes
  `AD=100,100`, `DP=200`, `GQ=0`, and zero-valued phasing/likelihood fields.
- `workflows/call_snvs.nf:296-314` concatenates these generated records with the
  primary gVCF and sorts them into the combined VCF used downstream.

During the audit, a synthetic low-support FreeBayes record (QUAL 3, AF 0.03,
DP 12) passed the script when marked pathogenic and emerged with fabricated
`AD=100,100` and `DP=200`.

### Clinical mechanism

A low-evidence oncology rescue call can appear in the clinical pipeline with
read depth and allelic balance that never existed. Downstream filters and human
review can therefore overestimate analytical support, creating a direct
false-positive/reporting risk. The text match can also select annotations
without a structured, exact clinical-significance policy.

### Required approach

1. Identify analyses with the FreeBayes branch enabled and flag rescue variants
   for review against the original FreeBayes VCF/BAM evidence.
2. Stop fabricating FORMAT evidence. Preserve caller-native AD, DP, quality, and
   filter fields, or leave unavailable fields missing.
3. Replace free-text `grep -i pathogenic` selection with parsed, versioned
   annotation fields and an approved significance policy.
4. Validate low-VAF performance with assay-specific truth material across depth,
   allele fraction, variant class, and difficult sequence context.
5. Retrospectively compare delivered rescue calls with native FreeBayes evidence
   and decide which cases need clinical reassessment.

Release acceptance criterion: every value presented as read support is derived
from the source reads/caller output, rescue eligibility is structured and
versioned, and validated analytical thresholds gate clinical inclusion.

## Immediate priority order

1. DPAF family-wide record deletion: direct diagnostic false negatives and the
   first issue selected for investigation.
2. Panel SV genotype fabrication: reference calls can become clinical
   heterozygous candidates.
3. Oncology FreeBayes evidence fabrication: low-support calls can acquire false
   high depth/allelic balance.
4. Mitochondrial filtering and quality fabrication: unvalidated calls and a
   coordinate-specific deletion enter the final SNV path.

All four require clinical ownership, retrospective impact assessment, and
release-level regression evidence. A code correction alone does not establish
that historical cases were unaffected.
