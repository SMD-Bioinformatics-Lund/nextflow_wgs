# Temporary TODO: gnomAD SV annotation

Revisit after the first comparison with the previous VEP 114.2 pipeline.
The behaviours below currently preserve VEP compatibility; they are follow-up
decisions, not confirmed bugs. Keep them separate from the bedtools migration
and the fix for sequence-resolved matching.

- [ ] Review missing-value formatting. Literal `.` values are retained, while
  absent INFO fields are omitted from annotation lists. Decide whether a more
  consistent representation is useful and check its effect on downstream scoring.
- [ ] Review reference IDs. Currently only the first semicolon-separated ID is
  reported. Missing IDs, including the literal ID `0`, become coordinates using
  the source chromosome name. Consider retaining all IDs and treating `0` as an ID.
- [ ] Review escaping. VEP converts commas and pipes to `&`, whitespace to `_`,
  semicolons to `%3B`, and a standalone `-` to an empty value. Consider lossless
  escaping of identifiers once downstream CSQ readers have been checked.
- [ ] Review INFO flags in annotation fields. Flags currently become `1`, as in
  VEP. Consider validating unexpected flag values in numeric fields such as AF.
- [ ] Review comma-separated INFO values. VEP's heuristic discards the first
  value when the list contains one more value than the number of ALT alleles.
  Check whether using the VCF header's Number definition would be more appropriate.
- [ ] Review the restored gnomAD INFO header declarations. They match VEP's
  previous headers, although the annotation values are carried inside CSQ.
  Decide whether to simplify these headers and add explicit matcher provenance.
- [ ] Decide whether replacing pre-existing gnomAD CSQ fields is needed. This
  could be simplified if fresh VEP input without those fields is guaranteed.
- [ ] Check which symbolic variant encodings occur in the actual inputs before
  simplifying the type mappings. Preserve the existing same-type matching rule.

For any later change, compare annotation values and resulting rank scores on the
same data before accepting an output difference.
