# ggmethylation 0.2.0

## Breaking changes

- `parse_mm_ml()` now returns a named list with two data frames (`$sites` and
  `$insertion_sites`) instead of a single data frame. Code calling this
  internal function directly must be updated.

## New features

- `read_methylation()` now retains modification calls on inserted bases (CIGAR
  `I` operations). These are stored in `methylation_data$insertion_sites` with
  columns `read_name`, `ref_anchor`, `query_pos`, `ins_offset`, `ins_length`,
  `mod_prob`, and `mod_code`. The `$cigar_features` data frame is also
  populated for all reads.
- `insertion_sites(m)` — convenience accessor for `m$insertion_sites`.
- `list_insertion_loci(m, tol_pos, tol_len, min_reads)` — discovers recurrent
  insertion loci across reads using a greedy single-pass clustering algorithm.
  Returns a summary data frame with carrier counts, length statistics, and
  mean modification probability per locus.
- `plot_insertion_locus(m, locus_id, ...)` — visualises modification
  probabilities at a single insertion locus in a stitched coordinate system
  (left reference flank | insertion bases | right reference flank), with an
  optional loess-smoothed comparison panel for carrier vs. non-carrier reads.

## Internal changes

- `.smooth_xy()` helper extracted from `smooth_methylation()` for reuse by
  `plot_insertion_locus()`.
