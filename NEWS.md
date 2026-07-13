# ggmethylation 0.3.0

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
- `parse_mm_ml()` now honours the trailing `?`/`.` flag on MM tag entries:
  `.` emits explicit `mod_prob = 0` rows for unlisted canonical positions
  (implicit-unmodified); `?` or no flag continues to omit them.
- `theme_ggmethylation()` exposes the shared panel theme.
- `plot_methylation()` gains `show_ci` (default `TRUE`) to draw a shaded loess
  confidence-interval ribbon behind the smoothed modification probability
  line(s) in the bottom panel, using the `lower`/`upper` columns from
  `smooth_methylation()`.
- `plot_methylation()` gains `call_mode` (`"continuous"` default or
  `"binary"`), `call_threshold` (default `0.5`), and `call_ambiguous`
  (default `NULL`) to render read-panel modification sites as discrete
  methylated/unmethylated/ambiguous calls instead of a continuous
  probability gradient.
- `plot_methylation()` gains `show_delta` (default `FALSE`). When grouping
  yields exactly two groups, setting `show_delta = TRUE` appends a bottom
  panel showing the signed difference in loess-smoothed modification
  probability between the two groups (group2 - group1), rendered as a
  diverging area coloured by sign. Not supported for multi-sample data.

## Bug fixes

- `.insert_deletion_breaks()` (internal) now also nulls the `lower`/`upper`
  confidence-interval columns inside consensus deletion gaps, matching the
  existing `mean_prob` masking, so the CI ribbon does not bridge across
  deletion breaks.

## Changes

- Default `group_colours` updated to a colorblind-safe Okabe-Ito blue/orange
  pair. Pass an explicit `group_colours` vector to restore prior colours.

## Internal changes

- `.smooth_xy()` helper extracted from `smooth_methylation()` for reuse by
  `plot_insertion_locus()`.
