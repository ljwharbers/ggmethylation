# ggmethylation 0.3.2

## Bug fixes

- Supplementary-alignment indicators were drawn on **both** ends of a read. The
  overlay chose its side from `$reads$clip_side`, which only asks whether the
  CIGAR begins or ends in `S`/`H` — and adapter trimming clips both ends of
  nearly every long read, so `clip_side` was `"both"` even for a read with a
  single supplementary partner. The side is now derived from read (query)
  coordinates by comparing the primary alignment's extent with each `SA`
  entry's, so an indicator appears only on the flank where a partner actually
  joins. `$reads` gains `sa_side` and per-flank `sa_chrom_left` /
  `sa_pos_left` / `sa_chrom_right` / `sa_pos_right`; `clip_side` is unchanged
  and still reports raw clipping. A read that genuinely spans two breakpoints
  still gets two indicators, each coloured by its own partner. Objects built
  by an earlier version fall back to the old clip-based sides with a warning.
- Reads split on large deletions (`show_cigar = TRUE`) received a
  supplementary indicator at every segment edge, including interior deletion
  boundaries. Indicators are now emitted only on the read's outer ends.
- The grouped read panel drew supplementary indicators 1.5× wider than the
  zone where modification ticks were suppressed, so ticks showed through the
  marker. Both now use one shared width.
- `pack_reads()`'s wider `clip_gap` was keyed on `clip_side`, so every
  adapter-trimmed read forced a 100 bp gap even though it carried no
  indicator. `plot_methylation()` now keys it on `sa_side`.

# ggmethylation 0.3.1

## Bug fixes

- `.gitignore`'s `test*` pattern also matched the `tests/` directory, so any
  **new** file under `tests/testthat/` was silently ignored by `git add`.
  Existing tests survived only because they predated the rule. The scratch-file
  patterns are now anchored to the repository root.
- Integration tests pointed at a BAM that no longer existed, so all three had
  been skipping every run and `read_methylation()` had no executing end-to-end
  coverage. They now run against a committed fixture at the MEG3 imprinted
  locus (`tests/testthat/fixtures/hg002_fiberseq_MEG3.bam`), built
  reproducibly by `data-raw/build_test_fixture.sh`.
- `group_colours` defaults to names `"1"`/`"2"` to match the `HP` tag, so any
  other grouping — SNV genotype (`"REF"`/`"ALT"`), a custom BAM tag — matched
  nothing: every group fell through to grey and ggplot2 emitted "No shared
  levels found between `names(values)` of the manual scale and the data's
  values". Grouping by SNV genotype therefore lost its colours by default. The
  palette is now resolved against the groups actually present, keeping exact
  name matches and filling the rest without ever reusing a colour.
- `sort_by` was never validated. An unknown column name produced `order(NULL)`
  → `integer(0)`, which dropped every read and rendered an empty plot with no
  error. Unknown columns now raise an error naming the valid ones.
- `show_ci` was accepted by the multi-sample renderer and silently discarded;
  the confidence ribbon is now drawn on the shared multi-sample smooth panel.
- The delta track and the multi-sample `show_delta` guard skipped via
  `message()`, which is invisible in scripts and knitr. Both now `warning()`.
- The smooth and delta panels' y-axis titles overflowed their (short) panels
  and collided at ordinary figure heights, rendering as
  "Δ methylationMean modification probability". Both titles are now shorter and
  set in a smaller size. (The 0.3.0 notes claimed this was fixed; it was not.)
- The vignette shipped with every chunk disabled because
  `inst/extdata/vignette_cache.rds` was never committed, so it rendered no
  figures at all. It now ships a cache and renders. Four latent defects in it
  were fixed: a `dot_size` argument removed in 0.2.x, `sort_by = "position"`
  (a `$sites` column, not a `$reads` one), inverted panel descriptions, and the
  wrong GitHub organisation in the install line.
- `data-raw/build_vignette_cache.R` used `meth_snv <- meth_hp` as a placeholder,
  which made the vignette's SNV-grouping section present haplotype-grouped data
  as allele-split. It now uses a real heterozygous SNV.

## Internal changes

- New `.validate_sort_by()` and `.resolve_group_colours()` helpers.
- Two tests were failing on `main` and now pass: `test-delta_track.R` assigned
  `c("gg", "ggplot")` over a plot's class vector to strip patchwork, which
  breaks dispatch under ggplot2 >= 4.0 (ggplot objects are S7); and
  `test-match_sa_to_vcf_bnd.R` still expected the `geom_vline` that `eab0c7d`
  removed from `build_bnd_layer()`.
- New test files for `write_methylation()`, `read_variants()`,
  `plot_insertion_locus()`, `.validate_sort_by()` and
  `.resolve_group_colours()`. `write_methylation()` and `read_variants()` had
  no direct coverage at all.
- `withr` added to `Suggests`.

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
- `plot_methylation()` gains `colour_ambiguous` (default `"#78909C"`) to set
  the read-panel colour of ambiguous calls, alongside the existing
  `colour_low`/`colour_high`. Only used when `call_mode = "binary"` and
  `call_ambiguous` is non-`NULL`.

## Bug fixes

- The delta panel's y-axis title was long enough to overflow its (short) panel
  and collide with the smooth panel's y-axis title above it. It is now
  `"Δ fraction methylated"` / `"Δ methylation"`; the `"(group2 - group1)"` line
  is dropped, since the fill colour now identifies which group is higher.
- Ambiguous calls were drawn in `grey75`, which is indistinguishable from the
  default `colour_low` (`"#BDBDBD"`) — setting `call_ambiguous` therefore did
  not actually make low-confidence calls visible. They now use a slate
  blue-grey deliberately off the `colour_low`/`colour_high` ramp.
- `.insert_deletion_breaks()` (internal) now also nulls the `lower`/`upper`
  confidence-interval columns inside consensus deletion gaps, matching the
  existing `mean_prob` masking, so the CI ribbon does not bridge across
  deletion breaks.

## Changes

- Default `group_colours` updated to a colorblind-safe Okabe-Ito blue/orange
  pair. Pass an explicit `group_colours` vector to restore prior colours.
- The `show_delta` panel now takes its fills from `group_colours` rather than
  a fixed diverging pair: the area is coloured by whichever group is higher at
  each position (positive = group 2, negative = group 1), so all panels agree
  on which colour means which group. Falls back to the previous diverging
  palette when `group_colours` is `NULL` or does not name both groups.

## Internal changes

- `.smooth_xy()` helper extracted from `smooth_methylation()` for reuse by
  `plot_insertion_locus()`.

# ggmethylation 0.2.0

Back-filled: this release shipped without a NEWS entry. Reconstructed from the
commit history for the record.

## New features

- `read_variants()` reads a VCF (via `VariantAnnotation`) and classifies each
  record as an SNV, insertion, deletion, or structural variant (`DEL`, `DUP`,
  `INV`, `BND`), returning a `variant_data` object. Passing it to
  `plot_methylation(variants = )` overlays the calls on the read panel: SNV
  asterisks on the reads carrying the ALT allele, spanning bars for structural
  variants, and mate labels for BND records. `bnd_match_tol` controls how far a
  read's supplementary alignment may sit from a VCF BND position and still
  count as a match.
- `read_annotations()` builds a gene model for the plotted region from a UCSC
  ncbiRefSeq GTF (downloaded and cached per genome, `"hg38"` or `"chm13"`), a
  user-supplied `gtf`, or a `TxDb`. `plot_methylation(annotations = )` draws it
  as an IGV-style gene track above the reads.
  `clear_annotation_cache()` drops the cached TxDb.
- `merge_methylation()` combines several `methylation_data` objects covering
  the same region into a `multi_methylation_data` object, which
  `plot_methylation()` renders as stacked per-sample read panels above one
  shared smooth panel. `print()` and `summary()` methods included.
- `write_methylation()` exports `$reads` and `$sites` to TSV or 6-column BED
  (0-based half-open, `mod_prob` as the score), optionally gzip-compressed.
- Structural variants taken from the CIGAR string are drawn on the read panel
  (`show_cigar`, `min_indel_size`), and reads with supplementary alignments are
  marked with arrowheads (`show_supplementary`).
- Consensus deletions break the smoothed curve rather than being interpolated
  across.
- The loess span adapts to the data when `smooth_span` is `NULL`.
- Insertion-aware modification calls: `read_methylation()` retains
  modifications on inserted bases in `$insertion_sites`, `list_insertion_loci()`
  clusters insertion events across reads, and `plot_insertion_locus()`
  visualises one locus in a stitched coordinate system.

## Bug fixes

- SNV asterisk overlay marked every read rather than only ALT carriers
  (issue #20).
- The SNV overlay used the wrong query position for soft-clipped reads.
- Soft-clipped reads and arrowhead placement were corrected (issues #9, #10).

# ggmethylation 0.1.0

Back-filled: initial release, never recorded in NEWS.

- `read_methylation()` parses MM/ML tags from a modBAM over a region, with
  MAPQ, strand, read-length and downsampling filters, and optional grouping by
  BAM tag or SNV genotype.
- `plot_methylation()` renders a `patchwork` composite of a packed read panel
  and a loess-smoothed per-group summary.
- `pack_reads()` assigns reads to display lanes by greedy interval scheduling.
- `smooth_methylation()` fits the per-group loess summary.
