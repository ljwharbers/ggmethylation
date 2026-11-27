# ggmethylation 0.6.0

A maintenance release: hard-clipped reads no longer get misplaced calls,
every exported function checks its arguments and names the offending one,
`plot_insertion_locus()` draws calls the way `plot_methylation()` does, and
the VCF overlay is fast enough for a few hundred variants.

## Breaking changes

- A hard-clipped alignment (CIGAR `H`) without an `MN` tag keeps no
  modification calls and is counted in the malformed-tags warning, which now
  says how many of the skipped reads were hard-clipped without `MN`. Without
  `MN` there is no way to tell whether MM/ML were rewritten for the clipped
  SEQ; when they were not, calls landed on the wrong bases without a warning.
  A hard-clipped record whose `MN` equals its SEQ length is still read. This
  is stricter than htslib, which accepts such records. minimap2 hard-clips
  supplementary alignments unless run with `-Y`, so supplementary reads from
  a pipeline that does not write `MN` lose their calls.
- `read_methylation()` validates its arguments: `bam` must be one file path,
  `group_tag` `NULL` or one two-character SAM tag name, `per_group_downsample`
  and `drop_na_group` `TRUE` or `FALSE`, and `max_reads` a whole number >= 1
  or `Inf` (`2.5` was accepted). These gave R or Rsamtools errors before.
- A `region` on a contig that is not in the BAM header is an error that names
  the contig, lists the header's contigs and hints at `chr1` vs `1` (was a raw
  Rsamtools error). A region without reads still gives a warning and an empty
  object.
- `summary()` on a `methylation_data` or `multi_methylation_data` object no
  longer prints as a side effect. It returns an object of class
  `summary.methylation_data` / `summary.multi_methylation_data` that prints
  the same report, so console use looks as before, but `s = summary(md)` and
  calls inside functions are silent.
- `plot_methylation()` checks its colour and layout arguments before drawing.
  `colour_low`, `colour_high` and `colour_ambiguous` must each be a single
  valid colour (was ggplot2's "Unknown colour name"), and `strand_colours`
  must hold valid colours when `colour_strand = TRUE` (was a render-time
  error). `panel_heights` must be `NULL` or finite positive numbers
  (`c(-1, 1)` was accepted). `sort_by` must be `NULL` or a non-empty
  character vector without `NA`s.
- `write_methylation()` validates `prefix`, `reads`, `sites` and `gzip`
  (`reads = NA` failed with "missing value where TRUE/FALSE needed");
  `read_annotations()` validates `collapse_transcripts`; `read_variants()`
  checks that `vcf` is one path to an existing file (or a URL) and rejects
  `bnd_match_tol` beyond the integer range.

## New features

- `read_methylation(include_duplicates = TRUE)` and `include_qc_fail = TRUE`
  read alignments flagged as duplicates (0x400) or QC-fail (0x200), which are
  skipped by default. Secondary alignments are never read.
- `plot_insertion_locus()` gains `call_threshold`, `call_ambiguous`,
  `colour_ambiguous` and `line_width`, with the same meaning and checks as in
  `plot_methylation()`. With `call_threshold` set, calls are drawn as
  methylated / unmethylated / ambiguous and the smoothed panel shows the
  fraction of methylated calls.
- `write_methylation()` accepts a `multi_methylation_data` object (from
  `merge_methylation()`) and writes one set of files per sample, named
  `<prefix>_<sample>_reads.tsv` etc. Sample names containing a path separator
  are an error.
- The help pages have examples that run on the shipped synthetic modBAM and
  VCF instead of `\dontrun{}` blocks with made-up file names.

## Behaviour changes

- `plot_insertion_locus()` draws modification calls the same way as
  `plot_methylation()`. With several modification codes each read bar is
  split into one strip per code, named on the y axis, so calls of different
  codes at the same base no longer hide each other. The read-panel legend
  title is now "Modification probability", and each code gets the same
  linetype in the smoothed panel as in `plot_methylation()` (the third and
  fourth codes were swapped).
- `print()` and `summary()` of a `methylation_data` object report the length
  of each read's full alignment (`bam_end - bam_pos + 1`), labelled "aligned
  length". They used the span clipped to the queried region, which
  understated every read longer than the window. An object without reads
  says so instead of printing `min=NA, median=NA, max=NA`.
- The "No reads found" warning of `read_methylation()` names only the kinds
  of alignment that were skipped.
- The help pages of `plot_methylation()` and `plot_insertion_locus()` show
  the default colours (`"#BDBDBD"`, `"#C62828"`, ...) instead of the names of
  internal constants.

## Bug fixes

- `read_annotations()` on a GTF without `transcript` rows (e.g. a UCSC Table
  Browser export) gave `NA` gene names, and with `collapse_transcripts = TRUE`
  an empty gene track. Gene names are now taken from the exon rows; a
  transcript without a gene name falls back to its transcript name and is
  kept when collapsing.
- `read_annotations()` skips the AnnotationDbi gene-name lookup when
  AnnotationDbi (a suggested package) is not installed.
- A VCF breakend without a parseable mate (a single breakend, `<BND>`) is
  labelled "⋈" instead of "⋈ NA:NA".
- The SNV overlay skips reads whose stored sequence or CIGAR is `NA` (was an
  error).

## Performance

- The VCF overlay maps all variants of a read in one step: 200 reads x 300
  SNVs take 0.2 s instead of 13 s, and the SV overlay pairs reads and SVs
  without a data frame per pair.
- Read bars are split at deletions and turned into polygons for all reads at
  once (about 25% faster plot build at ~1400 reads; output unchanged).

## Internal

- One set of argument validators (`R/checks.R`) replaces three local copies;
  error messages name the argument in backticks.
- `R/plot_methylation.R` is split: the smooth and delta panels are in
  `R/smooth_panel.R`, the consensus-deletion breaks in `R/deletion_breaks.R`.
- Removed the unused internal `seq_to_ref()` and `complement_base()`, and the
  tracked scratch files `test2.R`, `showcase.R` and `docs/superpowers/`.
  `tools/diagnose_site_coverage.R` defaults to the shipped fixture.
- Tests cover `panel_heights`, `variants`, the colour arguments and custom
  `sort_by` in `plot_methylation()`, `strand_filter` and
  `per_group_downsample` in `read_methylation()`, `collapse_transcripts`, and
  multi-sample `summary()` and export.
- CI runs R CMD check on Linux (release and oldrel-1) and macOS.

## Known limitations

- `$sites` does not record which strand a modification call is on
  (`C+m` vs `G-m`).
- Secondary alignments are never read.
- An MM entry without a `?` or `.` flag is treated like `?` (unlisted bases
  are not reported), deliberately deviating from the SAM specification,
  because older basecallers that write no flag list only CpG sites.

# ggmethylation 0.5.0

A correctness release from a full audit: MM/ML parsing follows the SAM
specification, insertion loci are clustered and their non-carriers defined
as documented, group colours are consistent across panels, VCF breakends are
matched against the right coordinate, and the docs describe what the code
does.

## Breaking changes

- `read_methylation()` skips secondary (FLAG 0x100), QC-fail (0x200) and
  duplicate (0x400) alignments. Supplementary alignments (0x800) are still
  read. There is no argument yet to include the skipped ones.
- MM entries with a `-` strand field (e.g. `G-m`) are placed on the correct
  base: bases are counted on the original read, and only the alignment
  strand decides which base is searched and in which direction.
- A read whose MM/ML tags are malformed or do not fit its sequence keeps no
  modification calls, and one warning gives the number of such reads. This
  covers a wrong number of ML values, a missing or non-numeric skip count,
  an unparseable entry, skip counts running past SEQ, an empty SEQ, and an
  `MN` tag that differs from the SEQ length. Previously these gave `NA`
  probabilities, calls at wrong positions, or an error.
- When `group_tag` is found in no read, the result is ungrouped
  (`group_tag = NULL`, no `group` column) with a warning.
- `write_methylation(format = "bed")` writes the integer score
  `round(mod_prob * 1000)`, `<read_name>:<mod_code>` in the name column, and
  rows sorted by chromosome and start.
- `read_methylation(snv = )`: a position outside `region`, alleles other
  than A/C/G/T, and `ref == alt` are errors.
- `list_insertion_loci()`: a non-carrier must span the locus on both sides
  (start <= `anchor_pos - tol_pos` and end >= `anchor_pos + tol_pos`);
  previously any overlap counted, so reads clipped at the breakpoint were
  non-carriers. The full alignment end (`bam_end`) is used, not the end
  clipped to the region, so loci near the region edge get their
  non-carriers. Reads with an insertion of a different length at the locus
  are neither carriers nor non-carriers.
- `list_insertion_loci()` and `plot_insertion_locus()` reject `NA`,
  non-numeric, negative or (where a count is expected) non-whole `tol_pos`,
  `tol_len`, `min_reads`, `min_length` and `flank`, and non-logical
  `show_smoothed` / `include_noncarriers`.
- `merge_methylation()` errors on duplicated or `NA` sample names, on no
  samples, and on `.list` being a single object instead of a list.
- `plot_methylation()` validates `call_threshold` (including `NA`),
  `call_ambiguous`, `smooth_span`, `group_colours`, `line_width`,
  `min_indel_size` and the `show_*` flags. An unnamed `group_colours` with
  fewer colours than groups is an error.
- `read_variants()` warns and returns an empty object when the region's
  chromosome is not among the VCF's contigs (was a cryptic `readVcf()`
  error). Equal-length multi-base variants are typed `"MNV"` (were
  `"deletion"`), and records whose first ALT is `*` are dropped (were
  `"SNV"`). `bnd_match_tol` must be a single finite non-negative number.

## Behaviour changes

- Groups that `group_colours` does not name (SNV `"REF"`/`"ALT"`, a third
  haplotype) get further colour-blind-safe colours instead of grey, with a
  message when you supplied the palette. An unnamed `group_colours` is
  matched to the sorted groups by position. The read bars, smooth lines and
  delta area use the same colours.
- Legends shared by several panels are collected and shown once; the smooth
  panel no longer repeats the group legend; supplementary-alignment (SA)
  partner colours are identical across samples.
- With several modification codes, each read bar is split into one strip
  per code, multi-sample smooth panels draw one line per code, and the delta
  panel is skipped with a message.
- The delta area is filled by sign as ribbons that meet at the zero
  crossings.
- Smooth lines (and the delta area and CI ribbon) break at a deletion when at
  least 75% of the group's reads spanning the whole deletion carry it
  (previously: of all reads in the window).
- Loess fits are clamped to `[0, 1]`.
- With a delta panel, the default `panel_heights` are `1, 0.35, 0.3` (were
  `1, 0.25, 0.2`), so the smooth and delta y-axis titles no longer overlap.
- `show_delta = TRUE` on ungrouped data gives a message.
- `plot_methylation()` warns when `variants` or `annotations` do not overlap
  the plotted region, or are on another chromosome.
- `plot_insertion_locus()` has a labelled x axis: reference coordinates in
  the flanks and `+1...+W bp` in the insertion. Non-carriers are drawn below
  the carriers as documented, with a centred separator and no lane numbers.
- The gene track uses `theme_ggmethylation()` like the other panels.

## Bug fixes

### Parsing and reading

- Combined-code MM entries (`C+mh`) are supported, with interleaved ML
  values and correct ML offsets for later entries; numeric ChEBI codes count
  as one code.
- Every MM entry carrying the requested code is used, not only the last.
- Canonical base `N` matches any base, also on reverse-strand reads.
- `summary()` works on an empty object and when every group is `NA`.
- An `snv` request with no reads left still reports group tag `"SNV"`.
- `parse_region()` handles contig names containing colons and rejects
  reversed, out-of-range and non-scalar regions with a clear message.
- `max_reads`, `min_mapq` and `min_read_length` are validated.
- `.csi` BAM indexes are accepted.

### Insertion loci

- Insertions are clustered by position first and then by length, so a
  single off-length insertion no longer splits a locus or turns its carriers
  into non-carriers.
- `mean_ins_mod_prob` uses only each carrier's own insertion.
- `plot_insertion_locus()`: the reference base right after the insertion is
  drawn in the right flank instead of on top of the first inserted base, and
  both flanks show exactly `flank` bases. Inserted bases beyond the locus
  median length are clipped at the insertion boundary. The read and smooth
  panels share one x range, and each modification code gets its own smooth
  line (line type) instead of being averaged.

### Main plot

- SNV and three-haplotype grouping are no longer drawn in grey with
  "No shared levels" warnings.
- An unnamed `group_colours` now reaches the delta panel.
- The delta area and the CI ribbon are interrupted at consensus deletions.
- Binary mode no longer loses smooth-line segments when the fit leaves
  `[0, 1]`.
- Multi-sample smooth panels no longer pool 5mC and 5hmC.
- A read with an SA tag and a large deletion no longer gets a second SA
  marker.
- Sites next to a clipped read end are no longer dropped.
- A group with a single position no longer crashes the delta track; `NA`
  probabilities are ignored when smoothing.
- Two ungrouped samples in a mixed merge no longer share one CI ribbon.
- More than four modification codes no longer give `NA` line
  types (up to 13 are supported, with a clear error beyond).

### Annotations and variants

- Reads' SA positions are matched against the BND mate position. Previously
  they were compared with the BND's own position, so correct breakpoints
  were never validated.
- The gene track's x axis lines up with the other panels (it was stretched
  by about 10%); strand arrows at the region edge are clamped instead of
  dropped.
- `SVTYPE`/`END`/`SVLEN` declared `Number=.` are handled; VCFs without
  `##contig` lines can be read; a declared contig with no records returns an
  empty object.
- `read_annotations()` warns clearly when the region's chromosome is not in
  the annotation, and uses `txdbmaker` to build TxDbs when it is installed.

## New

- `read_methylation(seed = )` for reproducible downsampling; the caller's
  random number stream is restored afterwards.
- `$reads$bam_end`: the alignment end, next to `bam_pos`, not clipped to the
  region.
- The `MN` tag is read and checked against the SEQ length.
- `read_methylation()` warns when no read carries MM/ML tags (mentioning
  legacy `Mm`/`Ml` tags when those are present) and when a requested
  `mod_code` is in no MM tag.
- `list_insertion_loci(min_length = )` ignores insertions shorter than
  `min_length` (default 1 keeps all; e.g. 50 for real ONT/PacBio data).
- `read_variants()` reads only the INFO fields it uses and no genotypes.
- The synthetic modBAM/VCF fixture is shipped in `inst/extdata/`, and the
  vignette runs on it instead of a cache built from a non-public BAM.

## Internal

- MM/ML tags are parsed once per read for all requested codes (new core
  `.parse_mods()`; `parse_mm_ml()` is a thin wrapper).
- One resolver, `.resolve_group_colours()`, provides the group colours for
  every panel.
- DESCRIPTION: `Depends: R (>= 4.0)`, `URL`/`BugReports`, `grDevices` and
  `tools` in Imports, the optional Bioconductor packages in Suggests; the
  licence file is the standard MIT stub (full text in `LICENSE.md`).
- GitHub Actions R-CMD-check workflow; shared test helper; integration tests
  can be pointed at another BAM with `GGMETHYLATION_TEST_BAM` /
  `GGMETHYLATION_TEST_REGION`.

## Known limitations

- Hard-clipped records without an `MN` tag whose skip counts still fit the
  clipped sequence are placed without a warning; this cannot be detected.
- `$sites` does not record which strand a modification call is on
  (`C+m` vs `G-m`).
- Secondary and duplicate alignments cannot be re-included yet.
- An MM entry without a `?` or `.` flag is treated like `?` (unlisted bases
  are not reported), deliberately deviating from the SAM specification,
  because older basecallers that write no flag list only CpG sites.

# ggmethylation 0.4.0

A simplification release: fewer arguments, one code path per concept. Plots
for equivalent calls are unchanged except where noted under "Behaviour
changes".

## Breaking changes

| 0.3.0 | 0.4.0 |
|---|---|
| `plot_methylation(call_mode = "binary", call_threshold = 0.5)` | `plot_methylation(call_threshold = 0.5)` |
| `plot_methylation(call_mode = "continuous")` (default) | `plot_methylation()` (`call_threshold = NULL`) |
| `read_methylation(snv_position = 1500, ref_base = "C", alt_base = "T")` | `read_methylation(snv = list(position = 1500, ref = "C", alt = "T"))` |
| `plot_methylation(variants = vd, bnd_match_tol = 100)` | `vd = read_variants(vcf, region, bnd_match_tol = 100)`; `plot_methylation(variants = vd)` |
| `plot_insertion_locus(m, loci$locus_id[1], tol_pos = 10, tol_len = 0.2)` | `plot_insertion_locus(data, loci[1, ])` |
| `insertion_sites(m)`, `list_insertion_loci(m, ...)` | `insertion_sites(data)`, `list_insertion_loci(data, ...)` |
| `plot_methylation(mod_code_shapes = ...)` | removed (it was never drawn) |
| `write_methylation(...)` | the unused `...` is removed |

- `call_ambiguous` without `call_threshold` is now an error instead of being
  silently ignored.
- `list_insertion_loci()` gains list-columns `carriers` (insertion length per
  carrier read, named by `read_name`) and `noncarriers`, plus the `tol_pos` /
  `tol_len` used. `plot_insertion_locus()` uses them directly instead of
  re-clustering, so its tolerances can no longer drift from the table's.
- `variant_data` objects carry a `bnd_match_tol` element.

## Behaviour changes

- Multi-sample plots (`merge_methylation()`) now pack reads with the same
  clip-aware lane spacing as single-sample plots, and draw the loess CI ribbon
  (`show_ci`) on the shared smooth panel.
- Reads without a group (e.g. no HP tag) no longer get a stray `"NA"` smooth
  line in grouped multi-code panels or in multi-sample panels. They were
  already excluded from the single-code smooth panel; the line only appeared
  because group names were pasted into a string key.

## Bug fixes

- `read_methylation(drop_na_group = TRUE)` combined with `min_mapq` or
  `min_read_length` filtered the wrong reads (or failed with "subscript
  contains NAs"): the filters indexed full-length BAM vectors against the
  already-shortened reads table.
- `plot_methylation()` failed when merging a grouped with an ungrouped sample.
- An unknown `sort_by` column silently dropped every read; it is now an error.
  The vignette example that hit this (`sort_by = "position"`) is fixed.

## Internal

- `plot_methylation()` and the multi-sample renderer share one pipeline
  (`.prepare_reads()`, `.build_smooth_panel()`, ...) instead of duplicated
  branches; `build_read_panel()`, `read_methylation()` and
  `smooth_methylation()` are split into small helpers.
- Unused `GenomicAlignments` import removed.
- A synthetic modBAM/VCF fixture (`tests/testthat/fixtures/`, generated by
  `data-raw/make_test_fixture.R`) gives `read_methylation()` and
  `read_variants()` unit tests that run without the staging BAM.

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
