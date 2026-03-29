# Plan B8: Multiple Modification Codes Simultaneously

## Goal

Allow `read_methylation()` to accept a character vector of modification codes
(e.g., `c("m", "h")`) and render all requested modifications on the same plot,
distinguishing them by dot shape.

---

## Affected Files

| File | Nature of change |
|---|---|
| `R/read_methylation.R` | Accept vector `mod_code`; loop `parse_mm_ml()` per code per read; combine results; update `$mod_code` slot; update `print.methylation_data()`; update `empty_methylation_data()` |
| `R/plot_methylation.R` | Map `mod_code` to shape aesthetic; add `mod_code_shapes` parameter; adapt smooth panel to produce one line per group-×-mod_code combination; update legend layout |
| `R/smooth_methylation.R` | Accept an optional additional grouping column (`mod_code`) so it can smooth per group-×-mod_code combination |

`R/parse_mm_ml.R` requires no changes — it already handles one code per call
and already returns a two-column data.frame (`position`, `mod_prob`) that the
caller labels with `mod_code`.

---

## Design Decision: Shape Aesthetic for Mod Code

**Chosen approach:** shared colour scale for `mod_prob`; `mod_code` mapped to
point `shape`.

**Rejected alternatives:**

- *Separate colour scales via `ggnewscale`* — already used for the group
  colour / mod_prob split. Adding a third scale for a second mod code would
  require two `new_scale_colour()` calls, producing an unwieldy legend and
  fragile layer ordering.
- *Faceted panels per mod code* — loses the key visual comparison: seeing 5mC
  and 5hmC at the same read and position simultaneously. Facets also double
  the vertical space.
- *Different colour palettes per mod code* — visually ambiguous when a
  position has both codes at different probabilities; colour already encodes
  probability.

**Rationale for shapes:** Shape is a perceptually distinct, low-noise channel
that does not conflict with the existing colour-for-probability encoding. Two
shapes (e.g., filled circle vs. filled square) are immediately distinguishable
and integrate naturally into a single ggplot2 `scale_shape_manual()` legend.
The approach adds zero new dependencies.

---

## New Parameters / Signatures

### `read_methylation()`

```r
read_methylation(bam, region,
                 mod_code = "m",        # now accepts character vector
                 group_tag = NULL,
                 snv_position = NULL,
                 ref_base = NULL,
                 alt_base = NULL,
                 max_reads = 200L)
```

`mod_code` changes from "a single character string" to "a character vector of
one or more modification codes". The default `"m"` is unchanged.

`$mod_code` in the returned object stores the full vector as supplied.

### `plot_methylation()`

```r
plot_methylation(data, sort_by = NULL,
                 colour_low = "#BDBDBD",
                 colour_high = "#C62828",
                 dot_size = 1.5,
                 group_colours = NULL,
                 smooth_span = 0.3,
                 panel_heights = c(3, 1),
                 mod_code_shapes = NULL)   # NEW
```

`mod_code_shapes`: named integer vector mapping each mod code to a ggplot2
shape number, e.g. `c(m = 16, h = 15)`. When `NULL` (default), shapes are
assigned automatically from a built-in default sequence: `c(16, 15, 17, 18)`
(filled circle, filled square, filled triangle, filled diamond).

---

## Implementation Steps

### Step 1 — `read_methylation()`: normalise and validate `mod_code`

In the input validation block (section 1), add:

```r
mod_code <- as.character(mod_code)          # coerce in case of factor
if (length(mod_code) == 0L || any(!nzchar(mod_code))) {
  stop("'mod_code' must be a non-empty character vector of modification codes.",
       call. = FALSE)
}
mod_code <- unique(mod_code)                # drop accidental duplicates
```

No other validation is needed at this stage; an unrecognised code simply
returns zero rows from `parse_mm_ml()`, which is handled in step 3.

### Step 2 — `read_methylation()`: inner per-read loop (section 7)

Replace the current single `parse_mm_ml()` call inside the `for (j in ...)` loop
with a nested loop over `mod_code`:

```r
sites_list <- vector("list", length(keep))

for (j in seq_along(keep)) {
  idx       <- keep[j]
  seq_str   <- as.character(bam_data$seq[[idx]])
  mm        <- bam_data$tag$MM[idx]
  ml        <- bam_data$tag$ML[[idx]]

  code_results <- vector("list", length(mod_code))
  for (ci in seq_along(mod_code)) {
    result <- parse_mm_ml(
      seq    = seq_str,
      mm_tag = mm,
      ml_tag = ml,
      mod_code = mod_code[ci],
      strand = reads$strand[j],
      cigar  = bam_data$cigar[idx],
      pos    = bam_data$pos[idx]
    )
    if (nrow(result) > 0L) {
      result$read_name <- reads$read_name[j]
      result$mod_code  <- mod_code[ci]      # label before combining
    }
    code_results[[ci]] <- result
  }
  sites_list[[j]] <- do.call(rbind, code_results)
}
```

This ensures `mod_code` is stamped on rows produced by `parse_mm_ml()` before
any `rbind`.

### Step 3 — `read_methylation()`: remove the post-hoc `mod_code` stamp (section 9)

Delete the current line:

```r
sites$mod_code <- mod_code
```

The column is now populated correctly per-code inside the loop (Step 2).
When `sites` has zero rows the empty data.frame construction (section 7
fallback) must already include the column:

```r
sites <- data.frame(
  position  = integer(0L),
  mod_prob  = numeric(0L),
  read_name = character(0L),
  mod_code  = character(0L),
  stringsAsFactors = FALSE
)
```

This matches the existing `empty_methylation_data()` schema — no change needed
there beyond ensuring `mod_code` is stored as a vector in the object slot.

### Step 4 — `read_methylation()`: store `mod_code` vector in returned object

The `structure()` call at the end already stores `mod_code = mod_code`. When
`mod_code` is now a vector, this naturally stores the full vector. No code
change needed here.

### Step 5 — `print.methylation_data()`: update for vector `mod_code`

Replace:

```r
cat(sprintf("Modification code: %s\n", x$mod_code))
```

with:

```r
cat(sprintf("Modification code(s): %s\n", paste(x$mod_code, collapse = ", ")))
```

### Step 6 — `smooth_methylation()`: add optional `mod_code_col` argument

The smooth function currently groups solely by `group_col`. When multiple mod
codes are present, we need separate lines per mod-code (within each group).
Add an optional `mod_code_col` argument:

```r
smooth_methylation <- function(sites, group_col = "group",
                               mod_code_col = NULL, span = 0.3)
```

When `mod_code_col` is not NULL, create a composite grouping key before the
existing loop:

```r
if (!is.null(mod_code_col) && mod_code_col %in% names(sites)) {
  sites$.smooth_group <- paste(sites[[group_col]], sites[[mod_code_col]], sep = ":::")
  effective_group_col <- ".smooth_group"
} else {
  effective_group_col <- group_col
}
```

Use `effective_group_col` throughout the internal loop. After `rbind`, split
the composite key back into the original columns:

```r
if (!is.null(mod_code_col) && mod_code_col %in% names(sites)) {
  parts <- strsplit(out[[".smooth_group"]], ":::", fixed = TRUE)
  out[[group_col]]    <- vapply(parts, `[[`, character(1L), 1L)
  out[[mod_code_col]] <- vapply(parts, `[[`, character(1L), 2L)
  out[[".smooth_group"]] <- NULL
}
```

Return columns `c("position", "mean_prob", group_col)` when single-code,
`c("position", "mean_prob", group_col, mod_code_col)` when multi-code.

### Step 7 — `plot_methylation()`: resolve `mod_code_shapes`

Early in the function, after input validation, compute the shape mapping:

```r
codes <- unique(data$sites$mod_code)
default_shapes <- c(16L, 15L, 17L, 18L)   # circle, square, triangle, diamond

if (is.null(mod_code_shapes)) {
  mod_code_shapes <- stats::setNames(
    default_shapes[seq_along(codes)],
    codes
  )
} else {
  # Validate that all codes present in data have an entry
  missing_codes <- setdiff(codes, names(mod_code_shapes))
  if (length(missing_codes) > 0L) {
    warning(
      "mod_code_shapes has no entry for code(s): ",
      paste(missing_codes, collapse = ", "),
      ". Using default shapes.", call. = FALSE
    )
    extras <- stats::setNames(
      default_shapes[seq_along(missing_codes)],
      missing_codes
    )
    mod_code_shapes <- c(mod_code_shapes, extras)
  }
}

multi_code <- length(codes) > 1L
```

### Step 8 — `plot_methylation()`: top panel dot layer

In both the grouped and ungrouped branches, change the `geom_point()` call
that renders modification dots. Replace:

```r
ggplot2::aes(x = .data$position, y = .data$lane, colour = .data$mod_prob)
```

with:

```r
ggplot2::aes(
  x      = .data$position,
  y      = .data$lane,
  colour = .data$mod_prob,
  shape  = if (multi_code) .data$mod_code else NULL
)
```

Add after the existing `scale_colour_gradient()` call (in both branches):

```r
if (multi_code) {
  p_top <- p_top +
    ggplot2::scale_shape_manual(
      values = mod_code_shapes,
      name   = "Modification"
    )
}
```

When `multi_code` is FALSE the shape aesthetic is NULL, so `scale_shape_manual`
is not added and the single-code path is visually identical to today.

### Step 9 — `plot_methylation()`: smooth panel for multiple codes

Replace the current `smooth_methylation()` call(s) with multi-code-aware
calls.

**Ungrouped, single code (current behaviour preserved):**

```r
sites_smooth$group <- "all"
smoothed <- smooth_methylation(sites_smooth, group_col = "group",
                               span = smooth_span)
# plot as before: one line, no colour
```

**Ungrouped, multi-code:**

```r
sites_smooth$group <- "all"
smoothed <- smooth_methylation(sites_smooth, group_col = "group",
                               mod_code_col = "mod_code", span = smooth_span)
# aes: x = position, y = mean_prob, linetype/colour = mod_code
```

Use `linetype` (or `colour`) mapped to `mod_code` for the smooth lines.
`colour` is recommended here because in the ungrouped case the read bars are
a fixed grey and the colour scale for mod_prob only applies to the top panel
dots (not the smooth panel). Using colour for mod_code in the smooth panel
is unambiguous.

**Grouped, single code (current behaviour preserved):**

```r
smoothed <- smooth_methylation(data$sites, group_col = "group",
                               span = smooth_span)
# plot as before: colour = group
```

**Grouped, multi-code:**

```r
smoothed <- smooth_methylation(data$sites, group_col = "group",
                               mod_code_col = "mod_code", span = smooth_span)
# aes: x = position, y = mean_prob, colour = group, linetype = mod_code
```

Map `group` to colour (as today) and `mod_code` to `linetype` in the smooth
panel. Add `scale_linetype_manual()` with sensible defaults (`"solid"`,
`"dashed"`, `"dotdash"`, `"dotted"`).

### Step 10 — `plot_methylation()`: mean_mod_prob per-read (step 3 in plot function)

The current `tapply` computes the mean of `mod_prob` across all sites for a
read, irrespective of code. This is used only for sorting and is acceptable
averaged across all codes. No change needed.

---

## Backwards Compatibility

- `mod_code = "m"` (a length-1 character) passes `as.character()` unchanged,
  `unique()` returns `"m"`, `length(mod_code) == 1` throughout. `multi_code`
  is `FALSE`. All shape and linetype layers are suppressed. The output
  `$mod_code` slot is `"m"` (a length-1 character), matching today's value.
  Serialised objects from before this change will continue to print and plot
  correctly.
- `empty_methylation_data()` schema unchanged: `mod_code` column already
  exists in its `sites` data.frame.
- All existing function signatures gain only one new optional argument
  (`mod_code_shapes` / `mod_code_col`) with a `NULL` default, so existing
  call sites continue to work.

---

## Edge Cases

| Situation | Handling |
|---|---|
| A requested code is absent from all reads | `parse_mm_ml()` returns zero rows; after `rbind` those reads contribute no rows for that code. The `sites` data.frame has no rows with that `mod_code` value. `plot_methylation()` calls `unique(data$sites$mod_code)` to detect what is actually present — absent codes never appear in the shape/linetype scale. |
| A code is absent from some reads but present in others | Those reads simply have no dots for that code. The shape scale shows only codes actually observed (`unique(data$sites$mod_code)`). |
| Both codes hit the same genomic position in the same read | Two rows at the same `(position, read_name)` with different `mod_code` values — this is valid and expected (e.g., cytosine can carry both 5mC and 5hmC scores). They render as two overlapping dots with different shapes, distinguishable in the legend. |
| `mod_code` is length > 4 (more than default shapes) | `default_shapes[seq_along(codes)]` will index out of bounds. Add a guard that recycles or errors: `if (length(codes) > length(default_shapes)) stop("More than 4 mod codes requested; supply mod_code_shapes explicitly.", call. = FALSE)`. |
| `mod_code` contains duplicates | `unique(mod_code)` in the normalisation step deduplicates silently. |
| MM tag has multi-character code (e.g., ChEBI "21839") | `parse_mm_ml()` already extracts `code <- substring(spec, 3)` which handles arbitrary-length codes. No change needed. |
| `smooth_methylation()` called with `mod_code_col = "mod_code"` but the column is absent (legacy `sites` object) | Guard with `mod_code_col %in% names(sites)` before using it; fall back to single-group smoothing with a warning. |
| Region has reads but none carry any requested mod code | `sites` has zero rows. The existing zero-row handling in `plot_methylation()` already renders a "No reads" placeholder — extend this check to `nrow(data$sites) == 0L` (already done) or add a separate "no modification sites" message. |

---

## Testing Approach

1. **Single-code regression** — run `read_methylation(bam, region, mod_code = "m")`
   on the existing PacBio HiFi test BAM (see `reference_test_bam.md`). Confirm
   output is structurally identical to pre-B8 output: `$mod_code` is `"m"`,
   `$sites$mod_code` column is all `"m"`, plot is visually identical to
   current baseline.

2. **Two-code data parsing** — run `read_methylation(bam, region, mod_code = c("m", "h"))`.
   Assert:
   - `$mod_code` is `c("m", "h")`.
   - `$sites$mod_code` contains only values in `c("m", "h")`.
   - Row counts for code `"m"` match those from the single-code call.
   - If 5hmC is present in the test BAM, code `"h"` rows also appear.

3. **Absent code** — run with `mod_code = c("m", "a")` where `"a"` (6mA) is
   absent from the test BAM. Assert no error; assert `unique(sites$mod_code)`
   does not include `"a"`; assert `print()` still works.

4. **Plot — ungrouped multi-code** — call `plot_methylation()` on the two-code
   object. Assert the returned patchwork contains a `scale_shape_manual` layer
   (inspect `p_top$scales$scales`); visually verify distinct dot shapes.

5. **Plot — grouped multi-code** — call with `group_tag = "HP"`. Assert smooth
   panel `geom_line` uses `linetype` mapped to `mod_code`.

6. **Custom shapes** — call `plot_methylation(md, mod_code_shapes = c(m = 21, h = 22))`.
   Assert the manual shape scale uses values 21 and 22.

7. **`smooth_methylation()` unit test** — construct a synthetic `sites`
   data.frame with two groups and two mod codes. Call `smooth_methylation()`
   with `mod_code_col = "mod_code"`. Assert the result has one row-group per
   unique `(group, mod_code)` combination.

8. **Edge: duplicate codes** — `read_methylation(bam, region, mod_code = c("m", "m"))`.
   After deduplication the call should behave identically to `mod_code = "m"`.
