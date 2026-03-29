# Plan D12: Richer `print`/`summary` for `methylation_data`

**Goal:** Improve `print.methylation_data` with per-strand and per-group breakdowns, and add a `summary.methylation_data` S3 method returning a structured list for programmatic use.

**Affected files:**
- `R/read_methylation.R` — modify `print.methylation_data`, add `summary.methylation_data`
- `NAMESPACE` — add `S3method(summary,methylation_data)` via `devtools::document()`

---

## Current `print` Output (for reference)

```
methylation_data object
Region: chr1:1000-2000
Reads: 42
Modification sites: 387
Modification code: m
Group tag: HP (2 groups)
```

---

## Improved `print.methylation_data`

### Revised output (grouped)

```
methylation_data object
Region: chr1:1000-2000
Reads: 42  (+ strand: 21,  - strand: 21)
Median read length: 4821 bp
Modification sites: 387
Modification code: m
Group tag: HP (2 groups)
  H1: 21 reads, mean methylation 0.72
  H2: 21 reads, mean methylation 0.41
```

### Revised output (ungrouped)

```
methylation_data object
Region: chr1:1000-2000
Reads: 42  (+ strand: 21,  - strand: 21)
Median read length: 4821 bp
Modification sites: 387
Modification code: m
```

### Implementation steps

#### Step 1 — Compute strand counts

```r
n_plus  <- sum(x$reads$strand == "+", na.rm = TRUE)
n_minus <- sum(x$reads$strand == "-", na.rm = TRUE)
```

#### Step 2 — Compute median read length

```r
read_lengths <- x$reads$end - x$reads$start + 1L
med_len <- if (length(read_lengths) > 0L) as.integer(median(read_lengths)) else NA_integer_
```

#### Step 3 — Replace "Reads:" line

```r
cat(sprintf("Reads: %d  (+ strand: %d,  - strand: %d)\n",
            nrow(x$reads), n_plus, n_minus))
```

#### Step 4 — Add median read length line

```r
if (!is.na(med_len))
  cat(sprintf("Median read length: %d bp\n", med_len))
```

#### Step 5 — Add per-group lines inside the `!is.null(x$group_tag)` block

```r
groups <- sort(unique(x$reads$group[!is.na(x$reads$group)]))
for (g in groups) {
  n_g   <- sum(x$reads$group == g, na.rm = TRUE)
  mods  <- x$sites$mod_prob[!is.na(x$sites$group) & x$sites$group == g]
  mean_g <- if (length(mods) > 0L) round(mean(mods, na.rm = TRUE), 2L) else NA_real_
  cat(sprintf("  %s: %d reads, mean methylation %.2f\n", g, n_g, mean_g))
}
```

#### Step 6 — Keep `invisible(x)` as return value

---

## New `summary.methylation_data` Method

### Design

`summary.methylation_data` **prints** a compact table and **returns** a named list invisibly — same pattern as `summary.lm()`. This supports both interactive inspection and programmatic access.

### Function signature

```r
#' @export
summary.methylation_data <- function(object, ...) { ... }
```

### Returned list structure

```r
list(
  region                = character(1),   # "chr1:1000-2000"
  n_reads               = integer(1),
  n_sites               = integer(1),
  mod_code              = character(1),
  group_tag             = character(1) or NULL,
  strand                = data.frame(strand = character, n_reads = integer),
  read_length           = list(median = integer, min = integer, max = integer),
  groups                = data.frame(          # NULL when ungrouped
    group            = character(),
    n_reads          = integer(),
    mean_mod_prob    = numeric(),
    median_mod_prob  = numeric()
  ),
  overall_mean_mod_prob = numeric(1)
)
```

### Printed output (grouped example)

```
methylation_data summary
========================
Region:            chr1:1000-2000
Reads:             42
  + strand:        21
  - strand:        21
  Read length:     min=800, median=4821, max=12043 bp
Sites:             387
Modification:      m
Overall mean mod:  0.58

Group breakdown (HP):
  group  n_reads  mean_mod  median_mod
  H1          21      0.72        0.75
  H2          21      0.41        0.38
```

### Implementation steps

#### Step 1 — Extract region string

```r
chrom <- as.character(GenomicRanges::seqnames(object$region))
reg_str <- sprintf("%s:%d-%d", chrom,
                   GenomicRanges::start(object$region),
                   GenomicRanges::end(object$region))
```

#### Step 2 — Strand table

```r
strand_df <- as.data.frame(table(strand = object$reads$strand),
                            stringsAsFactors = FALSE)
names(strand_df)[2] <- "n_reads"
```

#### Step 3 — Read length stats

```r
lens <- object$reads$end - object$reads$start + 1L
rl <- if (length(lens) > 0L)
  list(median = as.integer(median(lens)), min = min(lens), max = max(lens))
else
  list(median = NA_integer_, min = NA_integer_, max = NA_integer_)
```

#### Step 4 — Overall mean mod prob

```r
overall_mean <- if (nrow(object$sites) > 0L)
  round(mean(object$sites$mod_prob, na.rm = TRUE), 4L)
else NA_real_
```

#### Step 5 — Per-group table

```r
groups_df <- NULL
if (!is.null(object$group_tag)) {
  g_levels <- sort(unique(object$reads$group[!is.na(object$reads$group)]))
  groups_df <- do.call(rbind, lapply(g_levels, function(g) {
    n_g  <- sum(object$reads$group == g, na.rm = TRUE)
    mods <- object$sites$mod_prob[
      !is.na(object$sites$group) & object$sites$group == g]
    data.frame(
      group           = g,
      n_reads         = n_g,
      mean_mod_prob   = if (length(mods) > 0L) round(mean(mods,   na.rm=TRUE), 4L) else NA_real_,
      median_mod_prob = if (length(mods) > 0L) round(median(mods, na.rm=TRUE), 4L) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}
```

#### Step 6 — Build return list and print

```r
out <- list(region=reg_str, n_reads=nrow(object$reads), n_sites=nrow(object$sites),
            mod_code=object$mod_code, group_tag=object$group_tag,
            strand=strand_df, read_length=rl, groups=groups_df,
            overall_mean_mod_prob=overall_mean)

cat("methylation_data summary\n========================\n")
cat(sprintf("Region:            %s\n", reg_str))
cat(sprintf("Reads:             %d\n", out$n_reads))
for (i in seq_len(nrow(strand_df)))
  cat(sprintf("  %s strand:        %d\n", strand_df$strand[i], strand_df$n_reads[i]))
cat(sprintf("  Read length:     min=%d, median=%d, max=%d bp\n",
            rl$min, rl$median, rl$max))
cat(sprintf("Sites:             %d\n", out$n_sites))
cat(sprintf("Modification:      %s\n", out$mod_code))
cat(sprintf("Overall mean mod:  %.2f\n", overall_mean))
if (!is.null(groups_df)) {
  cat(sprintf("\nGroup breakdown (%s):\n", object$group_tag))
  print(groups_df, row.names = FALSE, digits = 4)
}

invisible(out)
```

#### Step 7 — Add roxygen2 block

`@param object`, `@param ...`, `@return` describing the list structure, `@examples`, `@export`.

#### Step 8 — Run `devtools::document()`

Confirm `S3method(summary,methylation_data)` in `NAMESPACE`.

---

## Edge cases

| Scenario | Handling |
|---|---|
| No reads | Strand counts 0; `read_length` all `NA`; guard `length(lens) > 0L`. |
| No sites | `overall_mean_mod_prob = NA`; per-group means `NA`. |
| Ungrouped | `groups` is `NULL`; group section skipped in output. |
| Reads with `group = NA` | Filter with `!is.na()` before `unique()` and group-level stats. |
| SNV-grouped data | Same code path; `"SNV"` appears as group tag label. |
| Strand values other than `+`/`-` (e.g. `*`) | `table()` includes them as a row; shown in output. |

---

## Testing approach

Create `tests/testthat/test-print-summary.R`.

**`print.methylation_data` tests:**
1. Capture output with `capture.output(print(mock_md))`. Check it contains `"+ strand:"`, `"- strand:"`, `"Median read length:"`.
2. Grouped mock: output contains per-group lines with `"mean methylation"`.
3. Ungrouped mock: no group lines appear.
4. Empty data: no error, output contains `"Reads: 0"`.

**`summary.methylation_data` tests:**
1. Return value is a `list` with all expected names.
2. `$strand` is a `data.frame` with correct counts.
3. `$read_length$median` is the correct integer.
4. `$overall_mean_mod_prob` is within numeric tolerance of expected value.
5. Grouped mock: `$groups` is a non-`NULL` `data.frame` with correct columns.
6. Ungrouped mock: `$groups` is `NULL`.
7. Empty data: no error; `$n_reads == 0`, `$overall_mean_mod_prob` is `NA`.
8. Return is invisible: `x <- summary(mock_md)` assigns without double-printing.
