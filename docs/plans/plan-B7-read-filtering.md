# Plan B7: Read Filtering Parameters

**Goal:** Add `min_mapq`, `strand_filter`, and `min_read_length` parameters to `read_methylation()` to filter reads before MM/ML parsing.

**Affected files:** `R/read_methylation.R` only.

---

## New parameters

Add to function signature:

```r
read_methylation <- function(
  bam, region,
  mod_code         = "m",
  group_tag        = NULL,
  snv_position     = NULL,
  ref_base         = NULL,
  alt_base         = NULL,
  max_reads        = 200L,
  min_mapq         = 0L,          # no filtering by default
  strand_filter    = c("+", "-"), # both strands by default
  min_read_length  = 0L           # no filtering by default
)
```

---

## Implementation steps

### Step 1 — Add "mapq" to `what` in `ScanBamParam`

Current:
```r
what <- c("qname", "pos", "cigar", "strand", "seq")
```
Change to:
```r
what <- c("qname", "pos", "cigar", "strand", "seq", "mapq")
```

`Rsamtools::scanBam()` returns MAPQ as `bam_data$mapq` (integer vector, same length as `bam_data$qname`).

### Step 2 — Validate new parameters (add to section 1, after existing input validation)

```r
strand_filter <- match.arg(strand_filter, choices = c("+", "-"), several.ok = TRUE)
if (!is.integer(min_mapq))       min_mapq <- as.integer(min_mapq)
if (!is.integer(min_read_length)) min_read_length <- as.integer(min_read_length)
```

### Step 3 — Insert new filtering section between group extraction and downsampling

Insert after the group column is added (section 4) and before the downsampling block (section 6). Initialise `bam_indices` before this block:

```r
bam_indices <- seq_len(nrow(reads))   # move here from current line 103

# --- 4b. Apply read-level filters ---
n_before <- nrow(reads)
filter_mask <- rep(TRUE, n_before)

# MAPQ filter
if (min_mapq > 0L) {
  mapq_vec <- bam_data$mapq
  filter_mask <- filter_mask & (!is.na(mapq_vec) & mapq_vec >= min_mapq)
}

# Strand filter
if (!identical(sort(strand_filter), c("+", "-"))) {
  filter_mask <- filter_mask & (reads$strand %in% strand_filter)
}

# Read length filter (ref_widths computed earlier from CIGAR)
if (min_read_length > 0L) {
  filter_mask <- filter_mask & (ref_widths >= min_read_length)
}

if (!all(filter_mask)) {
  n_removed <- sum(!filter_mask)
  frac_removed <- n_removed / n_before
  if (frac_removed > 0.5) {
    warning(sprintf(
      "%.0f%% of reads (%d/%d) were removed by filters (min_mapq=%d, strand_filter=c(%s), min_read_length=%d).",
      frac_removed * 100, n_removed, n_before, min_mapq,
      paste(sprintf('"%s"', strand_filter), collapse = ", "),
      min_read_length
    ), call. = FALSE)
  }
  reads       <- reads[filter_mask, , drop = FALSE]
  bam_indices <- bam_indices[filter_mask]
  rownames(reads) <- NULL
}
```

Note: `ref_widths` is computed from `bam_data$cigar` early in the function, so it is available here. Reference-space length is consistent with how `reads$start`/`end` are derived.

### Step 4 — Early return when no reads remain

```r
if (nrow(reads) == 0L) {
  warning("No reads remain after applying filters.", call. = FALSE)
  return(empty_methylation_data(gr, mod_code, group_tag))
}
```

### Step 5 — Remove the original `bam_indices <- seq_len(nrow(reads))` line

(Now initialised before the filter block in Step 3.)

### Step 6 — Roxygen documentation

```r
#' @param min_mapq Integer. Minimum mapping quality (MAPQ) threshold
#'   (default \code{0L}, no filtering). Reads with MAPQ below this value,
#'   or with missing MAPQ, are excluded before MM/ML parsing.
#' @param strand_filter Character vector. Which strands to include. One or
#'   both of \code{"+"} and \code{"-"} (default \code{c("+", "-")}, both
#'   strands). Use \code{"+"} or \code{"-"} alone to restrict to a single
#'   strand.
#' @param min_read_length Integer. Minimum read length in reference-space
#'   base pairs (default \code{0L}, no filtering). Reads shorter than this
#'   value are excluded before MM/ML parsing.
```

---

## Edge cases

| Scenario | Behaviour |
|---|---|
| All reads pass filters | No filtering, no warning. |
| All reads fail a filter | Returns `empty_methylation_data()` with a warning. |
| MAPQ is `NA` (unmapped/supplementary) | Treated as failing when `min_mapq > 0`. |
| `strand_filter = c("+", "-")` (default) | Strand filter branch skipped entirely. |
| Filtering reduces reads below `max_reads` | Downsampling section handles gracefully (no-op). |
| Combined with `group_tag` | Group column already on `reads` before filter block; row removal preserves group column correctly. |
| `min_read_length` vs soft-clipped reads | Uses `ref_widths` from `cigarWidthAlongReferenceSpace()`, ignoring soft clips — consistent with `reads$start`/`end`. |

---

## Testing approach

1. `min_mapq`: use a BAM with known MAPQ distribution; confirm low-MAPQ reads absent from `result$reads`.
2. `strand_filter = "+"`: confirm all returned reads have `reads$strand == "+"`.
3. `min_read_length = 5000L`: confirm all returned reads span >= 5000 bp (`reads$end - reads$start + 1 >= 5000`).
4. `> 50%` filtered: verify warning fires.
5. All reads filtered: verify `empty_methylation_data` structure is returned with a warning.
6. Default parameters produce identical output to current code.
