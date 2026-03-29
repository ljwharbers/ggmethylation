# Plan B6: Per-Group Downsampling

**Goal:** Add a `per_group_downsample` parameter so users can opt in to applying the `max_reads` cap independently within each group, preventing high-coverage groups from crowding out low-coverage groups.

**Affected files:** `R/read_methylation.R` only.

**New parameters:** `per_group_downsample = FALSE` — when `TRUE` and grouping is active, `max_reads` is applied per group; otherwise the existing global cap applies.

---

## Roxygen update

Add `@param per_group_downsample` documenting that when `TRUE` and grouping is active (via `group_tag` or `snv_position`), `max_reads` is applied independently per group. When `FALSE` (default), the existing global cap behaviour is unchanged.

---

## Implementation steps

### Step 1 — Add `per_group_downsample` to function signature

```r
read_methylation <- function(
  bam, region,
  mod_code             = "m",
  group_tag            = NULL,
  snv_position         = NULL,
  ref_base             = NULL,
  alt_base             = NULL,
  max_reads            = 200L,
  per_group_downsample = FALSE
)
```

### Step 2 — Replace section 6 (global downsample) with conditional logic

Current code does an unconditional global `sample(n_reads, max_reads)`. Replace with:

```r
# --- 6. Downsample if needed ---
n_reads <- nrow(reads)
keep_local <- seq_len(n_reads)

active_grouping <- per_group_downsample &&
                   !is.null(group_tag) && "group" %in% names(reads)

if (active_grouping) {
  groups <- unique(reads$group[!is.na(reads$group)])
  keep_local <- unlist(lapply(groups, function(g) {
    idx <- which(reads$group == g)
    if (length(idx) > max_reads) idx <- sort(sample(idx, max_reads))
    idx
  }), use.names = FALSE)
  keep_local <- sort(keep_local)
  reads <- reads[keep_local, , drop = FALSE]
  rownames(reads) <- NULL
} else if (n_reads > max_reads) {
  keep_local <- sort(sample(n_reads, max_reads))
  reads <- reads[keep_local, , drop = FALSE]
  rownames(reads) <- NULL
}

keep <- bam_indices[keep_local]
```

`active_grouping` is TRUE only when `per_group_downsample = TRUE` AND grouping is active (either `group_tag`-based or SNV-based, both of which set `group_tag` to non-NULL and add a `group` column to `reads`). When `per_group_downsample = FALSE` (default), the code falls through to the existing global cap — byte-for-byte identical to current behaviour.

---

## Edge cases

| Scenario | Behaviour |
|---|---|
| `per_group_downsample = FALSE` (default) | Global cap always applies; no change from current behaviour. |
| `per_group_downsample = TRUE`, ungrouped data | `active_grouping` is FALSE; falls back to global cap. |
| `per_group_downsample = TRUE`, Group A=300, Group B=50, `max_reads=200` | Sample 200 from A, keep all 50 from B → 250 total. |
| Group smaller than cap | All reads from that group kept. |
| SNV path with `per_group_downsample = TRUE` | Falls into per-group path correctly. |

---

## Testing approach

1. Default (`per_group_downsample = FALSE`) produces exactly `max_reads` reads globally (unchanged behaviour).
2. `per_group_downsample = TRUE` with two unequal groups: each group independently capped at `max_reads`.
3. `per_group_downsample = TRUE` with ungrouped data: falls back to global cap (no group column).
4. `result$sites` contains only `read_name` values present in `result$reads`.
