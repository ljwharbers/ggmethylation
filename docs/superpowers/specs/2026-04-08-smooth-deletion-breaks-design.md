# Design: Interrupt LOESS smooth line at consensus deletions

**Date:** 2026-04-08  
**Branch:** fix/issues-3-7-gene-display-cigar-variants  
**Status:** Approved

---

## Summary

When `show_cigar = TRUE`, the LOESS-smoothed methylation line in the bottom panel
should visually break (gap) at genomic positions where ≥75% of reads in a group
carry a deletion. This prevents the smooth line from interpolating across regions
that are deleted in the majority of reads, which would be misleading.

---

## Behaviour

- Triggered by: `show_cigar = TRUE` in `plot_methylation()`
- Threshold: hardcoded at 0.75 (75% of reads per group)
- Minimum deletion size: respects existing `min_indel_size` parameter — deletions
  smaller than `min_indel_size` are excluded from consensus computation
- Applies to: all smooth line variants (ungrouped, grouped, multi-sample)
- No new user-facing parameters

---

## Architecture

Two internal helpers added to `plot_methylation.R`:

### `.consensus_deletion_ranges(cigar_features, reads, group_col, threshold = 0.75)`

**Input:**
- `cigar_features` — data frame from `methylation_data$cigar_features`; rows with
  `type == "D"` and columns `read_name`, `ref_start`, `ref_end`, `length`
- `reads` — `data$reads` data frame with `read_name` and optionally `group`
- `group_col` — name of the group column (`"group"` or `"smooth_key"` for multi-sample)
- `threshold` — fraction of reads required (default 0.75)

**Logic:**
1. Join deletions to reads via `read_name` to get group membership per deletion
2. For each group, count total reads (`n_reads`)
3. For each deletion interval `[ref_start, ref_end]`, count how many distinct reads
   in the group carry it (`n_del`)
4. Merge overlapping deletion intervals within a group (sort by `ref_start`, union
   overlapping ranges)
5. For each merged interval, keep only those where `n_del / n_reads >= threshold`
6. Return data frame with columns: `group`, `del_start`, `del_end`

**Edge cases:**
- No deletions in `cigar_features` → return zero-row data frame, no breaks applied
- Group with a single read → threshold still applied (1/1 = 100% ≥ 75%)
- Deletions that span the entire region → entire smooth line becomes NA (line
  disappears for that group; acceptable)

### `.insert_deletion_breaks(smoothed, deletion_ranges, group_col)`

**Input:**
- `smoothed` — output of `smooth_methylation()`, columns: `position`, `mean_prob`,
  and the group column
- `deletion_ranges` — output of `.consensus_deletion_ranges()`, zero or more rows
- `group_col` — name of the group column

**Logic:**
1. For each row in `deletion_ranges`, find rows in `smoothed` for the matching group
   where `position >= del_start & position <= del_end`
2. Set `mean_prob = NA` for those rows
3. Additionally insert two explicit NA sentinel rows at `position = del_start - 0.5`
   and `position = del_end + 0.5` to guarantee a visible break even when no grid
   point falls inside a narrow deletion
4. Re-sort by `(group_col, position)`
5. Return modified `smoothed`

**Why half-integer offsets:** grid positions are integers (genomic coordinates);
half-integer sentinels cannot collide with real positions and will be clipped by
`scale_x_continuous(limits = ...)` if outside the plot region.

---

## Integration in `plot_methylation()`

Applied immediately after each `smooth_methylation()` call, before the ggplot
is constructed. Pseudocode for the single-sample path:

```r
smoothed <- smooth_methylation(sites_smooth, group_col = "group", span = smooth_span)

if (isTRUE(show_cigar) && !is.null(data$cigar_features) &&
    nrow(data$cigar_features) > 0L) {
  dels <- data$cigar_features[
    data$cigar_features$type == "D" &
      data$cigar_features$length >= min_indel_size, , drop = FALSE
  ]
  if (nrow(dels) > 0L) {
    # reads_grouped: data$reads with a "group" column added.
    # Ungrouped path: add group = "all" to mirror sites_smooth.
    # Grouped path: data$reads already has group column.
    reads_grouped <- data$reads
    if (is.null(data$group_tag)) reads_grouped$group <- "all"
    ranges <- .consensus_deletion_ranges(dels, reads_grouped, "group")
    smoothed <- .insert_deletion_breaks(smoothed, ranges, "group")
  }
}
```

The same pattern applies to:
- Multi-code path (group column = `".smooth_group"`): reads_grouped also needs the
  composite key column added, matching the key format used in `sites_smooth`
- Multi-sample path in `.plot_multi_methylation()` (group column = `"smooth_key"`):
  `combined_sites` already carries `smooth_key`; reads from each sample need to be
  tagged with the same key before passing to the helper

---

## Files Changed

| File | Change |
|------|--------|
| `R/plot_methylation.R` | Add two helper functions; call them after each `smooth_methylation()` invocation |

No changes to `smooth_methylation.R`, `build_read_panel.R`, or any data structures.

---

## Testing

- Unit test for `.consensus_deletion_ranges()`: verify threshold filtering, interval
  merging, ungrouped case, zero-deletion case
- Unit test for `.insert_deletion_breaks()`: verify NA rows inserted at correct
  positions, sentinel rows added, sort order preserved
- No integration test changes needed (BAM-dependent tests unaffected)
