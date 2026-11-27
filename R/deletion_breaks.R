# Smooth-line breaks at consensus deletions: deletion ranges per group and NA-break insertion into smoothed and delta data.

.ordered_plot_groups = function(groups) {
  ordered = sort(unique(groups[!is.na(groups)]))
  if (any(is.na(groups))) {
    c(ordered, NA_character_)
  } else {
    ordered
  }
}

.match_plot_group = function(values, group) {
  if (is.na(group)) {
    is.na(values)
  } else {
    !is.na(values) & values == group
  }
}

# For each group, returns genomic intervals [del_start, del_end] where at least
# `threshold` fraction of the group's reads *covering* the interval carry a
# deletion. A read covers the interval when its span (reads$start/end) spans
# the whole interval; carriers always count. Reads that do not reach the
# deletion say nothing about it, so in a wide window with short reads they
# must not dilute the fraction. No further minimum is applied: a lone carrier
# that is the only read there has no sites inside the deletion, so the smooth
# line has no data to draw across it either.
# cigar_features: pre-filtered to type == "D" and min_indel_size by caller
# reads: data.frame with read_name, start, end and the group column (without
#   start/end every read of the group counts as covering)
.consensus_deletion_ranges = function(cigar_features, reads, group_col,
                                       threshold = 0.75) {
  empty = data.frame(
    del_start = integer(0L),
    del_end   = integer(0L),
    stringsAsFactors = FALSE
  )
  empty[[group_col]] = character(0L)

  dels = cigar_features[cigar_features$type == "D", , drop = FALSE]
  if (nrow(dels) == 0L) return(empty)

  # Attach group membership to each deletion row
  dels = merge(
    dels,
    reads[, c("read_name", group_col), drop = FALSE],
    by = "read_name", all.x = FALSE
  )
  if (nrow(dels) == 0L) return(empty)

  groups = unique(dels[[group_col]])
  groups = groups[!is.na(groups)]
  result_list = vector("list", length(groups))

  for (k in seq_along(groups)) {
    grp  = groups[k]
    grp_mask = .match_plot_group(dels[[group_col]], grp)
    sub  = dels[grp_mask, , drop = FALSE]
    grp_reads = reads[.match_plot_group(reads[[group_col]], grp), , drop = FALSE]
    if (nrow(grp_reads) == 0L) next
    has_span = all(c("start", "end") %in% names(grp_reads))

    sub = sub[order(sub$ref_start), , drop = FALSE]

    # Sweep-merge overlapping deletion intervals, accumulating read names
    merged_starts = integer(0L)
    merged_ends   = integer(0L)
    merged_reads  = list()

    cur_start = sub$ref_start[1L]
    cur_end   = sub$ref_end[1L]
    cur_rds   = sub$read_name[1L]

    for (i in seq_len(nrow(sub))[-1L]) {
      if (sub$ref_start[i] <= cur_end + 1L) {
        cur_end = max(cur_end, sub$ref_end[i])
        cur_rds = c(cur_rds, sub$read_name[i])
      } else {
        merged_starts = c(merged_starts, cur_start)
        merged_ends   = c(merged_ends,   cur_end)
        merged_reads  = c(merged_reads,  list(unique(cur_rds)))
        cur_start = sub$ref_start[i]
        cur_end   = sub$ref_end[i]
        cur_rds   = sub$read_name[i]
      }
    }
    merged_starts = c(merged_starts, cur_start)
    merged_ends   = c(merged_ends,   cur_end)
    merged_reads  = c(merged_reads,  list(unique(cur_rds)))

    keep = vapply(seq_along(merged_reads), function(m) {
      carriers = merged_reads[[m]]
      covering = if (has_span) {
        # Reads that span the whole deletion: only they could show sequence
        # across it, so only they count against the carriers.
        grp_reads$read_name[grp_reads$start <= merged_starts[m] &
                              grp_reads$end >= merged_ends[m]]
      } else {
        grp_reads$read_name
      }
      length(carriers) / length(union(carriers, covering)) >= threshold
    }, logical(1L))

    if (!any(keep)) next

    df = data.frame(
      del_start = merged_starts[keep],
      del_end   = merged_ends[keep],
      stringsAsFactors = FALSE
    )
    df[[group_col]] = grp
    result_list[[k]] = df
  }

  out = do.call(rbind, Filter(Negate(is.null), result_list))
  if (is.null(out)) return(empty)
  rownames(out) = NULL
  out
}

# Insert NA breaks into a smoothed data frame at each consensus deletion interval.
# Sentinel NA rows are added at del_start - 0.5 and del_end + 0.5 to guarantee
# a visible gap even when no grid point falls inside the deletion.
.insert_deletion_breaks = function(smoothed, deletion_ranges, group_col) {
  if (nrow(deletion_ranges) == 0L) return(smoothed)

  id_cols = setdiff(names(smoothed), c("position", "mean_prob", "lower", "upper"))
  sentinel_list = vector("list", nrow(deletion_ranges))

  for (i in seq_len(nrow(deletion_ranges))) {
    grp_val   = deletion_ranges[[group_col]][i]
    del_start = deletion_ranges$del_start[i]
    del_end   = deletion_ranges$del_end[i]

    # Mask grid points within the deletion
    in_grp = smoothed[[group_col]] == grp_val
    in_del = smoothed$position >= del_start & smoothed$position <= del_end
    smoothed$mean_prob[in_grp & in_del] = NA_real_
    if ("lower" %in% names(smoothed)) smoothed$lower[in_grp & in_del] = NA_real_
    if ("upper" %in% names(smoothed)) smoothed$upper[in_grp & in_del] = NA_real_

    # Build sentinel rows — one pair per unique line-identity combo in this group
    grp_rows  = smoothed[in_grp, id_cols, drop = FALSE]
    templates = unique(grp_rows)

    if (nrow(templates) == 0L) next

    sentinels = vector("list", nrow(templates) * 2L)
    for (j in seq_len(nrow(templates))) {
      s1 = templates[j, , drop = FALSE]
      s1$position  = del_start - 0.5
      s1$mean_prob = NA_real_
      # NOTE: `templates` is derived from `id_cols`, which now excludes
      # lower/upper, so s1/s2 never carry those columns yet at this point.
      # Check against `smoothed` (the source of truth for which columns
      # exist) rather than `s1`/`s2`, then add the columns as NA so the
      # final rbind()/column-select against names(smoothed) succeeds.
      if ("lower" %in% names(smoothed)) { s1$lower = NA_real_; s1$upper = NA_real_ }
      s2 = templates[j, , drop = FALSE]
      s2$position  = del_end + 0.5
      s2$mean_prob = NA_real_
      if ("lower" %in% names(smoothed)) { s2$lower = NA_real_; s2$upper = NA_real_ }
      sentinels[[2L * j - 1L]] = s1
      sentinels[[2L * j]]      = s2
    }
    sentinel_list[[i]] = do.call(rbind, sentinels)
  }

  all_sentinels = do.call(rbind, Filter(Negate(is.null), sentinel_list))
  if (!is.null(all_sentinels) && nrow(all_sentinels) > 0L) {
    smoothed = rbind(smoothed, all_sentinels[, names(smoothed), drop = FALSE])
  }

  smoothed = smoothed[order(smoothed[[group_col]], smoothed$position), , drop = FALSE]
  rownames(smoothed) = NULL
  smoothed
}

# Deletions long enough to display (and to break smooth lines on).
.large_deletions = function(cigar_features, min_indel_size) {
  if (is.null(cigar_features)) return(data.frame())
  cigar_features[
    cigar_features$type == "D" & cigar_features$length >= min_indel_size,
    , drop = FALSE
  ]
}

# Break smooth lines at consensus deletions when show_cigar is TRUE.
.apply_deletion_breaks = function(smoothed, cigar_features, reads, group_col,
                                  min_indel_size, show_cigar) {
  if (!isTRUE(show_cigar)) return(smoothed)
  ranges = .consensus_deletion_ranges(
    .large_deletions(cigar_features, min_indel_size), reads, group_col
  )
  .insert_deletion_breaks(smoothed, ranges, group_col)
}

# Break the single delta line wherever *either* group has a consensus
# deletion. The ranges are keyed by the real per-read group values ("1"/"2"),
# so they are computed per group first and only then relabelled to the single
# "delta" line - relabelling the delta line up front would match no range.
.break_delta_on_deletions = function(delta_df, cigar_features, reads,
                                     min_indel_size) {
  ranges = .consensus_deletion_ranges(
    .large_deletions(cigar_features, min_indel_size), reads, "group"
  )
  if (nrow(ranges) == 0L) return(delta_df)
  ranges$group = "delta"
  line = data.frame(position = delta_df$position, mean_prob = delta_df$delta,
                    group = "delta", stringsAsFactors = FALSE)
  line = .insert_deletion_breaks(line, ranges, "group")
  data.frame(position = line$position, delta = line$mean_prob,
             sign = .delta_sign(line$mean_prob), stringsAsFactors = FALSE)
}
