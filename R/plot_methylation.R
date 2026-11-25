# Internal helpers for interrupting the smooth line at consensus deletion regions

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
# deletion. A read covers the interval when its aligned span overlaps it
# (reads$start/end); carriers always count. Reads that do not reach the
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

# One identity string per smooth line, built from the columns that tell lines
# apart (NA-safe: interaction() would turn every row with an NA group into a
# single NA level, merging e.g. two ungrouped samples' CI ribbons into one
# polygon). Used as the `group` of both the line and its CI ribbon.
.line_id = function(df, cols) {
  if (length(cols) == 0L || nrow(df) == 0L) return(rep("all", nrow(df)))
  parts = lapply(cols, function(col) {
    v = as.character(df[[col]])
    ifelse(is.na(v), "<NA>", v)
  })
  do.call(paste, c(parts, sep = "\r"))
}

# Common ggplot2 layers for the smoothed modification probability panel.
# Returns a list of scales/coords/theme components shared by every branch.
# `y_label` differs between call modes: continuous mode plots a mean probability,
# binary mode plots a fraction of methylated calls (see .smooth_y_label()).
# Fitted values are clamped to [0, 1] in .loess_predict(), so the y limits
# never drop part of a line.
.smooth_panel_base = function(region_start, region_end,
                               y_label = "Mean modification\nprobability") {
  list(
    ggplot2::scale_y_continuous(limits = c(0, 1), name = y_label),
    ggplot2::scale_x_continuous(labels = scales::comma_format()),
    ggplot2::coord_cartesian(xlim = c(region_start, region_end)),
    theme_ggmethylation()
  )
}

# y-axis label for the smooth panel, and for the delta panel below it. In binary
# mode the panel aggregates 0/1 calls, so the quantity is a fraction of
# methylated calls rather than a mean probability.
.smooth_y_label = function(binary) {
  if (binary) {
    "Fraction\nmethylated"
  } else {
    "Mean modification\nprobability"
  }
}

# The delta panel is the shortest in the stack, so this label has to stay short
# enough to fit its height - a longer one overflows vertically and collides with
# the smooth panel's y-axis title. The "(group2 - group1)" line it used to carry
# is now redundant: the area is filled with the colour of whichever group is
# higher, so the sign is readable off the plot itself.
.delta_y_label = function(binary) {
  if (binary) {
    "\u0394 fraction\nmethylated"
  } else {
    "\u0394 methylation"
  }
}

# Stop early on a `sort_by` column that doesn't exist: order() on a NULL column
# returns integer(0) and would silently drop every read.
.check_sort_by = function(sort_by, reads) {
  missing_cols = setdiff(sort_by, names(reads))
  if (length(missing_cols) > 0L) {
    stop("`sort_by` column(s) not found in `data$reads`: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
}

# Centred text placeholder used for empty panels.
.message_plot = function(label, size) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = label, size = size) +
    ggplot2::theme_void()
}

# Site frame used for every aggregate (smooth panel, delta panel, per-read
# means). In binary mode these aggregate 0/1 calls instead of raw
# probabilities, so they report a fraction of methylated calls. The raw
# `data$sites` is still what build_read_panel() colours, so it can keep the
# ambiguous category.
.sites_for_aggregation = function(sites, call_threshold, call_ambiguous) {
  if (is.null(call_threshold)) sites else .binarize_sites(sites, call_threshold, call_ambiguous)
}

# Sort reads and pack them into lanes (per group when grouped, with a blank
# separator lane between groups). Adds `mean_mod_prob` and `lane` to
# `data$reads`; returns list(data, separator_lanes).
.prepare_reads = function(data, sites_agg, sort_by = NULL) {
  reads = data$reads
  reads$mean_mod_prob = .read_mean_mod_prob(sites_agg, reads$read_name)
  reads$mean_mod_prob[is.na(reads$mean_mod_prob)] = 0

  if (is.null(sort_by)) {
    sort_by = if (is.null(data$group_tag)) "start" else c("start", "group", "mean_mod_prob")
  }
  .check_sort_by(sort_by, reads)
  reads = reads[do.call(order, lapply(sort_by, function(col) reads[[col]])), , drop = FALSE]
  rownames(reads) = NULL

  separator_lanes = numeric(0)
  if (is.null(data$group_tag)) {
    reads$lane = pack_reads(reads, clip_side = reads$clip_side)
  } else {
    reads$lane = integer(nrow(reads))
    lane_offset = 0L
    for (grp in .ordered_plot_groups(reads$group)) {
      idx = which(.match_plot_group(reads$group, grp))
      reads$lane[idx] = pack_reads(reads[idx, ], clip_side = reads$clip_side[idx]) + lane_offset
      lane_offset = max(reads$lane[idx]) + 2L
      separator_lanes = c(separator_lanes, lane_offset - 1L)
    }
    # No separator after the last group
    separator_lanes = utils::head(separator_lanes, -1L)
  }

  data$reads = reads
  list(data = data, separator_lanes = separator_lanes)
}

# Smooth panel for any combination of identity variables. `colour` and
# `linetype` name columns of `smoothed` (or NULL); with no `colour` a single
# line is drawn in `fixed_colour`. `labels` is passed to labs() (e.g.
# list(colour = "Group")); `colour_values` feeds a manual colour (and matching
# ribbon fill) scale. `colour_guide = FALSE` hides the colour legend (used when
# colour encodes the group, whose legend the read panel already shows).
#
# Every line gets its own `group` (see .line_id()), shared with its CI ribbon,
# so overlapping bands never collapse into one self-crossing polygon. The
# ribbon is layer 1 (drawn behind) and the lines layer 2; without a ribbon the
# lines are layer 1. NA rows are kept: they are the breaks that
# .insert_deletion_breaks() adds, and geom_line()/geom_ribbon() split there.
.build_smooth_panel = function(smoothed, region_start, region_end, y_label,
                               show_ci, colour = NULL, linetype = NULL,
                               labels = list(), colour_values = NULL,
                               linetype_values = NULL, linetype_name = NULL,
                               fixed_colour = .PROB_GRADIENT$high,
                               colour_guide = TRUE) {
  col_sym = .data_col
  smoothed$.line = .line_id(smoothed, c(colour, linetype))
  aes_args = list(x = quote(.data$position), y = quote(.data$mean_prob),
                  group = quote(.data$.line))
  if (!is.null(colour))   aes_args$colour   = col_sym(colour)
  if (!is.null(linetype)) aes_args$linetype = col_sym(linetype)

  p = ggplot2::ggplot(smoothed, do.call(ggplot2::aes, aes_args))

  has_ci = isTRUE(show_ci) && all(c("lower", "upper") %in% names(smoothed)) &&
    any(!is.na(smoothed$lower) & !is.na(smoothed$upper))
  if (has_ci) {
    rib_aes = list(x = quote(.data$position), ymin = quote(.data$lower),
                   ymax = quote(.data$upper), group = quote(.data$.line))
    if (!is.null(colour)) rib_aes$fill = col_sym(colour)
    rib_args = list(mapping = do.call(ggplot2::aes, rib_aes), alpha = 0.2,
                    colour = NA, inherit.aes = FALSE, show.legend = FALSE,
                    na.rm = TRUE)
    if (is.null(colour)) rib_args$fill = fixed_colour
    p = p + do.call(ggplot2::geom_ribbon, rib_args)
  }

  p = p + if (is.null(colour)) {
    ggplot2::geom_line(linewidth = 1, colour = fixed_colour, na.rm = TRUE)
  } else {
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE)
  }
  p = p +
    .smooth_panel_base(region_start, region_end, y_label) +
    do.call(ggplot2::labs, c(list(x = "Genomic position (bp)"), labels))
  if (!is.null(linetype_values)) {
    # Keys in the given order, and wide enough to show the dash pattern.
    p = p +
      ggplot2::scale_linetype_manual(values = linetype_values, name = linetype_name,
                                     limits = names(linetype_values)) +
      ggplot2::theme(legend.key.width = ggplot2::unit(1.2, "cm"))
  }
  guide = if (isTRUE(colour_guide)) "legend" else "none"
  if (!is.null(colour_values)) {
    p = p +
      ggplot2::scale_colour_manual(values = colour_values, na.value = .GROUP_NA_COLOUR,
                                   guide = guide) +
      ggplot2::scale_fill_manual(values = colour_values, na.value = .GROUP_NA_COLOUR,
                                 guide = "none")
  } else if (!is.null(colour) && !isTRUE(colour_guide)) {
    p = p + ggplot2::guides(colour = "none")
  }
  p
}

# Modification codes in display order: the object's `mod_code` order,
# restricted to codes that actually occur in `sites` (plus any stray codes).
.present_codes = function(mod_code, sites) {
  seen = unique(as.character(sites$mod_code))
  c(intersect(mod_code, seen), setdiff(seen, mod_code))
}

# Smooth panel for a single sample: one line per group (or one overall line),
# with a linetype per modification code when several codes are present.
.single_smooth_panel = function(data, sites_agg, region_start, region_end,
                                opts) {
  codes = .present_codes(data$mod_code, data$sites)
  multi_code = length(codes) > 1L
  grouped = !is.null(data$group_tag)

  reads = data$reads
  if (!grouped) {
    # `sites_agg` can be empty even though reads are present: in binary mode
    # every site may fall inside the ambiguous band and be dropped. Recycling
    # a length-1 value into a 0-row data.frame is an error, hence the guard.
    sites_agg$group = if (nrow(sites_agg) > 0L) "all" else character(0L)
    reads$group = "all"
  }
  smoothed = smooth_methylation(sites_agg, group_col = "group",
                                mod_code_col = if (multi_code) "mod_code",
                                span = opts$smooth_span)
  smoothed = .apply_deletion_breaks(smoothed, data$cigar_features, reads, "group",
                                    opts$min_indel_size, opts$show_cigar)

  y_label = .smooth_y_label(opts$binary)
  if (grouped) {
    .build_smooth_panel(
      smoothed, region_start, region_end, y_label, opts$show_ci,
      colour = "group", linetype = if (multi_code) "mod_code",
      labels = list(colour = "Group"), colour_values = opts$group_colours,
      linetype_values = if (multi_code) .linetype_values(codes, "modification codes"),
      linetype_name = "Modification",
      colour_guide = FALSE
    )
  } else if (multi_code) {
    .build_smooth_panel(smoothed, region_start, region_end, y_label, opts$show_ci,
                        colour = "mod_code", labels = list(colour = "Modification"))
  } else {
    .build_smooth_panel(smoothed, region_start, region_end, y_label, opts$show_ci)
  }
}

# Shared smooth panel for multi-sample data: one line per sample (and per
# modification code when several are present), or, when any sample is
# grouped, colour by group and linetype by sample (or sample x code).
.multi_smooth_panel = function(samples, prepared, region_start, region_end,
                               opts) {
  codes = unique(unlist(lapply(samples, function(s) .present_codes(s$mod_code, s$sites))))
  multi_code = length(codes) > 1L

  smoothed = do.call(rbind, Map(function(s, prep, nm) {
    if (is.null(prep) || nrow(prep$sites_agg) == 0L) return(NULL)
    # One line per group of a grouped sample (unphased reads get no line, as
    # in single-sample plots), or a single line for an ungrouped one.
    line_key = function(df) {
      if (is.null(s$group_tag)) rep(nm, nrow(df)) else as.character(df$group)
    }
    sites = prep$sites_agg
    reads = s$reads
    sites$smooth_key = line_key(sites)
    reads$smooth_key = line_key(reads)

    sm = smooth_methylation(sites, group_col = "smooth_key",
                            mod_code_col = if (multi_code) "mod_code",
                            span = opts$smooth_span)
    sm = .apply_deletion_breaks(sm, s$cigar_features, reads, "smooth_key",
                                opts$min_indel_size, opts$show_cigar)
    sm$sample = rep(nm, nrow(sm))
    sm$group = if (is.null(s$group_tag)) rep(NA_character_, nrow(sm)) else sm$smooth_key
    if (!multi_code) sm$mod_code = rep(codes[1L], nrow(sm))
    sm$smooth_key = NULL
    sm
  }, samples, prepared, names(samples)))

  if (is.null(smoothed) || nrow(smoothed) == 0L) {
    return(.message_plot("No methylation sites", size = 4))
  }
  rownames(smoothed) = NULL
  sample_levels = names(samples)
  smoothed$sample = factor(smoothed$sample, levels = sample_levels)

  y_label = .smooth_y_label(opts$binary)
  if (any(!vapply(samples, function(s) is.null(s$group_tag), logical(1L)))) {
    # Colour = group; linetype = sample, or sample x code with several codes.
    if (multi_code) {
      lt_levels = as.vector(outer(codes, sample_levels,
                                  function(cd, sm) paste(sm, cd)))
      smoothed$sample_code = factor(paste(smoothed$sample, smoothed$mod_code),
                                    levels = lt_levels)
      lt_col = "sample_code"
      lt_name = "Sample / modification"
    } else {
      lt_levels = sample_levels
      lt_col = "sample"
      lt_name = "Sample"
    }
    .build_smooth_panel(smoothed, region_start, region_end, y_label, opts$show_ci,
                        colour = "group", linetype = lt_col,
                        labels = list(colour = "Group"),
                        colour_values = opts$group_colours,
                        linetype_values = .linetype_values(lt_levels, "line types (samples x codes)"),
                        linetype_name = lt_name,
                        colour_guide = FALSE)
  } else {
    .build_smooth_panel(smoothed, region_start, region_end, y_label, opts$show_ci,
                        colour = "sample", linetype = if (multi_code) "mod_code",
                        labels = list(colour = "Sample"),
                        linetype_values = if (multi_code) .linetype_values(codes, "modification codes"),
                        linetype_name = "Modification")
  }
}

# Delta panel for exactly two groups and one modification code, or NULL (with
# a message) otherwise.
.delta_panel = function(data, sites_agg, region_start, region_end, opts) {
  codes = .present_codes(data$mod_code, sites_agg)
  if (length(codes) > 1L) {
    message(sprintf(
      "Delta track needs a single modification code but the data has %s; skipping. %s",
      paste(codes, collapse = ", "),
      "Read the data with one `mod_code` to show it."
    ))
    return(NULL)
  }
  delta_df = .compute_group_delta(sites_agg, "group", span = opts$smooth_span)
  if (is.null(delta_df)) return(NULL)
  # Read the group names off the attribute *now*: breaking the line on
  # deletions rebuilds delta_df, which drops attributes.
  delta_groups = attr(delta_df, "groups")
  if (isTRUE(opts$show_cigar)) {
    delta_df = .break_delta_on_deletions(delta_df, data$cigar_features, data$reads,
                                         opts$min_indel_size)
  }
  # Colour the delta area by whichever group is higher, using the same
  # resolved group colours as the read and smooth panels.
  group_colours = .resolve_group_colours(opts$group_colours, delta_groups)
  .build_delta_panel(delta_df, region_start, region_end,
                     .delta_y_label(opts$binary),
                     fill_pos = unname(group_colours[delta_groups[2L]]),
                     fill_neg = unname(group_colours[delta_groups[1L]]))
}

#' Plot read-level methylation data
#'
#' Creates a ggplot2 visualisation of read-level base modification data,
#' returned as a [patchwork::patchwork] composite: an optional gene track, one
#' read panel per sample (reads as bars with a coloured mark per modification
#' site), a panel with the loess-smoothed mean modification probability (one
#' line per group when the data is grouped, otherwise one overall line), and an
#' optional delta panel (`show_delta`). Setting `call_threshold` switches to
#' binary calls: sites are coloured methylated/unmethylated and the smooth
#' panel plots the fraction of methylated calls instead of the mean probability.
#'
#' Groups are coloured consistently in every panel (read bars, smooth lines,
#' delta area) by one palette resolved from `group_colours`; reads whose group
#' is `NA` (e.g. unphased) are drawn in grey. Legends shared by several panels
#' are collected once at the right of the figure.
#'
#' When the data holds several modification codes (e.g.
#' `mod_code = c("m", "h")`), each read bar is split into horizontal strips,
#' one per code, in the order of `data$mod_code` from top to bottom (the read
#' panel's y-axis title lists the order), and the smooth panel draws one line
#' per code (by line type).
#'
#' @param data A `methylation_data` object returned by [read_methylation()],
#'   or a `multi_methylation_data` object from [merge_methylation()].
#' @param sort_by Character vector of column names from `data$reads` used
#'   to sort reads before packing into lanes (reads are always packed per
#'   group, so group order is fixed). Default NULL uses `"start"` when
#'   ungrouped or `c("start", "group", "mean_mod_prob")` when grouped.
#'   `mean_mod_prob` is the per-read mean modification probability, or the
#'   per-read fraction of methylated calls when `call_threshold` is set.
#' @param colour_low Colour for low modification probability (default
#'   `"#BDBDBD"`).
#' @param colour_high Colour for high modification probability (default
#'   `"#C62828"`).
#' @param colour_ambiguous Colour for ambiguous calls (default `"#78909C"`, a
#'   slate blue-grey deliberately off the `colour_low`/`colour_high` ramp so it
#'   reads as "no confident call" rather than "intermediate methylation"). Only
#'   used when `call_threshold` and `call_ambiguous` are set.
#' @param line_width Positive number: linewidth of modification site lines
#'   (default 0.2).
#' @param colour_strand Logical. When `TRUE`, read bars are coloured by strand
#'   (`"+"` = forward, `"-"` = reverse). Ignored when data is grouped; group
#'   colour takes precedence. Default `FALSE`.
#' @param strand_colours Named character vector with `"+"` and `"-"` entries
#'   giving bar colours for each strand. Only used when `colour_strand = TRUE`
#'   and data is ungrouped. Default `c("+" = "#4393C3", "-" = "#D6604D")`.
#' @param group_colours Character vector of colours per group, used for the
#'   read bars, the smooth lines and the delta area alike. Defaults to
#'   `c("1" = "#0072B2", "2" = "#E69F00")` (the Okabe-Ito colorblind-safe
#'   blue/orange pair), matching HP haplotype tags. Groups the vector does not
#'   name (e.g. `"REF"`/`"ALT"` from SNV grouping, or a third haplotype) get
#'   further Okabe-Ito colours (with a message when you supplied the vector).
#'   An unnamed vector is matched to the sorted group names by position and
#'   needs at least one colour per group. Pass `NULL` to use ggplot2's default
#'   hue palette.
#' @param smooth_span `NULL` or a positive number: the loess smoothing span for
#'   the smooth panel. When `NULL` (default), an adaptive span is computed per
#'   group as `max(0.15, min(0.75, 15 / n_unique_sites))`, targeting ~15 data
#'   points per local fit regardless of region size or CpG density. Pass an
#'   explicit value (e.g. `0.3`) to use a fixed span.
#' @param panel_heights Numeric vector giving relative heights of the panels,
#'   top to bottom: the gene track (when `annotations` is given), one read
#'   panel per sample, the smooth panel, and the delta panel (when shown).
#'   Its length must equal the number of panels. Pass `NULL` (default) to use
#'   the built-in defaults.
#' @param annotations A `gene_annotations` object returned by
#'   [read_annotations()], or `NULL` (default). When provided, a gene
#'   annotation track is placed at the top of the composite figure (above the
#'   read panel and the smooth panel).
#' @param variants A `variant_data` object returned by [read_variants()], or
#'   `NULL` (default). When provided, VCF records are overlaid on the read
#'   panel(s):
#'   * **SNV**: a small black asterisk at the variant position on each read
#'     whose base there equals the ALT allele. Reads with the REF or another
#'     base are not marked; no base letters are drawn.
#'   * **DEL, DUP, INV** (symbolic SVs with an `END`): a styled segment on
#'     every read overlapping the SV span, clipped to the read and the region:
#'     DEL a dark red line with end ticks, DUP two teal parallel lines, INV a
#'     purple dashed line. Chevrons (« ») mark a span cut at the region
#'     edge.
#'   * **BND**: one text label `"⋈ mate_chrom:mate_pos"` at the bottom of
#'     the read panel, centred on the breakend position. When
#'     `show_supplementary = TRUE`, the SA marker of each read whose
#'     supplementary alignment lies on the BND mate chromosome within the
#'     object's `bnd_match_tol` of the mate position gets a dark magenta
#'     border.
#'   * Short insertions, deletions and MNVs, and symbolic types other than
#'     DEL/DUP/INV (e.g. `<CNV>`), are not drawn.
#'
#'   No vertical lines are drawn for variants. Requires that `data` was
#'   produced by a current version of [read_methylation()], which stores
#'   sequences and CIGARs.
#' @param show_cigar Logical. When `TRUE` (default), structural variants from
#'   CIGAR strings (insertions, deletions) are displayed on reads. Insertions
#'   appear as purple I-beam markers spanning the full height of the read bar;
#'   deletions open a gap in the read bar and draw a thin black line. Smooth
#'   lines (and the delta area) are interrupted where at least 75% of a
#'   group's reads covering a deletion carry it. Set to `FALSE` to hide these
#'   features.
#' @param min_indel_size Non-negative number. Minimum size (in bp) for
#'   insertions and deletions to be displayed when `show_cigar = TRUE`. Indels
#'   smaller than this threshold are hidden to reduce visual clutter from
#'   common small indels. Default `50`.
#' @param show_supplementary Logical. When `TRUE` (default), reads with a
#'   supplementary alignment (SA BAM tag) get a small marker at each clipped
#'   end, filled by the chromosome of the SA partner. The marker is drawn
#'   outside the read bar, so no modification site is hidden.
#' @param show_ci Logical. When `TRUE` (default), a shaded ribbon showing the
#'   loess confidence interval (`lower`/`upper` from [smooth_methylation()])
#'   is drawn behind each smooth line in the smooth panel. Has no effect when
#'   fewer than 4 unique positions are available for a group/code (no CI is
#'   computed in that case). Set to `FALSE` to hide the ribbon.
#' @param call_threshold `NULL` (default) or a number in `[0, 1]`. With
#'   `NULL`, modification sites are coloured by a continuous probability
#'   gradient (`colour_low` to `colour_high`) and the smooth panel plots the
#'   mean modification probability per position. A number switches to
#'   **binary calls**: sites with probability at or above the threshold are
#'   methylated, the rest unmethylated (or ambiguous, see `call_ambiguous`),
#'   coloured with a discrete scale. The same split drives the read panel and
#'   every aggregate, so the two always agree.
#'
#'   Binary calls change what every aggregate reports: the smooth panel, the
#'   delta panel (`show_delta`) and the per-read `mean_mod_prob` used for
#'   sorting all aggregate 0/1 calls rather than raw probabilities, so the
#'   smooth panel plots the **fraction of calls that are methylated** and its
#'   y-axis is relabelled accordingly. The two quantities genuinely differ: a
#'   position where every read reports probability 0.8 has a mean probability of
#'   0.8 but a methylated fraction of 1.0. The fraction is the conventional
#'   percent-methylation / beta value; the continuous mean is systematically
#'   pulled toward 0.5 by basecaller uncertainty.
#'
#'   Because binarised values are 0 or 1, positions with thin coverage
#'   contribute coarser values than with continuous probabilities, so the
#'   binary curve is noisier where few reads overlap. The loess span averages
#'   over neighbouring positions. Note that the loess fit weights every
#'   position equally regardless of how many reads cover it, so the confidence
#'   ribbon (`show_ci`) reflects the scatter and density of positions, not read
#'   coverage: it does not widen where coverage is thin. Fitted values are
#'   clamped to `[0, 1]`.
#' @param call_ambiguous `NULL` (default) for a hard threshold, or a single
#'   non-negative number `w`: sites whose probability lies in the open band
#'   `(call_threshold - w, call_threshold + w)` are labelled "ambiguous"
#'   rather than methylated/unmethylated, except a probability exactly equal
#'   to `call_threshold`, which is always methylated. Ambiguous sites are
#'   excluded from both the numerator and the denominator of the methylated
#'   fraction, so the smooth panel reports the fraction methylated *among
#'   confident calls*. The default `NULL` (or `0`) discards nothing — every
#'   site is called at `call_threshold`. Requires `call_threshold`.
#' @param show_delta Logical. When `TRUE`, and grouping yields exactly two
#'   groups, an additional panel is appended below the smooth panel showing
#'   the signed difference in loess-smoothed modification probability between
#'   the two groups (group2 - group1), filled with the colour of whichever
#'   group is higher. Requires `data` to be grouped (`group_tag` or `snv` in
#'   [read_methylation()]) with exactly two distinct group values and a single
#'   modification code; otherwise the panel is skipped with a message. Not
#'   supported for multi-sample data (`multi_methylation_data`); the argument
#'   is accepted there but ignored (with a message). Default `FALSE`.
#'
#' @return A [patchwork::patchwork] composite (also for ungrouped data), or,
#'   when a single-sample `data` has no reads, a placeholder [ggplot2::ggplot]
#'   saying so.
#'
#' @examples
#' \dontrun{
#' md = read_methylation("sample.bam", "chr1:1000-2000")
#' plot_methylation(md)
#'
#' md = read_methylation("sample.bam", "chr1:1000-2000", group_tag = "HP")
#' plot_methylation(md, group_colours = c("1" = "steelblue", "2" = "coral"))
#' }
#'
#' @export
plot_methylation = function(data, sort_by = NULL,
                            colour_low = .PROB_GRADIENT$low,
                            colour_high = .PROB_GRADIENT$high,
                            colour_ambiguous = .CALL_AMBIGUOUS_DEFAULT,
                            line_width = 0.2,
                            colour_strand = FALSE,
                            strand_colours = c("+" = "#4393C3", "-" = "#D6604D"),
                            group_colours = .GROUP_PALETTE_DEFAULT,
                            smooth_span = NULL,
                            panel_heights = NULL,
                            annotations = NULL,
                            variants = NULL,
                            show_cigar = TRUE,
                            min_indel_size = 50L,
                            show_supplementary = TRUE,
                            show_ci = TRUE,
                            call_threshold = NULL,
                            call_ambiguous = NULL,
                            show_delta = FALSE) {
  # --- 1. Validate input ---
  multi = inherits(data, "multi_methylation_data")
  if (!multi && !inherits(data, "methylation_data")) {
    stop(
      "`data` must be a `methylation_data` or `multi_methylation_data` object.",
      call. = FALSE
    )
  }
  .validate_plot_args(
    colour_strand = colour_strand, strand_colours = strand_colours,
    call_threshold = call_threshold, call_ambiguous = call_ambiguous,
    smooth_span = smooth_span, group_colours = group_colours,
    line_width = line_width, min_indel_size = min_indel_size,
    show_cigar = show_cigar, show_supplementary = show_supplementary,
    show_ci = show_ci, show_delta = show_delta
  )
  if (!is.null(variants) && !inherits(variants, "variant_data")) {
    stop("`variants` must be a `variant_data` object returned by read_variants().",
         call. = FALSE)
  }
  if (!is.null(annotations) && !inherits(annotations, "gene_annotations")) {
    stop("`annotations` must be a `gene_annotations` object from read_annotations().",
         call. = FALSE)
  }

  samples = if (multi) data$samples else list(data)
  grouped_samples = !vapply(samples, function(s) is.null(s$group_tag), logical(1L))
  if (isTRUE(colour_strand) && any(grouped_samples)) {
    message("'colour_strand = TRUE' is ignored when reads are grouped; group colour takes precedence.")
  }
  if (multi && isTRUE(show_delta)) {
    message("`show_delta` is not supported for multi-sample data; ignoring.")
  }
  if (!multi && isTRUE(show_delta) && is.null(data$group_tag)) {
    message("`show_delta` needs grouped data (`group_tag` or `snv` in read_methylation()); skipping the delta panel.")
  }
  if (!multi && nrow(data$reads) == 0L) {
    return(.message_plot("No reads in region", size = 5))
  }

  region_start = GenomicRanges::start(data$region)
  region_end   = GenomicRanges::end(data$region)

  # The gene panel only sees coordinates, so the chromosome is checked here.
  ann_region = annotations$region
  if (inherits(ann_region, "GRanges") && length(ann_region) == 1L) {
    ann_chrom  = as.character(GenomicRanges::seqnames(ann_region))
    data_chrom = as.character(GenomicRanges::seqnames(data$region))
    if (!identical(ann_chrom, data_chrom)) {
      warning(sprintf(
        "`annotations` were read for chromosome '%s' but the data are on '%s'.",
        ann_chrom, data_chrom
      ), call. = FALSE)
    }
  }

  # One group palette for every panel and sample, resolved over the union of
  # the groups present in any sample's reads.
  all_groups = unlist(lapply(samples[grouped_samples], function(s) s$reads$group))
  resolved_colours = .resolve_group_colours(group_colours, all_groups)
  # Legend keys shared by every read panel, so patchwork can merge them.
  legend_levels = list(
    group = names(resolved_colours),
    group_na = any(is.na(all_groups)),
    sa = .group_levels(unlist(lapply(samples, function(s) s$reads$sa_chrom)))
  )
  opts = list(smooth_span = smooth_span, min_indel_size = min_indel_size,
              show_cigar = show_cigar, show_ci = show_ci, binary = !is.null(call_threshold),
              group_colours = resolved_colours)

  # --- 2. Sort and pack each sample's reads ---
  prepared = lapply(samples, function(s) {
    if (nrow(s$reads) == 0L) return(NULL)
    sites_agg = .sites_for_aggregation(s$sites, call_threshold, call_ambiguous)
    c(.prepare_reads(s, sites_agg, sort_by), list(sites_agg = sites_agg))
  })

  # --- 3. Read panel(s) ---
  sample_names = if (multi) names(samples) else ""
  read_panels = Map(function(prep, nm) {
    if (is.null(prep)) {
      return(.message_plot(sprintf("No reads (%s)", nm), size = 4) + ggplot2::labs(title = nm))
    }
    s = prep$data
    p = build_read_panel(
      data               = s,
      separator_lanes    = prep$separator_lanes,
      region_start       = region_start,
      region_end         = region_end,
      colour_low         = colour_low,
      colour_high        = colour_high,
      colour_ambiguous   = colour_ambiguous,
      line_width         = line_width,
      colour_strand      = colour_strand,
      strand_colours     = strand_colours,
      group_colours      = resolved_colours,
      show_x_axis        = FALSE,
      variant_overlay    = build_variant_overlay(s, variants),
      show_cigar         = show_cigar,
      min_indel_size     = min_indel_size,
      show_supplementary = show_supplementary,
      call_threshold     = call_threshold,
      call_ambiguous     = call_ambiguous,
      legend_levels      = legend_levels
    )
    if (multi) p + ggplot2::labs(title = nm) else p
  }, prepared, sample_names)
  names(read_panels) = NULL

  # --- 4. Smooth panel, plus the delta panel for two groups ---
  p_delta = NULL
  if (multi) {
    p_smooth = .multi_smooth_panel(samples, prepared, region_start, region_end, opts)
  } else {
    prep = prepared[[1L]]
    p_smooth = .single_smooth_panel(prep$data, prep$sites_agg, region_start, region_end, opts)
    if (!is.null(data$snv_position)) {
      p_smooth = p_smooth +
        ggplot2::geom_vline(
          xintercept = data$snv_position,
          linetype = "dashed", colour = "black", linewidth = 0.5
        )
    }
    if (isTRUE(show_delta) && !is.null(data$group_tag)) {
      p_delta = .delta_panel(prep$data, prep$sites_agg, region_start, region_end, opts)
      if (!is.null(p_delta)) {
        # The smooth panel is no longer the bottom-most; hide its x-axis title
        p_smooth = p_smooth + ggplot2::theme(axis.title.x = ggplot2::element_blank())
      }
    }
  }

  # --- 5. Gene track on top, then combine ---
  p_gene = if (!is.null(annotations)) {
    build_gene_panel(annotations, region_start, region_end, bottom_label = FALSE)
  }
  panels = c(if (!is.null(p_gene)) list(p_gene), read_panels, list(p_smooth),
             if (!is.null(p_delta)) list(p_delta))

  heights = panel_heights
  if (is.null(heights)) {
    heights = if (multi) {
      c(if (!is.null(p_gene)) 0.6, rep(3, length(samples)), 1)
    } else {
      c(if (!is.null(p_gene)) 0.08, 1, 0.25, if (!is.null(p_delta)) 0.2)
    }
  } else if (length(heights) != length(panels)) {
    stop(sprintf(
      "`panel_heights` has length %d but there are %d panels.",
      length(heights), length(panels)
    ), call. = FALSE)
  }

  # Identical legends (modification scale, group, SA partner) from several
  # panels are merged into one.
  patchwork::wrap_plots(panels, ncol = 1, heights = heights, guides = "collect")
}

# Argument checks for plot_methylation(), each with a message naming the
# argument (instead of an R-level "missing value where TRUE/FALSE needed" or
# a silent fallback further down).
.validate_plot_args = function(colour_strand, strand_colours, call_threshold,
                               call_ambiguous, smooth_span, group_colours,
                               line_width, min_indel_size, show_cigar,
                               show_supplementary, show_ci, show_delta) {
  is_number = function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  flags = list(colour_strand = colour_strand, show_cigar = show_cigar,
               show_supplementary = show_supplementary, show_ci = show_ci,
               show_delta = show_delta)
  for (nm in names(flags)) {
    x = flags[[nm]]
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop(sprintf("`%s` must be TRUE or FALSE.", nm), call. = FALSE)
    }
  }
  if (isTRUE(colour_strand) && !all(c("+", "-") %in% names(strand_colours))) {
    stop("'strand_colours' must be a named vector with '+' and '-' entries.", call. = FALSE)
  }
  if (!is.null(call_threshold) &&
      !(is_number(call_threshold) && call_threshold >= 0 && call_threshold <= 1)) {
    stop("`call_threshold` must be NULL or a single number in [0, 1].", call. = FALSE)
  }
  if (!is.null(call_ambiguous)) {
    if (is.null(call_threshold)) {
      stop("`call_ambiguous` only applies to binary calls; set `call_threshold` too.",
           call. = FALSE)
    }
    if (!(is_number(call_ambiguous) && call_ambiguous >= 0)) {
      stop("`call_ambiguous` must be NULL or a single non-negative number.", call. = FALSE)
    }
  }
  if (!is.null(smooth_span) && !(is_number(smooth_span) && smooth_span > 0)) {
    stop("`smooth_span` must be NULL or a single positive number (e.g. 0.3).",
         call. = FALSE)
  }
  if (!(is_number(line_width) && line_width > 0)) {
    stop("`line_width` must be a single positive number.", call. = FALSE)
  }
  if (!(is_number(min_indel_size) && min_indel_size >= 0)) {
    stop("`min_indel_size` must be a single non-negative number.", call. = FALSE)
  }
  if (!is.null(group_colours)) {
    if (!is.character(group_colours) || length(group_colours) == 0L ||
        anyNA(group_colours)) {
      stop("`group_colours` must be NULL or a character vector of colours without NAs.",
           call. = FALSE)
    }
    bad = group_colours[!vapply(group_colours, function(col) {
      !inherits(tryCatch(grDevices::col2rgb(col), error = function(e) e), "error")
    }, logical(1L))]
    if (length(bad) > 0L) {
      stop(sprintf("`group_colours` contains invalid colour(s): %s.",
                   paste(bad, collapse = ", ")), call. = FALSE)
    }
    nms = names(group_colours)
    if (!is.null(nms) && (any(is.na(nms) | !nzchar(nms)) || anyDuplicated(nms))) {
      stop("`group_colours` must be either unnamed or named with unique, non-empty group names.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}
