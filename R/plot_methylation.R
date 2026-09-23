# Internal helpers for interrupting the smooth line at consensus deletion regions

.ordered_plot_groups <- function(groups) {
  ordered <- sort(unique(groups[!is.na(groups)]))
  if (any(is.na(groups))) {
    c(ordered, NA_character_)
  } else {
    ordered
  }
}

.match_plot_group <- function(values, group) {
  if (is.na(group)) {
    is.na(values)
  } else {
    !is.na(values) & values == group
  }
}

# For each group, returns genomic intervals [del_start, del_end] where at least
# `threshold` fraction of reads in that group carry a deletion.
# cigar_features: pre-filtered to type == "D" and min_indel_size by caller
# reads: data.frame with read_name and the group column
.consensus_deletion_ranges <- function(cigar_features, reads, group_col,
                                       threshold = 0.75) {
  empty <- data.frame(
    del_start = integer(0L),
    del_end   = integer(0L),
    stringsAsFactors = FALSE
  )
  empty[[group_col]] <- character(0L)

  dels <- cigar_features[cigar_features$type == "D", , drop = FALSE]
  if (nrow(dels) == 0L) return(empty)

  # Attach group membership to each deletion row
  dels <- merge(
    dels,
    reads[, c("read_name", group_col), drop = FALSE],
    by = "read_name", all.x = FALSE
  )
  if (nrow(dels) == 0L) return(empty)

  groups <- unique(dels[[group_col]])
  groups <- groups[!is.na(groups)]
  result_list <- vector("list", length(groups))

  for (k in seq_along(groups)) {
    grp  <- groups[k]
    grp_mask <- .match_plot_group(dels[[group_col]], grp)
    sub  <- dels[grp_mask, , drop = FALSE]
    n_reads <- sum(.match_plot_group(reads[[group_col]], grp))
    if (n_reads == 0L) next

    sub <- sub[order(sub$ref_start), , drop = FALSE]

    # Sweep-merge overlapping deletion intervals, accumulating read names
    merged_starts <- integer(0L)
    merged_ends   <- integer(0L)
    merged_reads  <- list()

    cur_start <- sub$ref_start[1L]
    cur_end   <- sub$ref_end[1L]
    cur_rds   <- sub$read_name[1L]

    for (i in seq_len(nrow(sub))[-1L]) {
      if (sub$ref_start[i] <= cur_end + 1L) {
        cur_end <- max(cur_end, sub$ref_end[i])
        cur_rds <- c(cur_rds, sub$read_name[i])
      } else {
        merged_starts <- c(merged_starts, cur_start)
        merged_ends   <- c(merged_ends,   cur_end)
        merged_reads  <- c(merged_reads,  list(unique(cur_rds)))
        cur_start <- sub$ref_start[i]
        cur_end   <- sub$ref_end[i]
        cur_rds   <- sub$read_name[i]
      }
    }
    merged_starts <- c(merged_starts, cur_start)
    merged_ends   <- c(merged_ends,   cur_end)
    merged_reads  <- c(merged_reads,  list(unique(cur_rds)))

    keep <- vapply(
      merged_reads,
      function(rds) length(rds) / n_reads >= threshold,
      logical(1L)
    )

    if (!any(keep)) next

    df <- data.frame(
      del_start = merged_starts[keep],
      del_end   = merged_ends[keep],
      stringsAsFactors = FALSE
    )
    df[[group_col]] <- grp
    result_list[[k]] <- df
  }

  out <- do.call(rbind, Filter(Negate(is.null), result_list))
  if (is.null(out)) return(empty)
  rownames(out) <- NULL
  out
}

# Insert NA breaks into a smoothed data frame at each consensus deletion interval.
# Sentinel NA rows are added at del_start - 0.5 and del_end + 0.5 to guarantee
# a visible gap even when no grid point falls inside the deletion.
.insert_deletion_breaks <- function(smoothed, deletion_ranges, group_col) {
  if (nrow(deletion_ranges) == 0L) return(smoothed)

  id_cols <- setdiff(names(smoothed), c("position", "mean_prob", "lower", "upper"))
  sentinel_list <- vector("list", nrow(deletion_ranges))

  for (i in seq_len(nrow(deletion_ranges))) {
    grp_val   <- deletion_ranges[[group_col]][i]
    del_start <- deletion_ranges$del_start[i]
    del_end   <- deletion_ranges$del_end[i]

    # Mask grid points within the deletion
    in_grp <- smoothed[[group_col]] == grp_val
    in_del <- smoothed$position >= del_start & smoothed$position <= del_end
    smoothed$mean_prob[in_grp & in_del] <- NA_real_
    if ("lower" %in% names(smoothed)) smoothed$lower[in_grp & in_del] <- NA_real_
    if ("upper" %in% names(smoothed)) smoothed$upper[in_grp & in_del] <- NA_real_

    # Build sentinel rows — one pair per unique line-identity combo in this group
    grp_rows  <- smoothed[in_grp, id_cols, drop = FALSE]
    templates <- unique(grp_rows)

    if (nrow(templates) == 0L) next

    sentinels <- vector("list", nrow(templates) * 2L)
    for (j in seq_len(nrow(templates))) {
      s1 <- templates[j, , drop = FALSE]
      s1$position  <- del_start - 0.5
      s1$mean_prob <- NA_real_
      # NOTE: `templates` is derived from `id_cols`, which now excludes
      # lower/upper, so s1/s2 never carry those columns yet at this point.
      # Check against `smoothed` (the source of truth for which columns
      # exist) rather than `s1`/`s2`, then add the columns as NA so the
      # final rbind()/column-select against names(smoothed) succeeds.
      if ("lower" %in% names(smoothed)) { s1$lower <- NA_real_; s1$upper <- NA_real_ }
      s2 <- templates[j, , drop = FALSE]
      s2$position  <- del_end + 0.5
      s2$mean_prob <- NA_real_
      if ("lower" %in% names(smoothed)) { s2$lower <- NA_real_; s2$upper <- NA_real_ }
      sentinels[[2L * j - 1L]] <- s1
      sentinels[[2L * j]]      <- s2
    }
    sentinel_list[[i]] <- do.call(rbind, sentinels)
  }

  all_sentinels <- do.call(rbind, Filter(Negate(is.null), sentinel_list))
  if (!is.null(all_sentinels) && nrow(all_sentinels) > 0L) {
    smoothed <- rbind(smoothed, all_sentinels[, names(smoothed), drop = FALSE])
  }

  smoothed <- smoothed[order(smoothed[[group_col]], smoothed$position), , drop = FALSE]
  rownames(smoothed) <- NULL
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

# Add a CI ribbon behind the smooth line(s) when requested and columns exist.
# `fill_aes` controls the visual fill colour (tracking the matching line
# colour aesthetic). `group_aes` controls polygon separation and defaults to
# `fill_aes`; pass an explicit interaction() when a branch has more than one
# identity variable (e.g. group + mod_code) so overlapping CI bands don't
# collapse into a single self-crossing polygon.
.add_ci_ribbon <- function(p, smoothed, show_ci, fill_aes = NULL, group_aes = NULL) {
  if (!isTRUE(show_ci)) return(p)
  if (!all(c("lower", "upper") %in% names(smoothed))) return(p)
  rib <- smoothed[!is.na(smoothed$lower) & !is.na(smoothed$upper), , drop = FALSE]
  if (nrow(rib) == 0L) return(p)
  aes_args <- list(x = quote(.data$position),
                   ymin = quote(.data$lower),
                   ymax = quote(.data$upper))
  if (!is.null(fill_aes)) aes_args$fill <- fill_aes
  if (!is.null(group_aes)) aes_args$group <- group_aes
  ribbon <- ggplot2::geom_ribbon(
    data = rib,
    mapping = do.call(ggplot2::aes, aes_args),
    alpha = 0.2, colour = NA,
    inherit.aes = FALSE,
    show.legend = FALSE
  )
  # Insert ribbon *before* existing line layers so it renders behind them.
  p$layers <- c(list(ribbon), p$layers)
  p
}

# Common ggplot2 layers for the smoothed modification probability panel.
# Returns a list of scales/coords/theme components shared by every branch.
# `y_label` differs between call modes: continuous mode plots a mean probability,
# binary mode plots a fraction of methylated calls (see .smooth_y_label()).
.smooth_panel_base <- function(region_start, region_end,
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
.smooth_y_label = function(call_mode) {
  if (identical(call_mode, "binary")) {
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
.delta_y_label = function(call_mode) {
  if (identical(call_mode, "binary")) {
    "Δ fraction\nmethylated"
  } else {
    "Δ methylation"
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
.sites_for_aggregation = function(sites, call_mode, call_threshold,
                                  call_ambiguous) {
  if (identical(call_mode, "binary")) {
    .binarize_sites(sites, call_threshold, call_ambiguous)
  } else {
    sites
  }
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
# ribbon fill) scale.
.build_smooth_panel = function(smoothed, region_start, region_end, y_label,
                               show_ci, colour = NULL, linetype = NULL,
                               labels = list(), colour_values = NULL,
                               linetype_values = NULL, linetype_name = NULL,
                               fixed_colour = .PROB_GRADIENT$high) {
  # `.data$<col>`, exactly as a literal aes(.data$group) would be written
  col_sym = function(col) rlang::call2("$", quote(.data), rlang::sym(col))
  aes_args = list(x = quote(.data$position), y = quote(.data$mean_prob))
  if (!is.null(colour))   aes_args$colour   = col_sym(colour)
  if (!is.null(linetype)) aes_args$linetype = col_sym(linetype)

  line = if (is.null(colour)) {
    ggplot2::geom_line(linewidth = 1, colour = fixed_colour)
  } else {
    ggplot2::geom_line(linewidth = 1)
  }
  p = ggplot2::ggplot(smoothed, do.call(ggplot2::aes, aes_args)) + line
  if (!is.null(linetype_values)) {
    p = p + ggplot2::scale_linetype_manual(values = linetype_values, name = linetype_name)
  }
  p = p +
    .smooth_panel_base(region_start, region_end, y_label) +
    do.call(ggplot2::labs, c(list(x = "Genomic position (bp)"), labels))
  if (!is.null(colour_values)) {
    p = p + ggplot2::scale_colour_manual(values = colour_values, na.value = "grey50")
  }

  # With two identity variables the ribbon needs both in `group`, or the
  # overlapping bands collapse into one self-crossing polygon.
  group_aes = if (!is.null(colour) && !is.null(linetype)) {
    rlang::expr(interaction(!!col_sym(colour), !!col_sym(linetype)))
  }
  p = .add_ci_ribbon(p, smoothed, show_ci,
                     fill_aes = if (!is.null(colour)) col_sym(colour),
                     group_aes = group_aes)
  if (!is.null(colour_values)) {
    p = p + ggplot2::scale_fill_manual(values = colour_values, na.value = "grey50",
                                       guide = "none")
  }
  p
}

# Smooth panel for a single sample: one line per group (or one overall line),
# with a linetype per modification code when several codes are present.
.single_smooth_panel = function(data, sites_agg, region_start, region_end,
                                opts) {
  codes = unique(data$sites$mod_code)
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

  y_label = .smooth_y_label(opts$call_mode)
  if (grouped) {
    .build_smooth_panel(
      smoothed, region_start, region_end, y_label, opts$show_ci,
      colour = "group", linetype = if (multi_code) "mod_code",
      labels = list(colour = "Group"), colour_values = opts$group_colours,
      linetype_values = if (multi_code) {
        stats::setNames(c("solid", "dashed", "dotdash", "dotted")[seq_along(codes)], codes)
      },
      linetype_name = "Modification"
    )
  } else if (multi_code) {
    .build_smooth_panel(smoothed, region_start, region_end, y_label, opts$show_ci,
                        colour = "mod_code", labels = list(colour = "Modification"))
  } else {
    .build_smooth_panel(smoothed, region_start, region_end, y_label, opts$show_ci)
  }
}

# Shared smooth panel for multi-sample data: one line per sample, or, when any
# sample is grouped, colour by group and linetype by sample.
.multi_smooth_panel = function(samples, prepared, region_start, region_end,
                               opts) {
  smoothed = do.call(rbind, Map(function(s, prep, nm) {
    if (is.null(prep) || nrow(prep$sites_agg) == 0L) return(NULL)
    # One line per group of a grouped sample (unphased reads get their own
    # "NA" line), or a single line for an ungrouped one.
    line_key = function(df) {
      if (is.null(s$group_tag)) rep(nm, nrow(df)) else
        ifelse(is.na(df$group), "NA", as.character(df$group))
    }
    sites = prep$sites_agg
    reads = s$reads
    sites$smooth_key = line_key(sites)
    reads$smooth_key = line_key(reads)

    sm = smooth_methylation(sites, group_col = "smooth_key", span = opts$smooth_span)
    sm = .apply_deletion_breaks(sm, s$cigar_features, reads, "smooth_key",
                                opts$min_indel_size, opts$show_cigar)
    sm$sample = rep(nm, nrow(sm))
    sm$group = if (is.null(s$group_tag)) rep(NA_character_, nrow(sm)) else sm$smooth_key
    sm$smooth_key = NULL
    sm
  }, samples, prepared, names(samples)))

  if (is.null(smoothed) || nrow(smoothed) == 0L) {
    return(.message_plot("No methylation sites", size = 4))
  }
  rownames(smoothed) = NULL

  y_label = .smooth_y_label(opts$call_mode)
  if (any(!vapply(samples, function(s) is.null(s$group_tag), logical(1L)))) {
    .build_smooth_panel(smoothed, region_start, region_end, y_label, opts$show_ci,
                        colour = "group", linetype = "sample",
                        labels = list(colour = "Group", linetype = "Sample"),
                        colour_values = opts$group_colours)
  } else {
    .build_smooth_panel(smoothed, region_start, region_end, y_label, opts$show_ci,
                        colour = "sample", labels = list(colour = "Sample"))
  }
}

# Delta panel for exactly two groups, or NULL (with a message) otherwise.
.delta_panel = function(data, sites_agg, region_start, region_end, opts) {
  delta_df = .compute_group_delta(sites_agg, "group", span = opts$smooth_span)
  if (is.null(delta_df)) return(NULL)
  # Read the group names off the attribute *now*: breaking the line on
  # deletions rebuilds delta_df, which drops attributes.
  delta_groups = attr(delta_df, "groups")
  if (isTRUE(opts$show_cigar)) {
    delta_df = .break_delta_on_deletions(delta_df, data$cigar_features, data$reads,
                                         opts$min_indel_size)
  }
  # Colour the delta area by whichever group is higher, so it reads against
  # the group colours in the panels above. Fall back to the standalone
  # diverging pair only if *both* groups can't be resolved — a partial
  # fallback would leave one half matching and one half not.
  fill_neg = .DELTA_DIVERGING$neg
  fill_pos = .DELTA_DIVERGING$pos
  group_colours = opts$group_colours
  if (!is.null(group_colours) && all(delta_groups %in% names(group_colours))) {
    fill_neg = group_colours[[delta_groups[1L]]]
    fill_pos = group_colours[[delta_groups[2L]]]
  }
  .build_delta_panel(delta_df, region_start, region_end,
                     .delta_y_label(opts$call_mode),
                     fill_pos = fill_pos, fill_neg = fill_neg)
}

#' Plot read-level methylation data
#'
#' Creates a ggplot2 visualisation of read-level base modification data. When
#' no grouping is present, produces a single panel showing reads as grey bars
#' with coloured modification dots. When groups are present, adds a bottom
#' panel with loess-smoothed mean modification probability per group and
#' combines the panels using patchwork. With `call_mode = "binary"` the smooth
#' panel plots the fraction of methylated calls instead of the mean probability.
#'
#' @param data A `methylation_data` object returned by [read_methylation()].
#' @param sort_by Character vector of column names from `data$reads` used
#'   to sort reads before packing into lanes (reads are always packed per
#'   group, so group order is fixed). Default NULL uses `"start"` when
#'   ungrouped or `c("start", "group", "mean_mod_prob")` when grouped.
#'   `mean_mod_prob` is the per-read mean modification probability, or the
#'   per-read fraction of methylated calls when `call_mode = "binary"`.
#' @param colour_low Colour for low modification probability (default
#'   `"#BDBDBD"`).
#' @param colour_high Colour for high modification probability (default
#'   `"#C62828"`).
#' @param colour_ambiguous Colour for ambiguous calls (default `"#78909C"`, a
#'   slate blue-grey deliberately off the `colour_low`/`colour_high` ramp so it
#'   reads as "no confident call" rather than "intermediate methylation"). Only
#'   used when `call_mode = "binary"` and `call_ambiguous` is non-`NULL`.
#' @param line_width Linewidth of modification site lines (default 0.2).
#' @param colour_strand Logical. When `TRUE`, read bars are coloured by strand
#'   (`"+"` = forward, `"-"` = reverse). Ignored when data is grouped; group
#'   colour takes precedence. Default `FALSE`.
#' @param strand_colours Named character vector with `"+"` and `"-"` entries
#'   giving bar colours for each strand. Only used when `colour_strand = TRUE`
#'   and data is ungrouped. Default `c("+" = "#4393C3", "-" = "#D6604D")`.
#' @param group_colours Named character vector of colours per group. Defaults
#'   to `c("1" = "#0072B2", "2" = "#E69F00")` (the Okabe-Ito colorblind-safe
#'   blue/orange pair, see `.GROUP_PALETTE_DEFAULT`), matching the typical HP
#'   haplotype tag output. Pass `NULL` to use ggplot2 defaults, or supply a
#'   fully named vector for other group names.
#' @param smooth_span Loess smoothing span for the bottom panel. When `NULL`
#'   (default), an adaptive span is computed per group as
#'   `max(0.15, min(0.75, 15 / n_unique_sites))`, targeting ~15 data points
#'   per local fit regardless of region size or CpG density. Pass an explicit
#'   numeric value (e.g. `0.3`) to use a fixed span.
#' @param panel_heights Numeric vector giving relative heights of the panels.
#'   When `annotations = NULL`, expects length 2 (`reads`, `smooth`); when
#'   annotations are provided, expects length 3 (`gene track`, `reads`,
#'   `smooth`) — the gene panel is placed at the top, reads in the middle, and
#'   the smoothed probability curve at the bottom. When `show_delta = TRUE`
#'   and grouping yields exactly two groups, an additional delta panel is
#'   appended last, adding one to the expected length in each case above.
#'   Pass `NULL` to use the built-in defaults.
#' @param annotations A `gene_annotations` object returned by
#'   [read_annotations()], or `NULL` (default). When provided, a gene
#'   annotation track is placed at the top of the composite figure (above the
#'   read panel and the smooth panel).
#' @param variants A `variant_data` object returned by [read_variants()], or
#'   `NULL` (default). When provided, per-read base letters are drawn at
#'   variant positions (coloured by match to ref/alt), and vertical dashed
#'   lines are added at every variant position across all panels. Requires
#'   that `data` was produced by a current version of [read_methylation()] that
#'   stores sequences and CIGARs in the object.
#' @param show_cigar Logical. When `TRUE` (default), structural variants from
#'   CIGAR strings (insertions, deletions) are displayed on reads. Insertions
#'   appear as purple I-beam markers spanning the full height of the read bar;
#'   deletions open a gap in the read bar and draw a thin black line. Set to
#'   `FALSE` to hide these features.
#' @param min_indel_size Integer. Minimum size (in bp) for insertions and
#'   deletions to be displayed when `show_cigar = TRUE`. Indels smaller than
#'   this threshold are hidden to reduce visual clutter from common small
#'   indels. Default `50`.
#' @param show_supplementary Logical. When `TRUE` (default), a coloured halo
#'   is drawn around read bars indicating the chromosome of the supplementary
#'   alignment partner (from the SA BAM tag). The original bar colouring
#'   (group, strand, or default grey) is preserved inside the halo. Reads with
#'   no supplementary alignment have no halo.
#' @param bnd_match_tol Integer. Position tolerance (bp) for matching
#'   supplementary-alignment breakpoints to VCF BND calls. The SA matching runs
#'   whenever `variants` is supplied and reads carry SA tags; the visual border
#'   marking requires `show_supplementary = TRUE`. Default 50.
#' @param show_ci Logical. When `TRUE` (default), a shaded ribbon showing the
#'   loess confidence interval (`lower`/`upper` from [smooth_methylation()])
#'   is drawn behind each smooth line in the bottom panel. Has no effect when
#'   fewer than 4 unique positions are available for a group/code (no CI is
#'   computed in that case). Set to `FALSE` to hide the ribbon.
#' @param call_mode Character. `"continuous"` (default) colours modification
#'   sites by a continuous probability gradient (`colour_low` to
#'   `colour_high`), and the smooth panel plots the mean modification
#'   probability per position. `"binary"` classifies each site as methylated,
#'   unmethylated, or ambiguous (see `call_threshold`/`call_ambiguous`) and
#'   colours them with a discrete scale instead.
#'
#'   `"binary"` also changes what every aggregate reports: the smooth panel, the
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
#'   contribute coarser values than in continuous mode, so the binary curve is
#'   noisier where few reads overlap. The loess span still averages over
#'   neighbouring positions and the confidence ribbon (`show_ci`) widens where
#'   evidence is sparse.
#' @param call_threshold Numeric in `[0, 1]`. Modification probability at or
#'   above which a site is classified as methylated when
#'   `call_mode = "binary"`. Used for both the read panel colours and the
#'   binarised aggregates, so the two panels always split at the same value.
#'   Default `0.5`. Ignored when `call_mode = "continuous"`.
#' @param call_ambiguous `NULL` (default) for a hard threshold, or a numeric
#'   half-width defining a band `[call_threshold - w, call_threshold + w)`
#'   around `call_threshold` within which sites are labelled "ambiguous"
#'   rather than methylated/unmethylated. Ambiguous sites are excluded from
#'   both the numerator and the denominator of the methylated fraction, so the
#'   smooth panel reports the fraction methylated *among confident calls*. The
#'   default `NULL` discards nothing — every site is called at
#'   `call_threshold`. Ignored when `call_mode = "continuous"`.
#' @param show_delta Logical. When `TRUE`, and grouping yields exactly two
#'   groups, an additional panel is appended below the smooth panel showing
#'   the signed difference in loess-smoothed modification probability between
#'   the two groups (group2 - group1), coloured by sign. Requires `data` to be
#'   grouped (`data$group_tag` non-`NULL`) with exactly two distinct group
#'   values; otherwise the panel is silently skipped (with a message). Not
#'   supported for multi-sample data (`multi_methylation_data`); the argument
#'   is accepted there but ignored (with a message). Default `FALSE`.
#'
#' @return A [ggplot2::ggplot] object (ungrouped) or a
#'   [patchwork::patchwork] composite (grouped).
#'
#' @examples
#' \dontrun{
#' md <- read_methylation("sample.bam", "chr1:1000-2000")
#' plot_methylation(md)
#'
#' md <- read_methylation("sample.bam", "chr1:1000-2000", group_tag = "HP")
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
                            bnd_match_tol = 50L,
                            show_ci = TRUE,
                            call_mode = c("continuous", "binary"),
                            call_threshold = 0.5,
                            call_ambiguous = NULL,
                            show_delta = FALSE) {
  call_mode = match.arg(call_mode)

  # --- 1. Validate input ---
  multi = inherits(data, "multi_methylation_data")
  if (!multi && !inherits(data, "methylation_data")) {
    stop(
      "`data` must be a `methylation_data` or `multi_methylation_data` object.",
      call. = FALSE
    )
  }
  if (!is.logical(colour_strand) || length(colour_strand) != 1L) {
    stop("'colour_strand' must be a single logical value.", call. = FALSE)
  }
  if (isTRUE(colour_strand) && !all(c("+", "-") %in% names(strand_colours))) {
    stop("'strand_colours' must be a named vector with '+' and '-' entries.", call. = FALSE)
  }
  if (!is.null(variants) && !inherits(variants, "variant_data")) {
    stop("`variants` must be a `variant_data` object returned by read_variants().",
         call. = FALSE)
  }
  if (!is.null(annotations) && !inherits(annotations, "gene_annotations")) {
    stop("`annotations` must be a `gene_annotations` object from read_annotations().",
         call. = FALSE)
  }

  samples = if (multi) data$samples else list(data)
  if (isTRUE(colour_strand) &&
      any(!vapply(samples, function(s) is.null(s$group_tag), logical(1L)))) {
    message("'colour_strand = TRUE' is ignored when reads are grouped; group colour takes precedence.")
  }
  if (multi && isTRUE(show_delta)) {
    message("`show_delta` is not supported for multi-sample data; ignoring.")
  }
  if (!multi && nrow(data$reads) == 0L) {
    return(.message_plot("No reads in region", size = 5))
  }

  region_start = GenomicRanges::start(data$region)
  region_end   = GenomicRanges::end(data$region)
  opts = list(smooth_span = smooth_span, min_indel_size = min_indel_size,
              show_cigar = show_cigar, show_ci = show_ci, call_mode = call_mode,
              group_colours = group_colours)

  # --- 2. Sort and pack each sample's reads ---
  prepared = lapply(samples, function(s) {
    if (nrow(s$reads) == 0L) return(NULL)
    sites_agg = .sites_for_aggregation(s$sites, call_mode, call_threshold, call_ambiguous)
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
      group_colours      = group_colours,
      show_x_axis        = FALSE,
      variant_overlay    = build_variant_overlay(s, variants, bnd_match_tol = bnd_match_tol),
      show_cigar         = show_cigar,
      cigar_features     = if (isTRUE(show_cigar)) s$cigar_features else NULL,
      min_indel_size     = min_indel_size,
      show_supplementary = show_supplementary,
      call_mode          = call_mode,
      call_threshold     = call_threshold,
      call_ambiguous     = call_ambiguous
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

  patchwork::wrap_plots(panels, ncol = 1, heights = heights)
}
