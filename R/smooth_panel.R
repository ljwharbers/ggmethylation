# Smooth and delta panels of plot_methylation(): line identities, axis labels, the shared panel builder and the single-/multi-sample/delta panels.

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
# The delta sign is a plotmath `Delta` (drawn from the symbol font) rather
# than the character U+0394, which the default pdf() device cannot draw.
.delta_y_label = function(binary) {
  if (binary) {
    expression(atop(Delta ~ "fraction", "methylated"))
  } else {
    expression(Delta ~ "methylation")
  }
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
