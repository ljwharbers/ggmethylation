# Default group palette. Named "1"/"2" to match HP haplotype tag output.
# Uses two well-separated hues (blue / orange) from the colorblind-safe
# Okabe-Ito palette.
.GROUP_PALETTE_DEFAULT = c("1" = "#0072B2", "2" = "#E69F00")

# Colour-blind-safe (Okabe-Ito) sequence used, in order, for groups that the
# palette in use does not name (SNV groups "REF"/"ALT", a third haplotype, ...).
# Starts with the two default colours so an unnamed two-group object looks the
# same as HP 1/2.
.GROUP_PALETTE_EXTENSION = c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
                             "#56B4E9", "#D55E00", "#F0E442", "#000000")

# Colour for reads/lines whose group is NA (unphased reads, ungrouped samples
# in a mixed multi-sample plot).
.GROUP_NA_COLOUR = "grey50"

# Sorted, de-duplicated non-NA group levels (the order every panel uses).
.group_levels = function(groups) {
  groups = as.character(groups)
  sort(unique(groups[!is.na(groups)]))
}

# Resolve `group_colours` (as passed to plot_methylation()) into a named
# colour vector covering exactly the observed group levels, in level order.
# This is the single source of group colours for the read, smooth and delta
# panels.
#
# * NULL: ggplot2's default hue palette over the levels.
# * Unnamed: matched to the sorted levels by position; too few colours is an
#   error.
# * Named: named levels keep their colour; levels it does not name get the
#   next unused colours of .GROUP_PALETTE_EXTENSION (then hue). A message is
#   emitted when a user-supplied palette misses a level (not for the default
#   palette, which only names HP "1"/"2").
.resolve_group_colours = function(group_colours, groups) {
  levels = .group_levels(groups)
  if (length(levels) == 0L) return(stats::setNames(character(0L), character(0L)))
  if (is.null(group_colours)) {
    return(stats::setNames(scales::hue_pal()(length(levels)), levels))
  }

  nms = names(group_colours)
  if (is.null(nms)) {
    if (length(group_colours) < length(levels)) {
      stop(sprintf(
        "`group_colours` has %d unnamed colour(s) but the data has %d groups (%s). %s",
        length(group_colours), length(levels), paste(levels, collapse = ", "),
        "Supply one colour per group, or name the colours by group."
      ), call. = FALSE)
    }
    return(stats::setNames(unname(group_colours[seq_along(levels)]), levels))
  }

  out = group_colours[intersect(levels, nms)]
  missing = setdiff(levels, nms)
  if (length(missing) > 0L) {
    if (!identical(group_colours, .GROUP_PALETTE_DEFAULT)) {
      message(sprintf(
        "`group_colours` does not name group(s) %s; using generated colours for them.",
        paste(missing, collapse = ", ")
      ))
    }
    pool = setdiff(.GROUP_PALETTE_EXTENSION, out)
    if (length(pool) < length(missing)) {
      pool = c(pool, scales::hue_pal()(length(missing) - length(pool)))
    }
    out = c(out, stats::setNames(pool[seq_along(missing)], missing))
  }
  out[levels]
}

# Line types for a discrete variable (mod codes, samples). The first six are
# the named, easily told apart types; the rest are ggplot2's default sequence.
.LINETYPES = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash",
               "1343", "73", "2262", "12223242", "F282", "F4448444", "224282F2")

.linetype_values = function(levels, what) {
  if (length(levels) > length(.LINETYPES)) {
    stop(sprintf(
      "Too many %s (%d) to tell apart by line type in the smooth panel (at most %d).",
      what, length(levels), length(.LINETYPES)
    ), call. = FALSE)
  }
  stats::setNames(.LINETYPES[seq_along(levels)], levels)
}

# Probability gradient endpoints (colorblind-safe grey -> red). Default for
# `colour_low`/`colour_high` in plot_methylation() and plot_insertion_locus().
.PROB_GRADIENT = list(low = "#BDBDBD", high = "#C62828")

# Default diverging palette for .build_delta_panel() when called on its own.
# plot_methylation() always passes the two resolved group colours instead
# (neg = group 1, pos = group 2), so the delta matches the panels above.
.DELTA_DIVERGING = list(neg = "#0072B2", pos = "#D55E00", zero = "grey60")

# Colour for low-confidence ("ambiguous") calls in binary call mode. Chosen to
# sit off the colour_low -> colour_high ramp so it reads as "no confident call"
# rather than "intermediate methylation", and desaturated enough not to be
# mistaken for a group colour on the read bars beneath.
.CALL_AMBIGUOUS_DEFAULT = "#78909C"

#' ggmethylation plot theme
#'
#' Shared minimal theme used across `ggmethylation` panels. Provides the base
#' `theme_minimal()` plus compact legend styling and suppressed minor grid.
#'
#' @return A [ggplot2::theme] object.
#' @export
theme_ggmethylation = function() {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.text      = ggplot2::element_text(size = ggplot2::rel(0.75)),
      legend.title     = ggplot2::element_text(size = ggplot2::rel(0.75)),
      legend.key.size  = ggplot2::unit(0.4, "cm")
    )
}
