# Okabe-Ito colorblind-safe qualitative palette
.OKABE_ITO <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)

# Default group palette. Named "1"/"2" to match HP haplotype tag output.
# Uses two well-separated Okabe-Ito hues (blue / orange).
.GROUP_PALETTE_DEFAULT <- c("1" = "#0072B2", "2" = "#E69F00")

# Map a group palette onto the group values actually present.
#
# `.GROUP_PALETTE_DEFAULT` is named "1"/"2" because those are the values an HP
# haplotype tag produces. Any other grouping -- SNV genotype ("REF"/"ALT"),
# a custom BAM tag, sample names -- shares no names with it, and
# scale_*_manual() then matches nothing: every group falls through to
# `na.value` grey and ggplot2 emits "No shared levels found between
# `names(values)` of the manual scale and the data's values".
#
# So resolve the palette against the real groups: keep an exact match, and
# otherwise assign the colours positionally over `.ordered_plot_groups()`
# order (which is also the order the delta track's sign convention uses).
# Warn only when the caller passed their own named palette that does not fit --
# the package default silently adapting is the intended behaviour.
.resolve_group_colours <- function(group_colours, groups) {
  if (is.null(group_colours)) return(NULL)

  levels <- .ordered_plot_groups(groups)
  levels <- levels[!is.na(levels)]
  if (length(levels) == 0L) return(group_colours)

  nms <- names(group_colours)
  if (!is.null(nms) && all(levels %in% nms)) {
    return(group_colours)
  }

  if (!is.null(nms) && !identical(group_colours, .GROUP_PALETTE_DEFAULT)) {
    warning(
      "`group_colours` names (", paste(nms, collapse = ", "),
      ") do not cover the groups present (", paste(levels, collapse = ", "),
      "); filling in the rest.",
      call. = FALSE
    )
  }

  out <- stats::setNames(rep(NA_character_, length(levels)), levels)

  # 1. Honour any exact name matches.
  matched <- character(0L)
  if (!is.null(nms)) {
    matched   <- levels[levels %in% nms]
    out[matched] <- unname(group_colours[matched])
  }

  # 2. Fill the rest from the palette entries that went unused, then from
  #    Okabe-Ito. Never reuse a colour already assigned: two groups sharing a
  #    colour is worse than falling off the requested palette, because the
  #    figure then cannot be read at all.
  needed <- names(out)[is.na(out)]
  if (length(needed) > 0L) {
    spare <- if (is.null(nms)) {
      unname(group_colours)
    } else {
      unname(group_colours[!(nms %in% matched)])
    }
    pool <- c(spare, .OKABE_ITO)
    pool <- pool[!duplicated(pool)]
    pool <- setdiff(pool, out[!is.na(out)])
    out[needed] <- rep(pool, length.out = length(needed))[seq_along(needed)]
  }

  out
}

# Probability gradient endpoints (colorblind-safe grey -> red)
.PROB_GRADIENT <- list(low = "#BDBDBD", high = "#C62828")

# Fallback diverging palette for the delta track. Used when `group_colours` is
# NULL or does not name both groups; otherwise the delta panel takes its fills
# straight from the group palette (neg = group 1, pos = group 2).
.DELTA_DIVERGING <- list(neg = "#0072B2", pos = "#D55E00", zero = "grey60")

# Colour for low-confidence ("ambiguous") calls in binary call mode. Chosen to
# sit off the colour_low -> colour_high ramp so it reads as "no confident call"
# rather than "intermediate methylation", and desaturated enough not to be
# mistaken for a group colour on the read bars beneath.
.CALL_AMBIGUOUS_DEFAULT <- "#78909C"

#' ggmethylation plot theme
#'
#' Shared minimal theme used across `ggmethylation` panels. Provides the base
#' `theme_minimal()` plus compact legend styling and suppressed minor grid.
#'
#' @return A [ggplot2::theme] object.
#' @export
theme_ggmethylation <- function() {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.text      = ggplot2::element_text(size = ggplot2::rel(0.75)),
      legend.title     = ggplot2::element_text(size = ggplot2::rel(0.75)),
      legend.key.size  = ggplot2::unit(0.4, "cm")
    )
}
