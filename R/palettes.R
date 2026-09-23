# Default group palette. Named "1"/"2" to match HP haplotype tag output.
# Uses two well-separated hues (blue / orange) from the colorblind-safe
# Okabe-Ito palette.
.GROUP_PALETTE_DEFAULT <- c("1" = "#0072B2", "2" = "#E69F00")

# Probability gradient endpoints (colorblind-safe grey -> red). Default for
# `colour_low`/`colour_high` in plot_methylation() and plot_insertion_locus().
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
