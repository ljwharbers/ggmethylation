# ggmethylation (development version)

## New features
* `theme_ggmethylation()` exposes the shared panel theme.
* `plot_methylation()` gains `show_ci` (default `TRUE`) to draw a shaded loess
  confidence-interval ribbon behind the smoothed modification probability
  line(s) in the bottom panel, using the `lower`/`upper` columns from
  `smooth_methylation()`.

## Bug fixes
* `.insert_deletion_breaks()` (internal) now also nulls the `lower`/`upper`
  confidence-interval columns inside consensus deletion gaps, matching the
  existing `mean_prob` masking, so the CI ribbon does not bridge across
  deletion breaks.

## Changes
* Default `group_colours` updated to a colorblind-safe Okabe-Ito blue/orange
  pair. Pass an explicit `group_colours` vector to restore prior colours.
