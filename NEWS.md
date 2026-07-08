# ggmethylation (development version)

## New features
* `theme_ggmethylation()` exposes the shared panel theme.
* `plot_methylation()` gains `show_ci` (default `TRUE`) to draw a shaded loess
  confidence-interval ribbon behind the smoothed modification probability
  line(s) in the bottom panel, using the `lower`/`upper` columns from
  `smooth_methylation()`.
* `plot_methylation()` gains `call_mode` (`"continuous"` default or
  `"binary"`), `call_threshold` (default `0.5`), and `call_ambiguous`
  (default `NULL`) to render read-panel modification sites as discrete
  methylated/unmethylated/ambiguous calls instead of a continuous
  probability gradient.
* `plot_methylation()` gains `show_delta` (default `FALSE`). When grouping
  yields exactly two groups, setting `show_delta = TRUE` appends a bottom
  panel showing the signed difference in loess-smoothed modification
  probability between the two groups (group2 - group1), rendered as a
  diverging area coloured by sign. Not supported for multi-sample data.

## Bug fixes
* `.insert_deletion_breaks()` (internal) now also nulls the `lower`/`upper`
  confidence-interval columns inside consensus deletion gaps, matching the
  existing `mean_prob` masking, so the CI ribbon does not bridge across
  deletion breaks.

## Changes
* Default `group_colours` updated to a colorblind-safe Okabe-Ito blue/orange
  pair. Pass an explicit `group_colours` vector to restore prior colours.
