# Visualization Improvements — Design

**Date:** 2026-07-06
**Status:** Approved (design)
**Author:** Luuk Harbers (with Claude)

## Context

`ggmethylation` visualises read-level base modification data from long-read
sequencing. The core plot (`plot_methylation()`) already supports: read bars
with per-site probability colouring, strand/group/probability colouring,
loess-smoothed group summary panel, gene annotation tracks, variant overlays
(SNV/SV/BND), CIGAR indels, supplementary-alignment halos, and multi-sample
layouts.

This work adds analytical power and reduces exploratory friction. The primary
consumer is **exploratory analysis** by the maintainer and collaborators, not
publication finish. Improvements are therefore weighted toward reading data
faster and comparing groups directly.

## Scope

Two milestones. Each feature is independently shippable and gated behind a new
argument that defaults to preserving current behaviour.

### Milestone 1 (this cycle)

- **A1 — Group-difference (delta) track**
- **A2 — Confidence ribbon on the smooth**
- **B1 — Binarized call mode**
- **C2 — Centralized theme + colorblind-safe defaults**

### Milestone 2 (later)

- **C1 — Rasterized dense layers**
- **A3 — Per-CpG fraction-methylated track**
- **B2 — Per-read summary sidebar**

Out of scope entirely for now: read × position heatmap, methylation-pattern
read clustering, interactivity (ggiraph/plotly). These may be revisited.

## Feature designs

### A1 — Group-difference (delta) track

**Goal.** When data is grouped, show *where* groups differ instead of making the
eye diff two loess lines.

- New arg `show_delta = FALSE` on `plot_methylation()`.
- Active only when the data is grouped into **exactly two** groups. With a
  different number of groups (including ungrouped, or >2), emit a `message()`
  and skip the panel silently otherwise.
- Computation: fit both group smooths on a **shared** position grid (currently
  each group's loess is predicted on its own min–max grid; the delta needs a
  common grid spanning the intersection or union of positions — use the union,
  predicting each loess over the shared grid, `NA` outside a group's support).
  Delta = `group2_value − group1_value` at each grid point. Group ordering
  follows `.ordered_plot_groups()` so the sign is deterministic.
- Rendering: a `geom_area`/`geom_ribbon` around a zero baseline with a diverging
  fill (positive vs negative), plus a horizontal zero line. Diverging palette is
  colorblind-safe (see C2).
- Layout: appears as an additional bottom sub-panel **below** the smooth panel.
  `panel_heights` handling must extend to the new optional panel count. Default
  delta panel height ~0.2 (relative to reads = 1).
- Interaction with deletion breaks: reuse the same shared grid, and apply the
  existing `.apply_deletion_breaks()` logic so the delta line also breaks over
  consensus deletions.

### A2 — Confidence ribbon on the smooth

**Goal.** Show where the smoothed curve is trustworthy vs sparse.

- New arg `show_ci = TRUE` (ribbon on by default — it is low-risk and
  informative; can be set FALSE to restore the bare line).
- `smooth_methylation()` gains SE output: `stats::predict(fit, ..., se = TRUE)`
  returns `$fit` and `$se.fit`. Add `lower` and `upper` columns computed as
  `fit ± 1.96 * se.fit`, clamped to `[0, 1]`. When the raw-means fallback path
  is taken (fewer than 4 unique positions), `lower`/`upper` are `NA`.
- Rendering: a `geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2)` drawn
  **behind** each existing `geom_line`, matching the line's colour (fill by
  group / mod_code as applicable). Ribbon rows with `NA` bounds are dropped.
- Deletion breaks: sentinel-NA insertion must also null out `lower`/`upper` so
  the ribbon breaks with the line.

### B1 — Binarized call mode

**Goal.** Reduce visual noise; show methylation *calls* rather than continuous
scores.

- New args `call_mode = c("continuous", "binary")` (default `"continuous"`,
  matched via `match.arg`) and `call_threshold = 0.5`.
- In binary mode, per-site read-panel segments are coloured by a two-state
  discrete scale: below threshold → `colour_low`, at/above → `colour_high`.
  Reuses the existing endpoint colour args so the palette stays consistent.
- Optional ambiguous band: `call_ambiguous = NULL` (default off). When set to a
  numeric half-width `w`, calls with `abs(mod_prob - call_threshold) < w` are
  drawn in a muted grey and excluded from any downstream binary aggregation.
- Scope: affects the **read panel** only in v1. The smooth panel continues to
  show mean probability. A later option (documented as future work) is a
  fraction-of-binary-calls smooth.
- Implementation lands in `.add_mod_prob_segments()` (or a sibling helper),
  switching the colour aesthetic/scale based on `call_mode`. The threshold and
  mode are threaded from `plot_methylation()` → `build_read_panel()` →
  `.add_mod_prob_segments()`.

### C2 — Centralized theme + colorblind-safe defaults

**Goal.** One place to tune look-and-feel; trustworthy default colours.

- New exported `theme_ggmethylation()` that returns the shared
  `theme_minimal()` + `theme(...)` block currently duplicated across
  `build_read_panel()` and `.smooth_panel_base()`. Both call sites refactor to
  use it. The read panel adds its y-axis-blanking `theme()` on top.
- Centralize palette constants in one internal file (e.g. `R/palettes.R`):
  - Group default palette → Okabe–Ito (`group_colours` default updated; keep it
    named so existing `c("1", "2")` haplotype output still maps). Verify the
    current teal/orange default is retained or replaced deliberately.
  - Probability gradient → keep the colorblind-safe grey→red
    (`colour_low`/`colour_high` defaults unchanged).
  - Diverging palette for the delta track (A1) defined here.
- No change to public defaults that would alter existing figures unless the
  change is a deliberate colorblind-safety upgrade, called out in NEWS.

## Data structures

No change to the `methylation_data` S3 contract. `smooth_methylation()` return
gains `lower`/`upper` columns (additive; downstream code selects columns
explicitly so this is backward-compatible). All new behaviour is argument-gated.

## Testing

- **A1:** unit test the shared-grid delta computation (deterministic sign,
  correct handling of non-overlapping support, message + skip for ≠2 groups).
  Snapshot/structure test that the composite gains one panel when `show_delta`.
- **A2:** test that `smooth_methylation()` returns `lower`/`upper`, clamped to
  `[0,1]`, `NA` on the raw-means fallback. Test ribbon layer presence.
- **B1:** test binary colour assignment at/above/below threshold and the
  ambiguous band exclusion; `match.arg` validation.
- **C2:** test `theme_ggmethylation()` returns a theme; test palette constants
  exist and are the expected length/names. Refactor must not change existing
  snapshot tests.
- Integration tests continue to skip gracefully when the BAM fixture is absent.

## Risks & mitigations

- **Layout regressions** from the new delta panel: extend `panel_heights`
  validation and defaults carefully; cover with a panel-count test.
- **Shared-grid loess** for A1 is the trickiest piece — predicting each group's
  loess outside its support yields `NA`; ensure the delta is only drawn where
  both groups have support.
- **Default palette change (C2)** could alter existing users' figures — only
  change defaults for a clear colorblind-safety reason and document in NEWS.
