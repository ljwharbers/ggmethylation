# ggmethylation — Design Spec

**Date:** 2026-03-27
**Status:** Approved for planning

---

## Overview

An R package for visualising read-level DNA base modifications from long-read sequencing data (Oxford Nanopore / PacBio). Given a BAM file and a genomic region, the package produces a two-panel ggplot2 composite:

- **Top panel:** individual reads with modification probability shown as coloured dots at each site
- **Bottom panel:** smoothed mean modification rate per group across the region (only when grouping is defined)

---

## Input

| Argument | Type | Description |
|----------|------|-------------|
| `bam` | `character` | Path to BAM file |
| `region` | `character` | Genomic region string, e.g. `"chr1:1000-2000"` |
| `mod_code` | `character` | MM tag modification code, e.g. `"m"` (5mC), `"h"` (5hmC), `"a"` (6mA) |
| `group_tag` | `character` or `NULL` | BAM tag to use for grouping reads, e.g. `"HP"`, `"RG"`. `NULL` disables grouping and hides the bottom panel. |
| `max_reads` | `integer` | Maximum number of reads to display. Default: `200`. Reads are randomly downsampled when the count exceeds this value. |

The BAI index file must exist alongside the BAM file (same path with `.bai` extension). The package will error informatively if the index is missing.

---

## Architecture: Two-Layer (Parse → Plot)

### Layer 1 — `read_methylation()`

Loads and parses only the reads overlapping the specified region using `Rsamtools::scanBam()` with a `ScanBamParam` `which` argument. The BAI index is used to jump directly to relevant reads — the full BAM is never loaded.

**Returns** a named list (S3 class `methylation_data`) with:

```
$reads        data.frame — one row per read
              columns: read_name, start, end, strand, group (from BAM tag)

$sites        data.frame — one row per modification site per read
              columns: read_name, position, mod_prob (0–1, from ML tag), mod_code

$region       GRanges or named vector — the parsed genomic region
$mod_code     character — the modification code used
$group_tag    character or NULL
```

**MM/ML tag parsing:**

- `MM` tag encodes modification positions as delta-encoded offsets from the relevant base (e.g. C for 5mC)
- `ML` tag encodes corresponding probabilities as uint8 values (0–255, scaled to 0–1 by dividing by 255)
- Only sites matching `mod_code` are extracted
- Rsamtools exposes raw tags; parsing will be implemented in pure R

**Downsampling:** Applied after loading, before returning. If `nrow(reads) > max_reads`, a random sample of `max_reads` rows is drawn (set seed via `set.seed()` before calling for reproducibility). When `group_tag` is set, downsampling is applied globally (not per group) at this stage.

---

### Layer 2 — `plot_methylation()`

Takes a `methylation_data` object and returns a `patchwork` composite of ggplot2 panels.

```r
meth <- read_methylation("sample.bam", region = "chr1:1000-2000",
                         mod_code = "m", group_tag = "HP")
plot_methylation(meth)
```

---

## Top Panel — Read-Level View

- Each read is drawn as a horizontal rectangle spanning its clipped start–end coordinates
- Reads extending beyond the region boundaries are **clipped** to the region
- Modification sites are drawn as filled circles on top of the read bar, coloured by `mod_prob` (0–1)
- Default colour scale: low probability → light grey (`#BDBDBD`), high probability → deep red (`#C62828`), using a continuous gradient (e.g. `scale_colour_gradient()`)
- The colour scale and dot size should be user-overridable via arguments

**Read layout (vertical stacking):**

Reads are packed into lanes to minimise vertical space (genome-browser style), then sorted. Default sort priority (applied before packing):

1. `start` position
2. `group` (if `group_tag` is set)
3. Mean modification probability (descending)

When the user provides `sort_by`, it takes priority over all defaults.

The user can override sort order via a `sort_by` argument accepting a character vector of column names from `$reads`, e.g. `sort_by = c("group", "start")`.

**Axes:**
- x-axis: genomic position (bp)
- y-axis: read lane index (no label, ticks hidden)
- No legend for y-axis

---

## Bottom Panel — Group Summary

Only rendered when `group_tag` is not `NULL`.

- One loess-smoothed line per group across the genomic region
- y-axis: mean modification probability (0–1)
- x-axis: genomic position (shared with top panel, aligned via patchwork)
- Each group gets a distinct colour; palette is user-overridable
- Smoothing span is user-overridable (default `0.3`)
- The raw per-site group means are computed from `$sites`, then smoothed using `stats::loess()`

---

## Output

`plot_methylation()` returns a `patchwork` object. Users can apply ggplot2 themes, modify individual panels, or save with `ggplot2::ggsave()`.

When `group_tag` is `NULL`, only the top panel is returned (a plain `ggplot` object, not patchwork).

---

## Package Structure

```
ggmethylation/
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── read_methylation.R      # Layer 1: BAM parsing, MM/ML tag extraction
│   ├── plot_methylation.R      # Layer 2: top + bottom panel assembly
│   ├── parse_mm_ml.R           # Internal: MM/ML tag parsing logic
│   ├── pack_reads.R            # Internal: read lane packing algorithm
│   ├── smooth_methylation.R    # Internal: per-group loess smoothing
│   └── utils.R                 # Internal: region string parsing, validation
├── man/
├── tests/
│   └── testthat/
└── docs/
```

---

## Dependencies

| Package | Role |
|---------|------|
| `Rsamtools` (Bioconductor) | BAM/BAI reading, region-specific read fetching |
| `GenomicRanges` (Bioconductor) | Region representation |
| `ggplot2` | Plot panels |
| `patchwork` | Compositing top + bottom panels |
| `stats` | `loess()` for smoothing |

All Bioconductor dependencies listed in `DESCRIPTION` under `Imports` with `BiocManager` noted in installation instructions.

---

## Error Handling

- Missing BAI index → informative error with suggested `samtools index` command
- Invalid region string → error with expected format shown
- No reads in region → warning, return empty plot with message annotation
- `group_tag` not found in any read → warning, fall back to `group_tag = NULL` behaviour
- `mod_code` not present in any MM tag → informative error

---

## Out of Scope (for now)

- Multi-sample / multi-BAM overlays
- Bisulfite sequencing BAMs (XM/XB tags)
- Interactive plots (plotly / shiny)
- Saving parsed data to disk (export functions)
- Strand-specific visualisation (reads not coloured by strand)
