# ggmethylation

Read-level base modification visualisation for long-read sequencing data (Oxford Nanopore / PacBio).

Given a modBAM file and a genomic region, `ggmethylation` produces a two-panel plot:

- **Top panel:** individual reads as horizontal bars, with coloured dots at each modification site (colour encodes modification probability from the ML tag)
- **Bottom panel:** loess-smoothed mean modification probability per group (only shown when grouping is specified)

## Installation

`ggmethylation` depends on the following (Bioconductor) packages. Install them first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Rsamtools", "GenomicRanges", "GenomicAlignments", "IRanges"))
install.packages(c("ggplot2", "patchwork"))
```

Then install `ggmethylation`:

```r
# From GitHub (once published):
devtools::install_github("ljwharbers/ggmethylation")

# Or from a local clone:
# devtools::install("path/to/ggmethylation")
```

## Usage

```r
library(ggmethylation)

# Parse a modBAM file (only loads reads in the specified region)
meth <- read_methylation(
  bam = "sample.bam",
  region = "chr1:1000-2000",
  mod_code = "m",         # 5mC (default)
  group_tag = "HP",       # group by haplotype tag
  max_reads = 200         # downsample if more reads overlap
)

# Plot
plot_methylation(meth)
```

### Without grouping

When `group_tag = NULL` (the default), only the read-level panel is shown:

```r
meth <- read_methylation("sample.bam", "chr1:1000-2000")
plot_methylation(meth)
```

### Customisation

`plot_methylation()` returns a ggplot2 / patchwork object, so you can theme and modify it:

```r
library(ggplot2)

plot_methylation(
  meth,
  colour_low = "#EEEEEE",
  colour_high = "#1B5E20",
  dot_size = 2,
  group_colours = c("1" = "steelblue", "2" = "coral"),
  smooth_span = 0.5,
  panel_heights = c(4, 1)
) + theme(text = element_text(size = 14))
```

### Sorting

Reads are sorted before packing into lanes. The default sort order is:

1. Start position
2. Group (if grouping is active)
3. Mean modification probability

Override with `sort_by`:

```r
plot_methylation(meth, sort_by = c("group", "start"))
```

## Supported modification types

Any base modification encoded in MM/ML tags (SAM spec). Common codes:

| Code | Modification |
|------|-------------|
| `m`  | 5-methylcytosine (5mC) |
| `h`  | 5-hydroxymethylcytosine (5hmC) |
| `a`  | 6-methyladenine (6mA) |

Specify via the `mod_code` argument.

## Requirements

- A BAM file with MM/ML base modification tags (modBAM format)
- A corresponding BAM index file (`.bai`)
- R >= 4.0

## License

MIT
