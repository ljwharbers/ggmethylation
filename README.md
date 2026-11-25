# ggmethylation

Read-level base modification visualisation for long-read sequencing data (Oxford Nanopore / PacBio).

Given a modBAM file (MM/ML tags) and a genomic region, `ggmethylation` returns a
`patchwork` figure with, from top to bottom:

- an optional **gene track** (`annotations`);
- the **read panel**: individual reads as horizontal bars, with a coloured mark
  at each modification site (colour encodes the modification probability from
  the ML tag);
- the **smooth panel**: loess-smoothed mean modification probability, one line
  per group when the reads are grouped and one overall line otherwise;
- an optional **delta panel** (`show_delta = TRUE`) with the difference between
  two groups.

## Installation

`ggmethylation` needs R >= 4.0 and a few Bioconductor packages. Install them first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("Rsamtools", "GenomicRanges", "IRanges"))
install.packages(c("ggplot2", "ggnewscale", "patchwork", "rlang", "scales"))

# Optional: VCF overlay (read_variants()) and gene track (read_annotations())
BiocManager::install(c("VariantAnnotation", "SummarizedExperiment",
                       "GenomicFeatures", "txdbmaker", "rtracklayer",
                       "AnnotationDbi"))
```

Then install `ggmethylation`:

```r
remotes::install_github("ljwharbers/ggmethylation")
```

## Usage

The package ships a small synthetic modBAM (chr1, ~70 reads with HP tags,
5mC/5hmC calls, a deletion, an insertion locus, supplementary alignments and an
SNV at chr1:3000) and a matching VCF, used in the examples below.

```r
library(ggmethylation)

bam = system.file("extdata", "synthetic.bam", package = "ggmethylation")
vcf = system.file("extdata", "synthetic.vcf.gz", package = "ggmethylation")

# Parse a modBAM file (only loads reads in the specified region)
meth = read_methylation(
  bam,
  region = "chr1:1000-7000",
  mod_code = "m",       # 5mC (default)
  group_tag = "HP",     # group by haplotype tag
  max_reads = 200,      # downsample if more reads overlap
  seed = 1              # reproducible downsampling
)

plot_methylation(meth)
```

Secondary, duplicate and QC-fail alignments are skipped; supplementary
alignments are kept.

### Without grouping

When `group_tag = NULL` (the default) the read bars are not coloured by group
and the smooth panel shows a single overall line:

```r
meth_all = read_methylation(bam, "chr1:1000-7000")
plot_methylation(meth_all)
```

Reads can also be grouped by the base they carry at an SNV:

```r
meth_snv = read_methylation(bam, "chr1:1000-7000",
                            snv = list(position = 3000, ref = "A", alt = "G"))
plot_methylation(meth_snv)
```

### Customisation

`plot_methylation()` returns a patchwork object, so you can theme and modify it:

```r
library(ggplot2)

plot_methylation(
  meth,
  colour_low = "#EEEEEE",
  colour_high = "#1B5E20",
  line_width = 0.4,
  group_colours = c("1" = "steelblue", "2" = "coral"),
  smooth_span = 0.5,
  call_threshold = 0.5,       # binary methylated/unmethylated calls
  show_delta = TRUE,          # adds a panel: difference between the 2 groups
  panel_heights = c(4, 1, 1)  # read, smooth and delta panel
) & theme(text = element_text(size = 14))
```

`panel_heights` needs one value per panel, top to bottom (gene track, one read
panel per sample, smooth panel, delta panel); `NULL` uses the defaults.

### Sorting

Reads are packed into lanes per group, in the order given by `sort_by`. The
default is `"start"` for ungrouped data and
`c("start", "group", "mean_mod_prob")` for grouped data. Override it with any
column of `$reads`, or `mean_mod_prob`:

```r
plot_methylation(meth, sort_by = "mean_mod_prob")
```

### Several samples, variants, genes

```r
# One read panel per sample above a shared smooth panel
merged = merge_methylation(tumour = meth_tumour, normal = meth_normal)
plot_methylation(merged)

# VCF overlay: SNV carriers, DEL/DUP/INV spans, BND mate labels
vars = read_variants(vcf, "chr1:1000-7000")
plot_methylation(meth, variants = vars)

# Gene track from a TxDb or GTF/GFF file
genes = read_annotations(gtf = "genes.gtf", region = "chr1:1000-7000")
plot_methylation(meth, annotations = genes)
```

## Insertion-aware modifications

`read_methylation()` retains modification calls on inserted bases (bases
present in the read sequence but absent from the reference). These are
available via the `$insertion_sites` field and the `insertion_sites()` accessor.

Use `list_insertion_loci()` to cluster recurrent insertion events across reads,
then `plot_insertion_locus()` to visualise modification patterns in a stitched
coordinate system (left flank | insertion | right flank):

```r
# Loci carried by >= 2 reads. min_length = 50 ignores small (sequencing-error)
# insertions, sensible for real ONT/PacBio data; the default, 1, keeps all.
loci = list_insertion_loci(meth, min_reads = 2L, min_length = 50L)
loci

# Plot the first locus
plot_insertion_locus(meth, loci[1L, ])
```

## Supported modification types

Any base modification encoded in MM/ML tags (SAM spec), including combined
entries such as `C+mh`. Common codes:

| Code | Modification |
|------|-------------|
| `m`  | 5-methylcytosine (5mC) |
| `h`  | 5-hydroxymethylcytosine (5hmC) |
| `a`  | 6-methyladenine (6mA) |

Specify via the `mod_code` argument; several codes can be read at once
(`mod_code = c("m", "h")`).

## Requirements

- A BAM file with MM/ML base modification tags (modBAM format)
- A corresponding BAM index file (`.bai` or `.csi`)
- R >= 4.0

## License

MIT
