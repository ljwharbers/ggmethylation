# ggmethylation

Read-level base modification visualisation for long-read sequencing data (Oxford Nanopore / PacBio).

Given a modBAM file and a genomic region, `ggmethylation` produces a composite plot:

- **Read panel:** individual reads as horizontal bars, with coloured marks at each modification site (colour encodes modification probability from the ML tag)
- **Smooth panel:** loess-smoothed mean modification probability per group, with an optional confidence ribbon

Optional panels sit above and below these: a **gene annotation track**, and a
**group-difference (delta) track**. Variant calls, CIGAR indels and
supplementary alignments can be overlaid on the read panel.

## Installation

`ggmethylation` depends on the following (Bioconductor) packages. Install them first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Rsamtools", "GenomicRanges", "GenomicAlignments", "IRanges", "txdbmaker"))
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
  line_width = 0.4,
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

## Grouping and group comparison

Reads can be split by a BAM tag or by genotype at an SNV, and the grouping
drives every panel.

```r
# By haplotype tag
meth_hp <- read_methylation("sample.bam", region, group_tag = "HP",
                            drop_na_group = TRUE)

# By genotype at a heterozygous SNV (works on unphased data)
meth_snv <- read_methylation("sample.bam", region,
                             snv_position = 95066012, ref_base = "C",
                             alt_base = "T")
```

### Confidence ribbon

`show_ci = TRUE` (the default) draws the loess confidence interval behind each
smoothed curve, so you can see where the summary is well supported.

### Group-difference (delta) track

With exactly two groups, `show_delta = TRUE` adds a panel showing the signed
difference between them, filled with the colour of whichever group is higher:

```r
plot_methylation(meth_hp, show_delta = TRUE)
```

### Binary call mode

`call_mode = "binary"` renders discrete methylated/unmethylated calls instead
of a continuous gradient. `call_ambiguous` adds a third state for
low-confidence calls, which are excluded from the aggregate:

```r
plot_methylation(meth_hp, call_mode = "binary", call_threshold = 0.5,
                 call_ambiguous = 0.2)
```

In binary mode the smooth panel switches to the *fraction of calls that are
methylated*.

## Gene annotation track

`read_annotations()` fetches a gene model for the region — from a cached UCSC
ncbiRefSeq GTF for `"hg38"` or `"chm13"`, or from your own `gtf`/`txdb`:

```r
ann <- read_annotations(genome = "chm13", region = region)
plot_methylation(meth_hp, annotations = ann)
```

## Variant overlays

`read_variants()` reads a VCF and classifies SNVs, indels and structural
variants (`DEL`, `DUP`, `INV`, `BND`):

```r
vars <- read_variants("calls.vcf.gz", region)
plot_methylation(meth_hp, variants = vars)
```

SNVs are marked only on the reads that actually carry the ALT allele.
Structural variants are drawn as spanning bars and BND records are labelled
with their mate location.

## CIGAR features and supplementary alignments

Large indels from the CIGAR string are shown by default (`show_cigar`,
`min_indel_size`), and reads with supplementary alignments get arrowheads
(`show_supplementary`). Consensus deletions also break the smoothed curve, so
the summary is not interpolated across regions with no read sequence.

## Multiple samples

`merge_methylation()` combines several objects over the same region into a
`multi_methylation_data`, plotted as stacked per-sample read panels above a
shared smooth panel:

```r
merged <- merge_methylation(tumour = md1, normal = md2)
plot_methylation(merged)
```

## Exporting

```r
write_methylation(meth, prefix = "out/sample")                    # TSV
write_methylation(meth, prefix = "out/sample",
                  format = "bed", gzip = TRUE)                    # BED + gzip
```

BED output is 0-based half-open with `mod_prob` in the score column.

## Insertion-aware modifications

`read_methylation()` retains modification calls on inserted bases (bases
present in the read sequence but absent from the reference). These are
available via the `$insertion_sites` field and the `insertion_sites()` accessor.

Use `list_insertion_loci()` to cluster recurrent insertion events across reads,
then `plot_insertion_locus()` to visualise modification patterns in a stitched
coordinate system (left flank | insertion | right flank):

```r
# Discover loci present in >= 2 reads
loci <- list_insertion_loci(meth, min_reads = 2L)
print(loci)

# Plot the first locus
plot_insertion_locus(meth, loci$locus_id[1L])
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
