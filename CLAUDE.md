# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`ggmethylation` is an R package for read-level base modification visualization from long-read sequencing data (ONT/PacBio). It parses modBAM files (BAM with MM/ML tags) and produces ggplot2-based composite figures with individual reads (top panel) and optional loess-smoothed modification probability per group (bottom panel).

## Common Commands

```r
# Install/rebuild package
devtools::install()

# Run all tests
devtools::test()

# Run a single test file
devtools::test_file("tests/testthat/test-pack_reads.R")

# Regenerate documentation (roxygen2)
devtools::document()

# Full R CMD check
devtools::check()
```

## Architecture

The package has three layers:

**Layer 1 — Input** (`R/read_methylation.R`)
`read_methylation()` is the main entry point. It queries a BAM file via Rsamtools, applies filters (MAPQ, strand, read length, downsampling), and calls `parse_mm_ml()` to extract per-base modification probabilities from MM/ML tags. It also handles optional grouping by BAM tag (e.g., HP for haplotype) or by SNV genotype. Returns an S3 object of class `methylation_data` with two data frames: `$reads` (one row per read) and `$sites` (one row per modified base).

**Layer 2 — Processing**
- `parse_mm_ml.R`: Parses MM/ML SAM auxiliary tags. `seq_to_ref()` maps query positions to reference coordinates by walking the CIGAR string.
- `pack_reads.R`: Greedy interval scheduling algorithm that assigns reads to horizontal display lanes (like a genome browser). Used internally by `plot_methylation()`.
- `smooth_methylation.R`: Aggregates per-site modification probabilities and fits a loess curve on a 200-point grid for the smoothed lower panel.

**Layer 3 — Visualization** (`R/plot_methylation.R`)
`plot_methylation()` takes a `methylation_data` object and returns a `patchwork` composite. Reads are drawn as horizontal bars; modification sites as coloured dots. When data is grouped, panels are produced per group and combined. Supports sorting by position, group, or mean modification probability; multi-modification codes via shape aesthetics; and custom colour gradients.

**Optional output**: `write_methylation.R` exports reads/sites to TSV or BED (with optional gzip).

## Key Data Structure

`methylation_data` (S3 list):
- `$reads`: data frame — `read_name`, `start`, `end`, `strand`, `lane`, `mean_mod_prob`, optional `group`
- `$sites`: data frame — `read_name`, `position`, `mod_prob`, `mod_code`, optional `group`
- `$region`: `GenomicRanges` object
- `$mod_code`: character vector of modification codes
- `$group_tag`: BAM tag name or `NULL`

## Test Data

Integration tests require a real BAM file at:
```
/staging/leuven/stg_00096/home/lharbers/repositories/ggmethylation/data/PTCL8_PB_tumor_chr21_22_subset.bam
```
Tests are automatically skipped if the file is unavailable. Unit tests (pack_reads, parse_mm_ml, smooth, utils) have no external dependencies.

## Documentation

All exported functions use roxygen2 with `markdown = TRUE`. After editing docs, run `devtools::document()` to regenerate `man/` and `NAMESPACE`. The vignette is at `vignettes/ggmethylation.Rmd` and uses cached data from `inst/extdata/vignette_cache.rds`.
