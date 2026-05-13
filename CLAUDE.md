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
`read_methylation()` is the main entry point. It queries a BAM file via Rsamtools, applies filters (MAPQ, strand, read length, downsampling), and calls `parse_mm_ml()` to extract per-base modification probabilities from MM/ML tags. It also handles optional grouping by BAM tag (e.g., HP for haplotype) or by SNV genotype. Returns an S3 object of class `methylation_data` with data frames for reference-aligned sites (`$sites`), insertion-based sites (`$insertion_sites`), and per-read summaries (`$reads`).

**Layer 2 — Processing**
- `parse_mm_ml.R`: Parses MM/ML SAM auxiliary tags. `seq_to_ref()` maps query positions to reference coordinates by walking the CIGAR string; positions inside CIGAR `I` operations return `NA` and are classified as insertion sites. `parse_mm_ml()` now returns a named list with `$sites` and `$insertion_sites`.
- `pack_reads.R`: Greedy interval scheduling algorithm that assigns reads to horizontal display lanes (like a genome browser). Used internally by `plot_methylation()` and `plot_insertion_locus()`.
- `smooth_methylation.R`: Aggregates per-site modification probabilities and fits a loess curve on a 200-point grid for the smoothed lower panel. The shared `.smooth_xy(x, y)` helper is used by both `plot_methylation()` and `plot_insertion_locus()`.
- `insertion_loci.R`: `list_insertion_loci()` clusters insertion events across reads into loci using a greedy single-pass algorithm. `insertion_sites()` is a convenience accessor for `$insertion_sites`.

**Layer 3 — Visualization**
- `plot_methylation()` (`R/plot_methylation.R`): takes a `methylation_data` object and returns a `patchwork` composite. Reads are drawn as horizontal bars; modification sites as coloured dots. When data is grouped, panels are produced per group and combined. Supports sorting by position, group, or mean modification probability; multi-modification codes via shape aesthetics; and custom colour gradients.
- `plot_insertion_locus()` (`R/plot_insertion_locus.R`): visualises a single insertion locus in a stitched coordinate system (left flank | insertion | right flank), showing carrier and non-carrier reads with an optional loess-smoothed comparison panel.

**Optional output**: `write_methylation.R` exports reads/sites to TSV or BED (with optional gzip).

## Key Data Structure

`methylation_data` (S3 list):
- `$reads`: data frame — `read_name`, `start`, `end`, `strand`, `bam_pos`, `mean_mod_prob`, optional `group`
- `$sites`: data frame — `read_name`, `position`, `mod_prob`, `mod_code`, optional `group`
- `$insertion_sites`: data frame — `read_name`, `ref_anchor`, `query_pos`, `ins_offset`, `ins_length`, `mod_prob`, `mod_code`, optional `group`. Zero rows when no insertion modifications were found.
- `$cigar_features`: data frame — one row per CIGAR operation per read; columns `type`, `ref_start`, `ref_end`, `query_start`, `query_end`, `length`, `read_name`. Used by `list_insertion_loci()`.
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
