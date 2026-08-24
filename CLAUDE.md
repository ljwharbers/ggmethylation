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

# Run a single test file (test_file() is defunct as of devtools 2.5.0)
devtools::test(filter = "pack_reads")

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

### Per-read site presence

`$sites` does **not** contain every canonical base in the region for every read. Five mechanisms explain why a reference position may be absent for a given read:

1. **MM tag flag** — the basecaller's MM entry may use `?` (or no flag), meaning unlisted canonical bases carry no information; only listed bases are emitted. When the entry uses `.`, unlisted canonical bases are emitted with `mod_prob = 0` (implicit-unmodified). `parse_mm_ml()` now honours the flag: `.` → emit zeros, `?` / none → omit.
2. **Basecaller convention** — dorado 5mC typically lists all CpGs (high and low ML values); 6mA basecallers typically only list modified adenines.
3. **Sequence differences** — a SNV or sequencing error makes the base non-canonical at that position in a specific read; the MM delta walk never reaches it.
4. **CIGAR effects** — positions inside `D`/`N` have no query base; positions inside `I` have no reference coordinate (appear in `$insertion_sites` instead).
5. **Different alignment span** — a read simply does not cover that reference position.

## Test Data

Integration tests run against a committed fixture — no external paths, no skipping:

```
tests/testthat/fixtures/hg002_fiberseq_MEG3.bam
```

HG002 PacBio Fiber-seq, longphase-phased, aligned to CHM13, subset to the MEG3
promoter/DMR (`chr14:95055000-95070000`). MEG3 is imprinted, so the two
haplotypes carry genuinely different methylation — the grouped tests assert that
real difference rather than merely that the code ran.

The fixture is deliberately chosen to exercise several parser paths at once:

| Property | Value | Exercises |
|---|---|---|
| reads | 52 | — |
| MM / ML tags | 52/52 | MM/ML parsing |
| HP tags | 51/52, split 30 / 21 | `group_tag = "HP"`, delta track, per-group smoothing |
| reads without HP | 1 | `drop_na_group` |
| mod codes | `C+m` (5mC) and `A+a` (6mA), both on every read | `mod_code` filtering, multi-code shapes |
| MM flag | none (bare `C+m,`) | the "omit unlisted canonical bases" branch |
| insertion sites | 54 | insertion-aware parsing, `list_insertion_loci()` |

Rebuild it with `data-raw/build_test_fixture.sh` (requires cluster access to the
source BAM). The exact counts in `test-integration.R` are regression anchors: the
fixture is frozen, so a change in any of them means parsing behaviour changed and
should be reviewed deliberately, not silently re-baselined.

Unit tests (pack_reads, parse_mm_ml, smooth, utils) have no external dependencies.

## Documentation

All exported functions use roxygen2 with `markdown = TRUE`. After editing docs, run `devtools::document()` to regenerate `man/` and `NAMESPACE`. The vignette is at `vignettes/ggmethylation.Rmd` and uses cached data from `inst/extdata/vignette_cache.rds`.
