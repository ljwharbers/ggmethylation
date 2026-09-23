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
`read_methylation()` is the main entry point. It queries a BAM file via Rsamtools, applies filters (MAPQ, strand, read length, downsampling), and calls `parse_mm_ml()` to extract per-base modification probabilities from MM/ML tags. It also handles optional grouping by BAM tag (`group_tag`, e.g. HP for haplotype) or by SNV genotype (`snv = list(position, ref, alt)`). Each step is a small helper (`.bam_reads()`, `.filter_reads()`, `.snv_genotype()`, `.downsample_reads()`, `.parse_read_mods()`); `reads$.idx` maps rows back to the `scanBam()` record until it is dropped on return. Returns an S3 object of class `methylation_data` with data frames for reference-aligned sites (`$sites`), insertion-based sites (`$insertion_sites`), and per-read summaries (`$reads`).

**Layer 2 — Processing**
- `parse_mm_ml.R`: Parses MM/ML SAM auxiliary tags. `seq_to_ref()` maps query positions to reference coordinates by walking the CIGAR string; positions inside CIGAR `I` operations return `NA` and are classified as insertion sites. `parse_mm_ml()` returns a named list with `$sites` and `$insertion_sites`.
- `pack_reads.R`: Greedy interval scheduling algorithm that assigns reads to horizontal display lanes (like a genome browser). Used internally by `plot_methylation()` and `plot_insertion_locus()`.
- `smooth_methylation.R`: Aggregates per-site modification probabilities and fits a loess curve (one per group, or per group × mod code) on a 200-point grid, or on a shared `grid` for the delta track. All fitting goes through `.smooth_group()` / `.loess_predict()`; `.smooth_xy(x, y)` is the single-group shortcut used by `plot_insertion_locus()`. Sites with an `NA` group are never smoothed. `.binarize_sites()` turns probabilities into 0/1 calls for binary mode.
- `insertion_loci.R`: `list_insertion_loci()` clusters insertion events across reads into loci using a greedy single-pass algorithm; each row also carries its `carriers`/`noncarriers` (list-columns) and `tol_pos`/`tol_len`. `insertion_sites()` is a convenience accessor for `$insertion_sites`.

**Layer 3 — Visualization**
- `plot_methylation()` (`R/plot_methylation.R`): takes a `methylation_data` (or `multi_methylation_data`) object and returns a `patchwork` composite: optional gene track, one read panel per sample, a smooth panel, and an optional delta panel (two groups, single sample only). Single- and multi-sample data share one pipeline: `.prepare_reads()` (per-read means, sorting, per-group lane packing) → `build_read_panel()` (`R/build_read_panel.R`) → `.single_smooth_panel()` / `.multi_smooth_panel()`, both built by `.build_smooth_panel()` (colour/linetype columns + CI ribbon) → `.delta_panel()`. Binary calls are on exactly when `call_threshold` is non-`NULL`; multiple mod codes get one smooth line (linetype) per code.
- `plot_insertion_locus(data, locus)` (`R/plot_insertion_locus.R`): visualises one row of `list_insertion_loci()` output in a stitched coordinate system (left flank | insertion | right flank), showing carrier and non-carrier reads with an optional loess-smoothed comparison panel.
- `read_variants()` returns a `variant_data` object that also stores `bnd_match_tol`, used by `build_variant_overlay()` to match reads' SA breakpoints to BND calls.

**Optional output**: `write_methylation.R` exports reads/sites to TSV or BED (with optional gzip).

## Key Data Structure

`methylation_data` (S3 list):
- `$reads`: data frame — `read_name`, `start`, `end`, `bam_pos`, `strand`, `is_supplementary`, `sa_chrom`, `sa_pos`, `clip_side`, optional `group`. (`mean_mod_prob` and `lane` are added at plot time, not stored.)
- `$sites`: data frame — `read_name`, `position`, `mod_prob`, `mod_code`, optional `group`
- `$insertion_sites`: data frame — `read_name`, `ref_anchor`, `query_pos`, `ins_offset`, `ins_length`, `mod_prob`, `mod_code`, optional `group`. Zero rows when no insertion modifications were found.
- `$cigar_features`: data frame — one row per CIGAR operation per read; columns `type`, `ref_start`, `ref_end`, `query_start`, `query_end`, `length`, `read_name`. Used by `list_insertion_loci()`.
- `$region`: `GenomicRanges` object
- `$mod_code`: character vector of modification codes
- `$group_tag`: BAM tag name or `NULL`

### Per-read site presence

`$sites` does **not** contain every canonical base in the region for every read. Five mechanisms explain why a reference position may be absent for a given read:

1. **MM tag flag** — the basecaller's MM entry may use `?` (or no flag), meaning unlisted canonical bases carry no information; only listed bases are emitted. When the entry uses `.`, unlisted canonical bases are emitted with `mod_prob = 0` (implicit-unmodified).
2. **Basecaller convention** — dorado 5mC typically lists all CpGs (high and low ML values); 6mA basecallers typically only list modified adenines.
3. **Sequence differences** — a SNV or sequencing error makes the base non-canonical at that position in a specific read; the MM delta walk never reaches it.
4. **CIGAR effects** — positions inside `D`/`N` have no query base; positions inside `I` have no reference coordinate (appear in `$insertion_sites` instead).
5. **Different alignment span** — a read simply does not cover that reference position.

## Test Data

Integration tests require a real BAM file at:
```
/staging/leuven/stg_00096/home/lharbers/repositories/ggmethylation/data/PTCL8_PB_tumor_chr21_22_subset.bam
```
Tests are automatically skipped if the file is unavailable. Unit tests (pack_reads, parse_mm_ml, smooth, utils) have no external dependencies.

`tests/testthat/fixtures/` holds a small synthetic modBAM (`synthetic.bam`, ~70 reads on chr1 with HP tags, 5mC/5hmC calls, a shared deletion, an insertion locus, soft clips with SA tags and an SNV at chr1:3000) and VCF (`synthetic.vcf.gz`) used by `test-read_methylation.R`, `test-plot_methylation.R` and `test-read_variants_fixture.R`. Regenerate them with `data-raw/make_test_fixture.R` (instructions in its header).

## Documentation

All exported functions use roxygen2 with `markdown = TRUE`. After editing docs, run `devtools::document()` to regenerate `man/` and `NAMESPACE`. The vignette is at `vignettes/ggmethylation.Rmd` and uses cached data from `inst/extdata/vignette_cache.rds`.
