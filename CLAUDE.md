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

The package has five layers:

**Layer 1 — Input** (`R/read_methylation.R`)
`read_methylation()` is the main entry point. It queries a BAM file via Rsamtools, applies filters (MAPQ, strand, read length, downsampling), and calls `parse_mm_ml()` to extract per-base modification probabilities from MM/ML tags. It also handles optional grouping by BAM tag (e.g., HP for haplotype) or by SNV genotype. Returns an S3 object of class `methylation_data` with data frames for reference-aligned sites (`$sites`), insertion-based sites (`$insertion_sites`), and per-read summaries (`$reads`).

**Layer 2 — Processing**
- `parse_mm_ml.R`: Parses MM/ML SAM auxiliary tags. `seq_to_ref()` maps query positions to reference coordinates by walking the CIGAR string; positions inside CIGAR `I` operations return `NA` and are classified as insertion sites. `parse_mm_ml()` returns a named list with `$sites` and `$insertion_sites`.
- `pack_reads.R`: Greedy interval scheduling algorithm that assigns reads to horizontal display lanes (like a genome browser). Used internally by `plot_methylation()` and `plot_insertion_locus()`.
- `smooth_methylation.R`: Aggregates per-site modification probabilities and fits a loess curve on a 200-point grid for the smoothed lower panel, with `lower`/`upper` confidence columns. The shared `.smooth_xy(x, y)` helper is used by both `plot_methylation()` and `plot_insertion_locus()`.
- `insertion_loci.R`: `list_insertion_loci()` clusters insertion events across reads into loci using a greedy single-pass algorithm. `insertion_sites()` is a convenience accessor for `$insertion_sites`.
- `merge_methylation.R`: `merge_methylation()` combines several `methylation_data` objects covering the same region into a `multi_methylation_data` object, with `print`/`summary` methods.
- `utils.R`: region parsing, CIGAR decomposition, and shared validation helpers such as `.validate_sort_by()`.

**Layer 3 — Annotation and variant ingestion**
- `read_annotations.R`: `read_annotations()` builds a gene model for the region from a UCSC ncbiRefSeq GTF (downloaded and cached per `genome`, `"hg38"` or `"chm13"`), a user-supplied `gtf`, or a `TxDb`. Returns a `gene_annotations` object. `clear_annotation_cache()` drops the cached TxDb.
- `read_variants.R`: `read_variants()` reads a VCF via VariantAnnotation and classifies each record as SNV, insertion, deletion, or a structural variant (`DEL`/`DUP`/`INV`/`BND`). Returns a `variant_data` object. `parse_bnd_alt()` decodes BND ALT syntax.

**Layer 4 — Panel builders**
These do the actual drawing; `plot_methylation()` orchestrates them.
- `build_read_panel.R`: the read panel — read bars, per-site modification segments, CIGAR indels, supplementary-alignment arrowheads. `.classify_calls()` implements `call_mode = "binary"`.
- `build_gene_panel.R`: the gene annotation track.
- `delta_track.R`: `.compute_group_delta()` fits both groups on a shared position grid; `.build_delta_panel()` renders the signed difference as a diverging area.
- `variant_overlay.R` plus `variant_overlay_snv.R`, `variant_overlay_sv.R`, `variant_overlay_bnd.R`: variant layers over the read panel. SNV asterisks are drawn only on reads carrying the ALT allele, which needs the read sequences stored on the object.
- `palettes.R`: centralised colour constants, `theme_ggmethylation()`, and `.resolve_group_colours()`, which maps a group palette onto the group values actually present.

**Layer 5 — Composition** (`R/plot_methylation.R`)
- `plot_methylation()`: assembles the panels into a `patchwork` composite — gene track (top), reads, smooth, delta (bottom) — and applies sorting, packing, and scales. Optional panels are argument-gated: `annotations`, `variants`, `show_ci`, `show_delta`, `call_mode`, `show_cigar`, `show_supplementary`.
- `.plot_multi_methylation()`: the `multi_methylation_data` path — stacked per-sample read panels above one shared smooth panel. `show_delta` is not supported here and warns.
- `plot_insertion_locus()` (`R/plot_insertion_locus.R`): visualises a single insertion locus in a stitched coordinate system (left flank | insertion | right flank), showing carrier and non-carrier reads with an optional loess-smoothed comparison panel.

**Optional output**: `write_methylation.R` exports reads/sites to TSV or BED (with optional gzip).

## Key Data Structure

`methylation_data` (S3 list):
- `$reads`: data frame — `read_name`, `start`, `end`, `strand`, `bam_pos`, `mean_mod_prob`, optional `group`
- `$sites`: data frame — `read_name`, `position`, `mod_prob`, `mod_code`, optional `group`
- `$insertion_sites`: data frame — `read_name`, `ref_anchor`, `query_pos`, `ins_offset`, `ins_length`, `mod_prob`, `mod_code`, optional `group`. Zero rows when no insertion modifications were found.
- `$cigar_features`: data frame — one row per CIGAR operation per read; columns `type`, `ref_start`, `ref_end`, `query_start`, `query_end`, `length`, `read_name`. Used by `list_insertion_loci()`.
- `$sequences` / `$cigars`: per-read query sequences and CIGAR strings, keyed by read name. Required by the SNV overlay to determine which reads carry the ALT allele; objects created before these existed lose the SNV layer with a warning.
- `$region`: `GenomicRanges` object
- `$mod_code`: character vector of modification codes
- `$group_tag`: BAM tag name (e.g. `"HP"`), `"SNV"` for genotype grouping, or `NULL`
- `$snv_position`: the SNV used for grouping, or `NULL`

`multi_methylation_data` (from `merge_methylation()`) holds `$samples` (a named list of `methylation_data`), `$region`, and `$mod_code`.

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
