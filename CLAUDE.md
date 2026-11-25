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

# Run a single test file (matches test-pack_reads.R)
devtools::test(filter = "pack_reads")

# Regenerate documentation (roxygen2)
devtools::document()

# Full R CMD check
devtools::check()
```

## Code style

Use `=` for assignment, not `<-`. The one exception is an assignment inside a call's arguments (e.g. `expect_silent(p <- plot_methylation(md))`), where `=` would pass a named argument instead.

## Architecture

The package has three layers:

**Layer 1 — Input** (`R/read_methylation.R`)
`read_methylation()` is the main entry point. It queries a BAM file via Rsamtools (`scanBamFlag()` drops secondary, QC-fail and duplicate alignments; supplementary alignments are kept), applies filters (MAPQ, strand, read length, downsampling with optional `seed`), and parses the MM/ML tags of each read once for all requested codes via `.parse_mods()`. It also handles optional grouping by BAM tag (`group_tag`, e.g. HP for haplotype; a tag found in no read gives a warning and an ungrouped object) or by SNV genotype (`snv = list(position, ref, alt)`). Each step is a small helper (`.bam_reads()`, `.filter_reads()`, `.snv_genotype()`, `.downsample_reads()`, `.parse_read_mods()`); `reads$.idx` maps rows back to the `scanBam()` record until it is dropped on return. Returns an S3 object of class `methylation_data` with data frames for reference-aligned sites (`$sites`), insertion-based sites (`$insertion_sites`), and per-read summaries (`$reads`).

**Layer 2 — Processing**
- `parse_mm_ml.R`: Parses MM/ML SAM auxiliary tags. `.parse_mods()` is the core (several codes per call); `parse_mm_ml()` is a thin single-code wrapper. `seq_to_ref()` maps query positions to reference coordinates by walking the CIGAR string; positions inside CIGAR `I` operations return `NA` and are classified as insertion sites. `parse_mm_ml()` returns a named list with `$sites` and `$insertion_sites`. Decisions:
  - Bases are counted on the original read; only the alignment strand (FLAG 0x10) decides the searched base and scan direction. The MM `+`/`-` field does not (a `G-m` entry is counted on G of the original read). `$sites` does not record that strand.
  - Combined codes (`C+mh`) interleave their ML values per base; numeric ChEBI codes are one code; every entry carrying the requested code contributes; canonical `N` matches any base.
  - An entry with **no** `?`/`.` flag behaves like `?` (unlisted bases not emitted). This deliberately deviates from the SAM spec (missing = `.`), because old unflagged basecallers list only CpGs.
  - Malformed or inconsistent tags (ML count, bad skip count, unparseable entry, skips past SEQ, empty SEQ, `MN` != SEQ length) reject the whole read: it keeps no calls, and `read_methylation()` gives one aggregated warning with the count.
- `pack_reads.R`: Greedy interval scheduling algorithm that assigns reads to horizontal display lanes (like a genome browser). Used internally by `plot_methylation()` and `plot_insertion_locus()`.
- `smooth_methylation.R`: Aggregates per-site modification probabilities and fits a loess curve (one per group, or per group × mod code) on a 200-point grid, or on a shared `grid` for the delta track. All fitting goes through `.smooth_group()` / `.loess_predict()`; `.smooth_xy(x, y)` is the single-group shortcut used by `plot_insertion_locus()`. Sites with an `NA` group are never smoothed. `.binarize_sites()` turns probabilities into 0/1 calls for binary mode.
- `insertion_loci.R`: `list_insertion_loci()` clusters insertion events (length >= `min_length`) across reads into loci with `.cluster_insertions()`: by anchor position (gap > `tol_pos` splits), then by length within each position group (running median, `tol_len`), then by position again. A carrier has an insertion within `tol_len` of the locus median length; a non-carrier spans the locus on both sides (start <= `anchor_pos - tol_pos`, end >= `anchor_pos + tol_pos`, using `bam_pos`/`bam_end`, not the region-clipped span) with no insertion there; reads with an other-length insertion are neither. Each row also carries its `carriers`/`noncarriers` (list-columns) and `tol_pos`/`tol_len`. `insertion_sites()` is a convenience accessor for `$insertion_sites`.

**Layer 3 — Visualization**
- `plot_methylation()` (`R/plot_methylation.R`): takes a `methylation_data` (or `multi_methylation_data`) object and returns a `patchwork` composite: optional gene track, one read panel per sample, a smooth panel, and an optional delta panel (two groups, single sample only). Single- and multi-sample data share one pipeline: `.prepare_reads()` (per-read means, sorting, per-group lane packing) → `build_read_panel()` (`R/build_read_panel.R`) → `.single_smooth_panel()` / `.multi_smooth_panel()`, both built by `.build_smooth_panel()` (colour/linetype columns + CI ribbon) → `.delta_panel()`. Binary calls are on exactly when `call_threshold` is non-`NULL`; multiple mod codes get one strip per code in each read bar and one smooth line (linetype) per code. Group colours for every panel come from one resolver, `.resolve_group_colours()` (`R/palettes.R`). Smooth lines break at a deletion carried by >= 75% of the group's reads spanning the whole deletion.
- `plot_insertion_locus(data, locus)` (`R/plot_insertion_locus.R`): visualises one row of `list_insertion_loci()` output in a stitched coordinate system (left flank | insertion | right flank), showing carrier and non-carrier reads (non-carriers below) with an optional loess-smoothed comparison panel. With anchor `A` and `W` = median length: reference bases `A - flank` .. `A - 1` keep their coordinate, inserted bases 1..`W` sit at `A` .. `A + W - 1`, reference bases from `A` on are shifted right by `W`; inserted bases beyond `W` are clipped.
- `read_variants()` returns a `variant_data` object that also stores `bnd_match_tol`, used by `build_variant_overlay()` to match reads' SA breakpoints to BND calls.

**Optional output**: `write_methylation.R` exports reads/sites to TSV or BED (with optional gzip).

## Key Data Structure

`methylation_data` (S3 list):
- `$reads`: data frame, one row per read — `read_name`, `start`, `end` (clipped to the region), `bam_pos`, `bam_end` (full alignment span: BAM POS and last reference base of the CIGAR), `strand`, `is_supplementary`, `sa_chrom`, `sa_pos`, `clip_side`, optional `group`. (`mean_mod_prob` and `lane` are added at plot time, not stored.)
- `$sites`: data frame, one row per read per call — `read_name`, `position`, `mod_prob`, `mod_code`, optional `group`
- `$insertion_sites`: data frame — `read_name`, `ref_anchor`, `query_pos`, `ins_offset`, `ins_length`, `mod_prob`, `mod_code`, optional `group`. Zero rows when no insertion modifications were found.
- `$cigar_features`: data frame — one row per `I`/`D`/`N` operation per read overlapping the region; columns `type`, `ref_start`, `ref_end`, `query_start`, `query_end`, `length`, `read_name`. Used by `list_insertion_loci()` and the deletion/insertion markers.
- `$sequences`, `$cigars`: per-read SEQ and CIGAR (named by `read_name`), for the variant overlay.
- `$region`: `GenomicRanges` object
- `$mod_code`: character vector of modification codes
- `$group_tag`: BAM tag name or `NULL`

### Per-read site presence

`$sites` does **not** contain every canonical base in the region for every read. Five mechanisms explain why a reference position may be absent for a given read:

1. **MM tag flag** — the basecaller's MM entry may use `?` (or no flag, deliberately treated like `?`), meaning unlisted canonical bases carry no information; only listed bases are emitted. When the entry uses `.`, unlisted canonical bases are emitted with `mod_prob = 0` (implicit-unmodified).
2. **Basecaller convention** — dorado 5mC typically lists all CpGs (high and low ML values); 6mA basecallers typically only list modified adenines.
3. **Sequence differences** — a SNV or sequencing error makes the base non-canonical at that position in a specific read; the MM delta walk never reaches it.
4. **CIGAR effects** — positions inside `D`/`N` have no query base; positions inside `I` have no reference coordinate (appear in `$insertion_sites` instead).
5. **Different alignment span** — a read simply does not cover that reference position.

## Test Data

Integration tests (`test-integration.R`) use a real BAM, by default:
```
/staging/leuven/stg_00096/home/lharbers/repositories/ggmethylation/data/PTCL8_PB_tumor_chr21_22_subset.bam
```
Point them elsewhere with the `GGMETHYLATION_TEST_BAM` and `GGMETHYLATION_TEST_REGION` environment variables. They are skipped when the file is unavailable. Unit tests have no external dependencies.

`tests/testthat/fixtures/` holds a small synthetic modBAM (`synthetic.bam`, ~70 reads on chr1:1000-7000 with HP tags, 5mC/5hmC calls, a shared deletion, an insertion locus at 6100, soft clips with SA tags and an A>G SNV at chr1:3000) and VCF (`synthetic.vcf.gz`: that SNV, a `<DEL>`, a BND to chr7:5000), plus indexes. An identical copy lives in `inst/extdata/` for the vignette and README (`system.file("extdata", "synthetic.bam", package = "ggmethylation")`); keep both in sync. Regenerate them with `data-raw/make_test_fixture.R` (instructions in its header). Shared test helpers (`make_test_data()`, `fixture_bam()`, `fixture_vcf()`) are in `tests/testthat/helper-make_test_data.R`.

## Documentation

All exported functions use roxygen2 with `markdown = TRUE`. After editing docs, run `devtools::document()` to regenerate `man/` and `NAMESPACE`. The vignette (`vignettes/ggmethylation.Rmd`) runs live on the `inst/extdata/` fixture (no cache); chunks needing optional packages are guarded with `eval = requireNamespace(...)`. Render it into `$TMPDIR`, not the repo. CI: `.github/workflows/R-CMD-check.yaml`.
