# ggmethylation — New Features Design Spec

**Date:** 2026-04-02
**Status:** Draft

---

## Context

ggmethylation currently supports single-sample, single-region read-level methylation visualization. Compared to tools like methylartist, three key capabilities are missing: gene annotation context, multi-sample comparison, and variant overlay. This spec adds all three as composable layer functions that integrate with the existing `plot_methylation()` pipeline.

## Architecture: Composable Layers

New data sources are handled by standalone reader functions that return S3 objects. These objects are passed to `plot_methylation()` as optional parameters. This keeps the existing API stable while making each feature independently testable and extensible.

```r
# Full workflow with all three features
md1   <- read_methylation("tumor.bam", "chr1:1000-5000", group_tag = "HP")
md2   <- read_methylation("normal.bam", "chr1:1000-5000", group_tag = "HP")
merged <- merge_methylation(tumor = md1, normal = md2)

annot <- read_annotations(txdb = TxDb.Hsapiens.UCSC.hg38.knownGene, region = "chr1:1000-5000")
vars  <- read_variants("variants.vcf.gz", "chr1:1000-5000")

plot_methylation(merged, annotations = annot, variants = vars)
```

All new dependencies (`GenomicFeatures`, `rtracklayer`, `VariantAnnotation`) go in `Suggests`, not `Imports`. Functions check availability at runtime with `requireNamespace()`.

---

## Feature 1: Gene Annotation Track

### New exported function: `read_annotations()`

```r
read_annotations(txdb = NULL, gtf = NULL, region)
```

- Exactly one of `txdb` or `gtf` must be provided
- `txdb`: a `TxDb` object (e.g. `TxDb.Hsapiens.UCSC.hg38.knownGene`)
- `gtf`: path to a GTF/GFF file — loaded via `rtracklayer::import()`, converted to TxDb via `GenomicFeatures::makeTxDbFromGRanges()`
- `region`: character region string (same format as `read_methylation()`)
- Extracts transcripts and exons overlapping the region
- Returns an S3 object of class `gene_annotations`

### `gene_annotations` S3 object

```
$transcripts  — data.frame: tx_name, gene_name, strand, tx_start, tx_end
$exons        — data.frame: tx_name, exon_start, exon_end
$region       — GRanges object
```

### New internal function: `build_gene_panel()`

```r
build_gene_panel(annotations, region_start, region_end)
```

Renders gene models in a ggplot2 panel:
- One row per transcript
- `geom_segment()`: thin horizontal lines for introns
- `geom_rect()`: thick boxes for exons
- Strand arrows: small arrowhead segments along intron lines indicating transcription direction
- `geom_text()`: gene name labels on the left margin
- Shared x-axis limits with other panels
- This panel carries the x-axis label ("Genomic position (bp)") as the new bottom panel

### Integration in `plot_methylation()`

- New parameter: `annotations = NULL`
- When non-NULL, `build_gene_panel()` is called and appended below the smooth panel in the patchwork
- `panel_heights` becomes dynamic (see Panel Heights section below)

### Dependencies (Suggests)

- `GenomicFeatures` — for TxDb queries
- `rtracklayer` — for GTF/GFF import

---

## Feature 2: Multi-Sample Support (Faceted Panels)

### New exported function: `merge_methylation()`

```r
merge_methylation(..., .list = NULL)
```

- Accepts named arguments: `merge_methylation(tumor = md1, normal = md2)`
- Or a named list: `merge_methylation(.list = list(tumor = md1, normal = md2))`
- Validates: all objects are `methylation_data`, all share the same region and mod_code
- Preserves existing within-sample grouping (e.g. haplotype groups) — the `group` column is untouched
- Returns an S3 object of class `multi_methylation_data`

### `multi_methylation_data` S3 object

```
$samples   — named list of methylation_data objects (one per sample)
$region    — GRanges (shared)
$mod_code  — character vector (shared)
```

### Integration in `plot_methylation()`

- `plot_methylation()` detects `multi_methylation_data` via `inherits()` and takes a faceted code path
- Each sample gets its own read panel, stacked vertically with a sample label
- Within each sample panel, within-sample grouping works as before (separate lane packing, separator lines, group colours)
- One shared smooth panel at the bottom with lines per sample (or per sample × group if within-sample grouping exists)
- Gene track (if annotations provided) at the very bottom, shared

**Layout with 2 samples + haplotype grouping + gene track:**

```
┌─────────────────────────┐
│ Tumor (HP1 + HP2)       │  ← per-sample read panel
├─────────────────────────┤
│ Normal (HP1 + HP2)      │  ← per-sample read panel
├─────────────────────────┤
│ Smooth (all)            │  ← shared, lines per sample/group
├─────────────────────────┤
│ Gene track              │  ← shared gene annotations
└─────────────────────────┘
```

### `print` and `summary` methods

- `print.multi_methylation_data()`: region, sample names, per-sample read counts
- `summary.multi_methylation_data()`: per-sample breakdown with grouping stats

---

## Feature 3: VCF Variant Overlay

### Changes to `read_methylation()`

Read sequences and CIGARs are now always stored in the `methylation_data` object:

- `$sequences`: named character vector keyed by `read_name` (the full read sequence)
- `$cigars`: named character vector keyed by `read_name`

These are already available from the `scanBam()` call — no additional I/O. The `seq` and `cigar` fields are retained after the existing processing instead of being discarded.

Also update `empty_methylation_data()` to include empty `$sequences` and `$cigars` fields.

### New exported function: `read_variants()`

```r
read_variants(vcf, region)
```

- `vcf`: path to a bgzipped, tabix-indexed VCF file (`.vcf.gz` + `.vcf.gz.tbi`)
- `region`: character region string
- Uses `VariantAnnotation::readVcf()` with `ScanVcfParam` restricted to the region
- Extracts: positions, ref alleles, alt alleles, variant type (SNV/indel)
- Returns an S3 object of class `variant_data`

### `variant_data` S3 object

```
$variants  — data.frame: position, ref, alt, type (SNV/insertion/deletion)
$region    — GRanges object
```

### Per-read base extraction (inside `plot_methylation()`)

When `variants` is provided and the data has stored sequences/CIGARs:

1. For each variant position and each read overlapping it:
   - Use `ref_to_seq()` (existing in `utils.R`) to map the genomic position to the query (read) position via CIGAR
   - Extract the base from `$sequences` at that query position
   - Classify: `ref` (matches reference allele), `alt` (matches alternate allele), `other`, or `del` (deletion at that position in CIGAR)
2. Produces a data.frame: `read_name`, `position`, `base`, `variant_class`

This logic lives in a new internal function `extract_variant_bases()` in a new file `R/variant_overlay.R`.

### Rendering in the read panel (top)

- At each variant position on each read: a small colored letter (`geom_text()`) overlaid on the read bar showing the actual base
- Color by variant class: ref = grey, alt = red/orange, other = yellow, deletion = black
- `geom_vline()` vertical dashed lines at all variant positions spanning the full panel

### Rendering in smooth and gene panels

- `geom_vline()` vertical dashed lines at variant positions (extending the existing `snv_position` pattern to multiple positions)

### Integration in `plot_methylation()`

- New parameter: `variants = NULL`
- Existing `snv_position` logic is preserved (backward compatible) — `snv_position` marker still works independently

### Dependencies (Suggests)

- `VariantAnnotation` — for VCF reading

---

## Cross-Cutting: Dynamic Panel Heights

`panel_heights` currently accepts a numeric vector of length 2 (`c(3, 1)` for reads + smooth). With the new features, the number of panels varies:

- Single sample, no gene track: `c(3, 1)` (unchanged)
- Single sample + gene track: `c(3, 1, 0.5)`
- Multi-sample (2), no gene track: `c(3, 3, 1)`
- Multi-sample (2) + gene track: `c(3, 3, 1, 0.5)`

**Behavior:**
- If user provides `panel_heights`, it must match the actual panel count — error if not
- If user provides `NULL` (new default), heights are auto-computed:
  - Each sample read panel: weight 3
  - Smooth panel: weight 1
  - Gene track: weight 0.5

---

## New Files

| File | Type | Purpose |
|------|------|---------|
| `R/read_annotations.R` | Exported | `read_annotations()` + `print.gene_annotations()` |
| `R/build_gene_panel.R` | Internal | `build_gene_panel()` — ggplot2 gene track renderer |
| `R/read_variants.R` | Exported | `read_variants()` + `print.variant_data()` |
| `R/variant_overlay.R` | Internal | `extract_variant_bases()` — per-read base extraction |
| `R/merge_methylation.R` | Exported | `merge_methylation()` + print/summary methods |

## Modified Files

| File | Changes |
|------|---------|
| `R/read_methylation.R` | Store `$sequences` and `$cigars` in the S3 object; update `empty_methylation_data()` |
| `R/plot_methylation.R` | Add `annotations`, `variants` parameters; add `multi_methylation_data` dispatch; dynamic panel heights; variant overlay layers |
| `DESCRIPTION` | Add `GenomicFeatures`, `rtracklayer`, `VariantAnnotation` to Suggests |
| `NAMESPACE` | Export `read_annotations`, `read_variants`, `merge_methylation` + S3 methods |

---

## Verification Plan

### Unit tests (no BAM needed)

1. `test-read_annotations.R`: mock TxDb with known transcripts/exons, verify `gene_annotations` object structure
2. `test-merge_methylation.R`: merge two synthetic `methylation_data` objects, verify `multi_methylation_data` structure, region/mod_code validation errors
3. `test-variant_overlay.R`: test `extract_variant_bases()` with known sequences, CIGARs, and variant positions — verify correct base extraction and classification
4. `test-build_gene_panel.R`: verify `build_gene_panel()` returns a ggplot object with expected layers

### Integration tests (require BAM + VCF on staging)

5. `read_annotations()` with a real TxDb for the BAM region → verify transcripts are found
6. `read_variants()` with a real VCF for the BAM region → verify variants are parsed
7. Full pipeline: `read_methylation()` → `merge_methylation()` → `plot_methylation(annotations = ..., variants = ...)` → verify patchwork composite has expected panel count
8. Existing tests still pass (backward compatibility)

### Manual verification

9. Generate a plot with all three features and visually inspect: gene models align with known genes, variant bases match expected genotypes, multi-sample panels are correctly separated
