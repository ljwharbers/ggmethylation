# Plan C10: Unit Tests

## Goal

Add a comprehensive testthat test suite covering all internal functions and the public API, using only synthetic data for unit tests and the real modBAM for an optional integration test.

## Affected files

**Created:**
- `tests/testthat/test-parse_mm_ml.R`
- `tests/testthat/test-utils.R`
- `tests/testthat/test-pack_reads.R`
- `tests/testthat/test-smooth_methylation.R`
- `tests/testthat/test-integration.R`

**Modified:**
- `DESCRIPTION` — add `withr` to Suggests (needed for `local_` helpers if used; optional)
- `tests/testthat.R` — create the standard testthat entry-point if not already present

No production R files need changes. `testthat (>= 3.0.0)` is already in Suggests and `Config/testthat/edition: 3` is already set in DESCRIPTION.

---

## Implementation steps

### Step 1 — Create `tests/testthat.R`

Create the standard file that devtools/R CMD check uses to discover and run the suite:

```r
library(testthat)
library(ggmethylation)
test_check("ggmethylation")
```

### Step 2 — `test-utils.R`: parse_region, ref_to_seq, complement_base

Write tests for the three pure utility functions that need no BAM.

### Step 3 — `test-parse_mm_ml.R`: seq_to_ref + parse_mm_ml

These two functions are tightly coupled. Test `seq_to_ref` independently first, then test the full `parse_mm_ml` pipeline with synthetic MM/ML strings.

### Step 4 — `test-pack_reads.R`: pack_reads

Test the greedy interval scheduler with known inputs whose lane assignments can be computed by hand.

### Step 5 — `test-smooth_methylation.R`: smooth_methylation

Test the aggregation + loess path and the fallback path (< 4 unique positions).

### Step 6 — `test-integration.R`: read_methylation + plot_methylation

Skip automatically when the BAM is absent. When it is present, validate both ungrouped and HP-grouped calls end-to-end.

---

## Specific test cases

### `test-utils.R`

#### `parse_region()`

```
test_that("parse_region parses standard region", { ... })
```
- Input: `"chr1:1000-2000"` → `list(chrom="chr1", start=1000L, end=2000L)`

```
test_that("parse_region strips commas", { ... })
```
- Input: `"chr1:1,000-2,000"` → same result as above

```
test_that("parse_region rejects missing dash", { ... })
```
- Input: `"chr1:1000"` → `expect_error(..., "Invalid region format")`

```
test_that("parse_region rejects missing colon", { ... })
```
- Input: `"chr11000-2000"` → `expect_error(..., "Invalid region format")`

```
test_that("parse_region rejects empty string", { ... })
```
- Input: `""` → `expect_error`

```
test_that("parse_region accepts dotted chromosome names", { ... })
```
- Input: `"chr1.1:100-200"` → chrom `"chr1.1"`, start `100L`, end `200L`

#### `complement_base()`

```
test_that("complement_base returns correct complements", { ... })
```
- A→T, T→A, C→G, G→C (both upper and lower via `toupper` in production code)

```
test_that("complement_base errors on unknown base", { ... })
```
- Input: `"N"` → `expect_error(..., "Unknown base")`

#### `ref_to_seq()`

```
test_that("ref_to_seq maps M-only alignment", { ... })
```
- CIGAR: `"10M"`, ref_start: `100`, target: `103` → query pos `4`
- target: `109` → query pos `10`
- target: `110` (beyond read) → `NA_integer_`

```
test_that("ref_to_seq accounts for leading soft clip", { ... })
```
- CIGAR: `"5S10M"`, ref_start: `100`, target: `100` → query pos `6` (soft clip shifts q_pos by 5)
- target: `104` → query pos `10`

```
test_that("ref_to_seq returns NA in deletion", { ... })
```
- CIGAR: `"5M2D5M"`, ref_start: `1`, target: `7` (inside the 2D) → `NA_integer_`

```
test_that("ref_to_seq handles insertion (query only)", { ... })
```
- CIGAR: `"5M2I5M"`, ref_start: `1`, target: `8` → query pos `9` (2I consumed 2 query positions)

```
test_that("ref_to_seq handles hard clip (H)", { ... })
```
- CIGAR: `"3H5M"`, ref_start: `1`, target: `3` → query pos `3`
  (H does not consume query in BAM SEQ, so q_pos starts at 1)

---

### `test-parse_mm_ml.R`

#### `seq_to_ref()`

The function is defined in `parse_mm_ml.R` and is internal. Access it via `ggmethylation:::seq_to_ref`.

```
test_that("seq_to_ref M-only alignment", { ... })
```
- CIGAR: `"10M"`, pos: `1000`, query_positions: `c(1L, 5L, 10L)`
  → ref positions: `c(1000L, 1004L, 1009L)`

```
test_that("seq_to_ref soft clip shifts query indices", { ... })
```
- CIGAR: `"4S6M"`, pos: `100`, query_positions: `c(1L, 4L, 5L, 10L)`
  → `c(NA, NA, 100L, 105L)` (first 4 are soft-clipped → NA)

```
test_that("seq_to_ref deletion skips reference positions", { ... })
```
- CIGAR: `"5M3D5M"`, pos: `100`, query_positions: `c(5L, 6L)`
  → `c(104L, 108L)` (deletion consumes ref 105–107, next match starts at 108)

```
test_that("seq_to_ref insertion positions map to NA", { ... })
```
- CIGAR: `"3M2I4M"`, pos: `1`, query_positions: `c(3L, 4L, 5L, 6L)`
  → `c(3L, NA, NA, 4L)` (positions 4-5 are in the insertion)

#### `parse_mm_ml()`

Access via `ggmethylation:::parse_mm_ml`.

**Synthetic data construction rationale:** The MM tag format is `BASE+code,delta0,delta1,...;`. Delta values are skip counts between successive modified bases among all occurrences of the canonical base in the read. `delta=0` means the next canonical base is modified; `delta=1` means skip 1 canonical base then modify the next.

**Test case 1 — Forward strand, simple M CIGAR, CpG context**

```r
# Sequence: ACGACGACG (9 bases)
# C positions (1-based): 2, 5, 8
# MM tag: "C+m,0,1;" — delta=0 means modify 1st C (pos 2),
#                       delta=1 means skip one C then modify (pos 8)
# ML values: c(230L, 128L)  → mod_prob = 230/255, 128/255
# CIGAR: "9M", pos: 1000
# Expected: position=c(1001, 1007), mod_prob=c(230/255, 128/255)

seq      <- "ACGACGACG"
mm_tag   <- "C+m,0,1;"
ml_tag   <- as.integer(c(230, 128))
result   <- ggmethylation:::parse_mm_ml(seq, mm_tag, ml_tag, "m", "+", "9M", 1000L)
expect_equal(result$position, c(1001L, 1007L))
expect_equal(result$mod_prob, c(230/255, 128/255))
```

```
test_that("parse_mm_ml forward strand simple CpG", { ... })  # as above
```

**Test case 2 — Mod code not in MM tag returns empty df**

```
test_that("parse_mm_ml returns empty df when mod_code absent", { ... })
```
- `mm_tag = "C+m,0;"`, `ml_tag = c(200L)`, `mod_code = "h"`
  → nrow == 0, columns `position` (integer) and `mod_prob` (numeric)

**Test case 3 — Empty / NA mm_tag**

```
test_that("parse_mm_ml returns empty df on NULL mm_tag", { ... })
```
- `mm_tag = NULL` → empty df

```
test_that("parse_mm_ml returns empty df on NA mm_tag", { ... })
```
- `mm_tag = NA_character_` → empty df

```
test_that("parse_mm_ml returns empty df on empty ml_tag", { ... })
```
- `ml_tag = integer(0)` → empty df

**Test case 4 — Multiple modifications in MM tag, correct ML offset**

```r
# MM: "C+m,0;C+h,0;" — two mods, each with 1 modified base
# ML: c(200L, 100L)  — first for m, second for h
# Sequence "ACGACG", pos=1, CIGAR="6M"
# C positions: 2, 5
# mod_code "m": delta=0 → 1st C → pos 2 → ref 2
# mod_code "h": delta=0 → 1st C → pos 2 → ref 2
```
```
test_that("parse_mm_ml multi-mod tag correct ML offset for h", { ... })
```
- With `mod_code="h"`, expect `position=2L, mod_prob=100/255`

**Test case 5 — Reverse strand alignment**

On a reverse-strand read the MM tag reports deltas relative to the 5'→3' direction of the original molecule (i.e., right-to-left in the reference). When `strand="-"` and `mm_strand="+"`, the function complements the canonical base (C→G) and scans the sequence in reverse.

```r
# Reverse-strand read, original molecule 5'→3': ...G C G...
# BAM SEQ is the reverse-complement, so if original is "GCGACG"
# the BAM SEQ is "CGTCGC"
# G positions in reverse of "CGTCGC": we look for G and reverse scan
#
# Simpler construction:
# seq = "CGCGCG" (6 bases), strand = "-"
# MM: "C+m,0,0;" → mod strand "+", canonical C
# For rev-strand: search_base = complement("C") = "G"
# G positions in "CGCGCG": 2, 4, 6 → reversed: 6, 4, 2
# delta=0 → index 1 → pos 6; delta=0 → index 2 → pos 4
# CIGAR "6M", pos=1000
# ref positions: seq_to_ref("6M", 1000, c(6,4)) = c(1005, 1003)
```
```
test_that("parse_mm_ml reverse strand scans complement in reverse", { ... })
```
- Expected positions: `c(1005L, 1003L)` (order from the delta walk)

**Test case 6 — Soft clip: positions in clip map to NA (filtered from output)**

```r
# seq = "SSSSCCCCC" (9 bases), first 4 are soft-clipped
# CIGAR "4S5M", pos=100
# C at positions 5,6,7,8,9 (all in the 5M span)
# C at positions 1,2,3 — soft-clipped (these would be in "SSSS" which
#   happen to be S chars, but let's use real C's in clip to test filtering)
# Better: seq="CCCCAAAAA", only C's in soft clip
# MM: "C+m,0,0;" → mods at 1st and 2nd C (query pos 1 and 2 → ref NA)
# Expected: empty df (all positions filtered as NA)
```
```
test_that("parse_mm_ml filters sites in soft clip (NA ref pos)", { ... })
```

**Test case 7 — Delta walk with mixed CIGAR (deletion)**

```r
# seq = "ACGACG" (6 bases, no clip)
# CIGAR "3M2D3M", pos=10
# C at positions: 2, 5
# MM: "C+m,0,1;" → modify pos 2 (ref 11) and pos 8 (skip pos 5 → next C
#   but pos 5 is in the match after deletion)
# Actually: delta=0 → 1st C (query pos 2) → ref = 10 + 2 - 1 = 11
#           delta=1 → skip 1 C → 3rd C would be needed, but there are only 2
#           → valid check filters it out → only one row returned
# Let's verify: delta=0 → index 1; delta=1 → index 3; but length(C positions)=2
#   so index 3 is out of bounds → filtered → 1 row
# OR use "C+m,0;" (just the first C) to test deletion position mapping:
#   query pos 2 → ref 11 (first 3M starts at 10, query 2 → ref 10+1=11) ✓
```
```
test_that("parse_mm_ml with deletion CIGAR maps positions correctly", { ... })
```

**Test case 8 — All deltas out-of-bounds returns empty df**

```r
# seq = "ACGT", 1 C at pos 2
# MM: "C+m,1;" → delta=1 → index 2 → but only 1 C → out of bounds
```
```
test_that("parse_mm_ml returns empty when all deltas out of bounds", { ... })
```

**Test case 9 — mod_prob values are scaled correctly (0–1)**

```
test_that("parse_mm_ml mod_prob is ml_value / 255", { ... })
```
- ML value 255 → mod_prob 1.0; ML value 0 → mod_prob 0.0; ML value 128 → 128/255

**Test case 10 — Trailing semicolon in MM tag**

```
test_that("parse_mm_ml handles trailing semicolon in MM tag", { ... })
```
- `mm_tag = "C+m,0;"` and `mm_tag = "C+m,0"` should give identical results

---

### `test-pack_reads.R`

Access via `ggmethylation:::pack_reads`.

```
test_that("pack_reads returns integer(0) for empty input", { ... })
```
- Input: zero-row data.frame with columns `read_name`, `start`, `end`
  → `integer(0)`

```
test_that("pack_reads assigns lane 1 to single read", { ... })
```
- One read: start=100, end=200 → `c(1L)`

```
test_that("pack_reads packs non-overlapping reads into lane 1", { ... })
```
- Three reads: [100,200], [300,400], [500,600], gap=10
  → all in lane 1 (300 > 200+10=210, 500 > 400+10=410)
  → `c(1L, 1L, 1L)`

```
test_that("pack_reads assigns overlapping reads to different lanes", { ... })
```
- Two reads: [100,300], [200,400], gap=10
  → read 2 start (200) is NOT > 300+10=310, so new lane
  → `c(1L, 2L)`

```
test_that("pack_reads respects gap parameter", { ... })
```
- Reads [100,200] and [205,300], gap=10:
  → 205 > 200+10=210? No → different lanes → `c(1L, 2L)`
- Same reads, gap=0:
  → 205 > 200+0=200? Yes → same lane → `c(1L, 1L)`

```
test_that("pack_reads recycles lanes efficiently", { ... })
```
- Four reads sorted: [1,100], [1,100], [101,200], [101,200]
  → reads 1,2 go to lanes 1,2; read 3 start=101>100+10=110? No→new lane...
  → Actually with gap=10: 101>100+10=110 is FALSE → lanes 3,4
  → But with gap=0: 101>100 → reads 3,4 fit back in lanes 1,2
  → Test with gap=0: `c(1L, 2L, 1L, 2L)`

```
test_that("pack_reads result length matches nrow(reads)", { ... })
```
- n reads in → n integers out

---

### `test-smooth_methylation.R`

Access via `ggmethylation:::smooth_methylation`.

```
test_that("smooth_methylation returns empty df on NULL input", { ... })
```
- `smooth_methylation(NULL)` → 0-row df with columns `position`, `mean_prob`, `group`

```
test_that("smooth_methylation returns empty df on 0-row input", { ... })
```
- Input: 0-row df → same as above

```
test_that("smooth_methylation returns raw means when < 4 unique positions", { ... })
```
- Three unique positions: `position=c(1,1,2,3)`, `mod_prob=c(0.8,0.6,0.5,0.2)`, `group="A"`
  → nrow=3, mean at pos 1 = 0.7, at pos 2 = 0.5, at pos 3 = 0.2
  → result$position == c(1,2,3), result$mean_prob == c(0.7, 0.5, 0.2)

```
test_that("smooth_methylation fits loess when >= 4 unique positions", { ... })
```
- Five positions with clear trend: `position=1:5*100`, `mod_prob=c(0.1,0.2,0.3,0.4,0.5)`, `group="A"`
  → returns 200 rows (the grid), all `mean_prob` in [0,1] range
  → result$position range matches input range

```
test_that("smooth_methylation handles multiple groups independently", { ... })
```
- Groups "A" (5 positions) and "B" (5 positions)
  → result has rows for both groups
  → `unique(result$group)` == `c("A","B")` (or same set)

```
test_that("smooth_methylation respects custom group_col name", { ... })
```
- Input df has column `"haplotype"` instead of `"group"`
  → `smooth_methylation(sites, group_col="haplotype")` → result column is `"haplotype"`

```
test_that("smooth_methylation result has expected columns", { ... })
```
- Columns are exactly `c("position","mean_prob","group")` (or the custom group_col)

```
test_that("smooth_methylation drops NA group values", { ... })
```
- Rows with `group=NA` are excluded; only non-NA groups appear in output

---

## Integration test (`test-integration.R`)

```r
bam_path <- "/staging/leuven/stg_00096/home/lharbers/repositories/ggmethylation/data/PTCL8_PB_tumor_chr21_22_subset.bam"

skip_if_not(file.exists(bam_path), "Real BAM not available")
```

```
test_that("read_methylation returns methylation_data for real BAM (ungrouped)", { ... })
```
- `read_methylation(bam_path, "chr21:10000000-10010000")`
- Result inherits `"methylation_data"`, has `reads` and `sites` data.frames
- `reads` has columns `read_name`, `start`, `end`, `strand`
- `sites` has columns `position`, `mod_prob`, `read_name`, `mod_code`
- All `mod_prob` values in [0,1]

```
test_that("read_methylation returns methylation_data with HP grouping", { ... })
```
- `read_methylation(bam_path, "chr21:10000000-10010000", group_tag="HP")`
- `reads$group` is non-null (may be NA if HP tag absent, but column exists)

```
test_that("plot_methylation returns a ggplot/patchwork for real data", { ... })
```
- Call `plot_methylation(md)` on an ungrouped result
- Result inherits from `"gg"` or `"patchwork"`

```
test_that("read_methylation warns on empty region", { ... })
```
- Use a region with no reads (e.g., a gap chromosome or out-of-range coordinates)
- `expect_warning(read_methylation(...), "No reads found")`

---

## Edge cases to handle

- `parse_mm_ml`: MM tag with `?` or `.` flag characters (already stripped by `gsub("[?.]", "", entry)`) — add a test confirming they don't break parsing.
- `seq_to_ref`: query_positions vector longer than the read (out-of-bounds indexing into `query_to_ref`) — confirm `NA` is returned rather than an error.
- `pack_reads`: reads that are in non-start-sorted order — the function processes them in the order provided, so test that reads given in reverse order still produce valid (if suboptimal) lane assignments.
- `smooth_methylation`: single unique position (< 4 threshold) — 1-row output with the mean, no loess.
- `parse_region`: chromosome names with underscores (e.g., `"chr1_random:1-100"`) — should parse correctly via `[\w.]+`.

---

## Dependencies

No new packages required. The suite uses only:
- `testthat (>= 3.0.0)` — already in `Suggests`
- All production-code dependencies are already imported

If `withr` helpers (e.g., `withr::local_tempfile`) are needed for the integration test, add `withr` to `Suggests`. It is not strictly necessary since the BAM path is absolute.
