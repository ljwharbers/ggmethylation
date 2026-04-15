# Tests for parse_region, ref_to_seq, and complement_base in R/utils.R

# --- parse_region ---

test_that("parse_region parses standard region", {
  result <- ggmethylation:::parse_region("chr1:1000-2000")
  expect_equal(result$chrom, "chr1")
  expect_equal(result$start, 1000L)
  expect_equal(result$end, 2000L)
})

test_that("parse_region strips commas", {
  result <- ggmethylation:::parse_region("chr1:1,000-2,000")
  expect_equal(result$chrom, "chr1")
  expect_equal(result$start, 1000L)
  expect_equal(result$end, 2000L)
})

test_that("parse_region rejects missing dash", {
  expect_error(ggmethylation:::parse_region("chr1:1000"), "Invalid region format")
})

test_that("parse_region rejects missing colon", {
  expect_error(ggmethylation:::parse_region("chr11000-2000"), "Invalid region format")
})

test_that("parse_region rejects empty string", {
  expect_error(ggmethylation:::parse_region(""))
})

test_that("parse_region accepts dotted chromosome names", {
  result <- ggmethylation:::parse_region("chr1.1:100-200")
  expect_equal(result$chrom, "chr1.1")
  expect_equal(result$start, 100L)
  expect_equal(result$end, 200L)
})

# --- ref_to_seq ---
# ref_to_seq(cigar, ref_start, target_ref_pos) -> single integer query position or NA

test_that("ref_to_seq M-only alignment maps correctly", {
  # 10M from ref 1000; ref pos 1000 -> query 1, ref pos 1004 -> query 5
  expect_equal(ggmethylation:::ref_to_seq("10M", 1000L, 1000L), 1L)
  expect_equal(ggmethylation:::ref_to_seq("10M", 1000L, 1004L), 5L)
  expect_equal(ggmethylation:::ref_to_seq("10M", 1000L, 1009L), 10L)
})

test_that("ref_to_seq returns NA for position outside alignment", {
  expect_equal(ggmethylation:::ref_to_seq("10M", 1000L, 999L), NA_integer_)
  expect_equal(ggmethylation:::ref_to_seq("10M", 1000L, 1010L), NA_integer_)
})

test_that("ref_to_seq returns NA for position in deletion", {
  # 5M3D5M from ref 100; deletion covers ref 105-107
  expect_equal(ggmethylation:::ref_to_seq("5M3D5M", 100L, 105L), NA_integer_)
  expect_equal(ggmethylation:::ref_to_seq("5M3D5M", 100L, 107L), NA_integer_)
})

test_that("ref_to_seq maps correctly after deletion", {
  # 5M3D5M from ref 100; ref 108 is first base of second 5M -> query pos 6
  expect_equal(ggmethylation:::ref_to_seq("5M3D5M", 100L, 108L), 6L)
})

test_that("ref_to_seq soft clip shifts effective query positions", {
  # 4S6M from ref 100; ref 100 -> query 5, ref 105 -> query 10
  expect_equal(ggmethylation:::ref_to_seq("4S6M", 100L, 100L), 5L)
  expect_equal(ggmethylation:::ref_to_seq("4S6M", 100L, 105L), 10L)
  # Soft-clipped region has no reference position -> query 1-4 are not reachable by ref coord
  expect_equal(ggmethylation:::ref_to_seq("4S6M", 100L, 99L), NA_integer_)
})

# --- parse_sa_tag ---

test_that("parse_sa_tag returns empty data.frame for NULL input", {
  result <- ggmethylation:::parse_sa_tag(NULL)
  expect_equal(nrow(result), 0L)
  expect_equal(names(result), c("rname", "pos", "strand", "cigar", "mapq", "nm"))
})

test_that("parse_sa_tag returns empty data.frame for NA input", {
  result <- ggmethylation:::parse_sa_tag(NA_character_)
  expect_equal(nrow(result), 0L)
})

test_that("parse_sa_tag returns empty data.frame for empty string", {
  result <- ggmethylation:::parse_sa_tag("")
  expect_equal(nrow(result), 0L)
})

test_that("parse_sa_tag parses a single SA entry", {
  result <- ggmethylation:::parse_sa_tag("chr5,45000,+,50M,60,0;")
  expect_equal(nrow(result), 1L)
  expect_equal(result$rname,  "chr5")
  expect_equal(result$pos,    45000L)
  expect_equal(result$strand, "+")
  expect_equal(result$cigar,  "50M")
  expect_equal(result$mapq,   60L)
  expect_equal(result$nm,     0L)
})

test_that("parse_sa_tag parses multiple SA entries", {
  result <- ggmethylation:::parse_sa_tag("chr5,45000,+,50M,60,0;chr8,12000,-,30M,0,1;")
  expect_equal(nrow(result), 2L)
  expect_equal(result$rname,  c("chr5", "chr8"))
  expect_equal(result$pos,    c(45000L, 12000L))
  expect_equal(result$strand, c("+", "-"))
  expect_equal(result$mapq,   c(60L, 0L))
})

test_that("parse_sa_tag silently skips malformed entries", {
  # Only second entry is valid
  result <- ggmethylation:::parse_sa_tag("bad;chr8,12000,-,30M,0,1;")
  expect_equal(nrow(result), 1L)
  expect_equal(result$rname, "chr8")
})

# --- complement_base ---

test_that("complement_base returns correct complements", {
  expect_equal(ggmethylation:::complement_base("A"), "T")
  expect_equal(ggmethylation:::complement_base("T"), "A")
  expect_equal(ggmethylation:::complement_base("C"), "G")
  expect_equal(ggmethylation:::complement_base("G"), "C")
})

test_that("complement_base errors on unknown base", {
  expect_error(ggmethylation:::complement_base("N"), "Unknown base")
})

# --- cigar_ref_width ---

test_that("cigar_ref_width computes pure M width", {
  expect_equal(ggmethylation:::cigar_ref_width("10M"), 10L)
})

test_that("cigar_ref_width ignores soft clips", {
  expect_equal(ggmethylation:::cigar_ref_width("5S10M3S"), 10L)
})

test_that("cigar_ref_width ignores hard clips", {
  expect_equal(ggmethylation:::cigar_ref_width("5H10M5H"), 10L)
})

test_that("cigar_ref_width includes deletions", {
  expect_equal(ggmethylation:::cigar_ref_width("5M3D5M"), 13L)
})

test_that("cigar_ref_width excludes insertions", {
  expect_equal(ggmethylation:::cigar_ref_width("5M2I5M"), 10L)
})

test_that("cigar_ref_width returns 0 for star", {
  expect_equal(ggmethylation:::cigar_ref_width("*"), 0L)
})

test_that("cigar_ref_width returns 0 for NA", {
  expect_equal(ggmethylation:::cigar_ref_width(NA_character_), 0L)
})

test_that("cigar_ref_width handles complex CIGAR", {
  expect_equal(ggmethylation:::cigar_ref_width("100S500M200D800M100S"), 1500L)
})

# --- detect_clip_side ---

test_that("detect_clip_side returns NA for no clips", {
  expect_equal(ggmethylation:::detect_clip_side("10M"), NA_character_)
})

test_that("detect_clip_side detects left soft clip", {
  expect_equal(ggmethylation:::detect_clip_side("5S10M"), "left")
})

test_that("detect_clip_side detects right soft clip", {
  expect_equal(ggmethylation:::detect_clip_side("10M3S"), "right")
})

test_that("detect_clip_side detects both soft clips", {
  expect_equal(ggmethylation:::detect_clip_side("5S10M3S"), "both")
})

test_that("detect_clip_side detects hard clips", {
  expect_equal(ggmethylation:::detect_clip_side("5H10M"), "left")
})

test_that("detect_clip_side returns NA for star", {
  expect_equal(ggmethylation:::detect_clip_side("*"), NA_character_)
})

test_that("detect_clip_side returns NA for NA input", {
  expect_equal(ggmethylation:::detect_clip_side(NA_character_), NA_character_)
})

# --- region_to_granges ---

test_that("region_to_granges returns GRanges with correct coordinates", {
  gr <- ggmethylation:::region_to_granges("chr1:1000-2000")
  expect_s4_class(gr, "GRanges")
  expect_equal(as.character(GenomicRanges::seqnames(gr)), "chr1")
  expect_equal(GenomicRanges::start(gr), 1000L)
  expect_equal(GenomicRanges::end(gr), 2000L)
})

test_that("region_to_granges handles comma-formatted positions", {
  gr <- ggmethylation:::region_to_granges("chr2:1,000-2,000")
  expect_equal(GenomicRanges::start(gr), 1000L)
  expect_equal(GenomicRanges::end(gr), 2000L)
})

test_that("region_to_granges rejects malformed region", {
  expect_error(ggmethylation:::region_to_granges("chr1:1000"), "Invalid region format")
})
