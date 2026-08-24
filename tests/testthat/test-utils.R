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

# --- .query_extent ---
# Extents are 0-based, in ORIGINAL read orientation, and count hard clips
# towards the read length so a primary and its supplementary counterpart share
# one axis.

test_that(".query_extent places a forward alignment after its leading clip", {
  expect_equal(
    ggmethylation:::.query_extent("20S100M400S", "+"),
    list(start = 20L, end = 119L, qlen = 520L)
  )
})

test_that(".query_extent mirrors the offsets for a reverse alignment", {
  # Stored reference-oriented, so the TRAILING clip is the read-5' offset.
  expect_equal(
    ggmethylation:::.query_extent("20S100M400S", "-"),
    list(start = 400L, end = 499L, qlen = 520L)
  )
})

test_that(".query_extent counts hard clips towards read length", {
  expect_equal(
    ggmethylation:::.query_extent("100H400M", "+"),
    list(start = 100L, end = 499L, qlen = 500L)
  )
})

test_that(".query_extent returns NULL for unusable CIGARs", {
  expect_null(ggmethylation:::.query_extent("*", "+"))
  expect_null(ggmethylation:::.query_extent(NA_character_, "+"))
  expect_null(ggmethylation:::.query_extent("500S", "+"))
})

# --- sa_partner_sides ---
# Which REFERENCE flank of this alignment each SA partner joins, derived from
# read coordinates. Test CIGAR/SA pairs are kept physically consistent: the two
# extents of a real chimeric read never overlap.

test_that("sa_partner_sides puts a downstream partner on the right (+ strand)", {
  res <- ggmethylation:::sa_partner_sides(
    "100M400S", "+", "chr7,500,+,100H400M,60,0"
  )
  expect_equal(res$side, "right")
  expect_equal(res$rname, "chr7")
  expect_equal(res$pos, 500L)
})

test_that("sa_partner_sides mirrors the flank on the reverse strand", {
  # Primary covers read bases 0..99, which for a reverse alignment sit at the
  # right reference edge, so the downstream partner joins on the left.
  res <- ggmethylation:::sa_partner_sides(
    "400S100M", "-", "chr7,500,+,100H400M,60,0"
  )
  expect_equal(res$side, "left")
})

test_that("sa_partner_sides puts an upstream partner on the left (+ strand)", {
  res <- ggmethylation:::sa_partner_sides(
    "400S100M", "+", "chr7,500,+,400M100H,60,0"
  )
  expect_equal(res$side, "left")
})

test_that("sa_partner_sides handles a reverse-strand partner", {
  # The partner is stored reverse-complemented, so its trailing hard clip is
  # the read-5' offset: it covers read 100..499, i.e. downstream.
  res <- ggmethylation:::sa_partner_sides(
    "100M400S", "+", "chr7,500,-,400M100H,60,0"
  )
  expect_equal(res$side, "right")
})

test_that("sa_partner_sides gives one side when both ends are clipped", {
  # The regression this helper exists for: adapter trimming clips both ends,
  # so detect_clip_side() reports "both" for a read with a single partner.
  cig <- "20S100M400S"
  expect_equal(ggmethylation:::detect_clip_side(cig), "both")
  res <- ggmethylation:::sa_partner_sides(cig, "+", "chr7,500,+,120H400M,60,0")
  expect_equal(res$side, "right")
})

test_that("sa_partner_sides reports both flanks for two partners", {
  res <- ggmethylation:::sa_partner_sides(
    "200S100M200S", "+",
    "chr7,500,+,200M300H,60,0;chr9,900,+,300H200M,55,0"
  )
  expect_equal(res$side, c("left", "right"))
  expect_equal(res$rname, c("chr7", "chr9"))
})

test_that("sa_partner_sides returns zero rows when there is no SA tag", {
  expect_equal(
    nrow(ggmethylation:::sa_partner_sides("100M400S", "+", NA_character_)), 0L
  )
  expect_equal(
    nrow(ggmethylation:::sa_partner_sides("100M400S", "+", "")), 0L
  )
})

test_that("sa_columns records the flank and its partner", {
  out <- ggmethylation:::sa_columns(
    cigar   = c("100M400S", "400S100M", "500M"),
    strand  = c("+", "+", "+"),
    sa_tags = list("chr7,500,+,100H400M,60,0",
                   "chr9,900,+,400M100H,60,0",
                   NA_character_)
  )
  expect_equal(out$sa_side, c("right", "left", NA))
  expect_equal(out$sa_chrom_right, c("chr7", NA, NA))
  expect_equal(out$sa_chrom_left,  c(NA, "chr9", NA))
  expect_equal(out$sa_pos_right,   c(500L, NA, NA))
  # sa_chrom / sa_pos stay populated for BND matching.
  expect_equal(out$sa_chrom, c("chr7", "chr9", NA))
})

test_that("sa_columns keeps the two partners of a double breakpoint apart", {
  out <- ggmethylation:::sa_columns(
    cigar   = "200S100M200S",
    strand  = "+",
    sa_tags = list("chr7,500,+,200M300H,60,0;chr9,900,+,300H200M,55,0")
  )
  expect_equal(out$sa_side, "both")
  expect_equal(out$sa_chrom_left,  "chr7")
  expect_equal(out$sa_chrom_right, "chr9")
})

test_that("sa_columns picks the highest-MAPQ partner on a flank", {
  out <- ggmethylation:::sa_columns(
    cigar   = "100M400S",
    strand  = "+",
    sa_tags = list("chr7,500,+,100H400M,20,0;chr9,900,+,100H400M,60,0")
  )
  expect_equal(out$sa_side, "right")
  expect_equal(out$sa_chrom_right, "chr9")
})

test_that("sa_columns leaves sa_side NA but keeps sa_chrom when undecidable", {
  out <- ggmethylation:::sa_columns(
    cigar = "100M400S", strand = "+",
    sa_tags = list("chr7,500,+,100M400H,60,0")   # exact tie in read space
  )
  expect_true(is.na(out$sa_side))
  expect_equal(out$sa_chrom, "chr7")
})

test_that("sa_columns returns all-NA columns when the BAM has no SA tag", {
  out <- ggmethylation:::sa_columns(c("100M", "100M"), c("+", "-"), NULL)
  expect_equal(nrow(out), 2L)
  expect_true(all(vapply(out, function(x) all(is.na(x)), logical(1L))))
})

test_that("sa_partner_sides keeps the entry but drops the side when undecidable", {
  # Unparseable primary, unparseable entry, and an exact tie in read space all
  # leave `side` as NA rather than guessing a flank.
  expect_true(is.na(
    ggmethylation:::sa_partner_sides("*", "+", "chr7,500,+,100H400M,60,0")$side
  ))
  expect_true(is.na(
    ggmethylation:::sa_partner_sides("100M400S", "+", "chr7,500,+,*,60,0")$side
  ))
  expect_true(is.na(
    ggmethylation:::sa_partner_sides("100M400S", "+", "chr7,500,+,100M400H,60,0")$side
  ))
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
