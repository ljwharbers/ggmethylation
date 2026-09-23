# Tests for extract_variant_bases() in R/variant_overlay.R
# Internal function — call via ggmethylation:::extract_variant_bases()
#
# Signature:
#   extract_variant_bases(reads, sequences, cigars, variants)
#
# reads    : data.frame with columns read_name, start, end, strand, lane
# sequences: named character vector keyed by read_name
# cigars   : named character vector keyed by read_name
# variants : data.frame with columns position, ref, alt, type
#
# Returns: data.frame with columns read_name, position, base,
#          variant_class, lane

# Helper: build reads data.frame
make_reads = function(read_name, start, end, lane = 1L, bam_pos = start) {
  data.frame(
    read_name = read_name,
    start     = start,
    end       = end,
    bam_pos   = bam_pos,
    strand    = "+",
    lane      = lane,
    stringsAsFactors = FALSE
  )
}

# Helper: build variants data.frame
make_variants = function(position, ref, alt, type = "SNV") {
  data.frame(
    position = as.integer(position),
    ref      = ref,
    alt      = alt,
    type     = type,
    stringsAsFactors = FALSE
  )
}

# --- Test 1: basic 16M CIGAR, variant at position 103, base matches ref ---

test_that("extract_variant_bases extracts correct base — ref classification", {
  # Read: starts at 100, 16M CIGAR -> covers 100..115
  # Sequence: 16 chars; position 103 is 0-indexed offset 3 -> query pos 4 (1-indexed)
  seq_str = "ACGTACGTACGTACGT"  # 16 bp; query pos 4 = 'T'
  reads    = make_reads("r1", 100L, 115L, lane = 1L)
  seqs     = c(r1 = seq_str)
  cigs     = c(r1 = "16M")
  variants = make_variants(position = 103L, ref = "T", alt = "A")

  result = ggmethylation:::extract_variant_bases(reads, seqs, cigs, variants)

  expect_equal(nrow(result), 1L)
  expect_equal(result$read_name, "r1")
  expect_equal(result$position, 103L)
  expect_equal(result$base, "T")
  expect_equal(result$variant_class, "ref")
  expect_equal(result$lane, 1L)
})

test_that("extract_variant_bases classifies alt correctly", {
  seq_str  = "ACGAACGTACGTACGT"  # position 103 = query pos 4 = 'A' (0-based offset 3)
  reads    = make_reads("r1", 100L, 115L, lane = 2L)
  seqs     = c(r1 = seq_str)
  cigs     = c(r1 = "16M")
  variants = make_variants(position = 103L, ref = "T", alt = "A")

  result = ggmethylation:::extract_variant_bases(reads, seqs, cigs, variants)

  expect_equal(nrow(result), 1L)
  expect_equal(result$variant_class, "alt")
  expect_equal(result$lane, 2L)
})

# --- Test 2: deletion in CIGAR, variant at deleted position -> "del" ---

test_that("extract_variant_bases classifies deleted position as del", {
  # 3M1D12M from ref 100:
  #   ref 100-102 -> query 1-3 (3M)
  #   ref 103     -> DELETED
  #   ref 104-115 -> query 4-15 (12M)
  seq_str  = "ACGACGTACGTABC"   # 15 bp (only 15 query bases because 1 deleted)
  reads    = make_reads("r1", 100L, 115L, lane = 1L)
  seqs     = c(r1 = seq_str)
  cigs     = c(r1 = "3M1D12M")
  variants = make_variants(position = 103L, ref = "T", alt = "A")

  result = ggmethylation:::extract_variant_bases(reads, seqs, cigs, variants)

  expect_equal(nrow(result), 1L)
  expect_equal(result$variant_class, "del")
  expect_equal(result$base, "-")
})

# --- Test 3: variant outside read span -> 0-row result ---

test_that("extract_variant_bases returns 0 rows when variant is outside read span", {
  reads    = make_reads("r1", 100L, 115L)
  seqs     = c(r1 = "ACGTACGTACGTACGT")
  cigs     = c(r1 = "16M")
  variants = make_variants(position = 200L, ref = "A", alt = "T")

  result = ggmethylation:::extract_variant_bases(reads, seqs, cigs, variants)

  expect_equal(nrow(result), 0L)
  expect_true(all(c("read_name", "position", "base", "variant_class", "lane")
                  %in% names(result)))
})

# --- Test 4: empty variants data.frame -> 0-row result ---

test_that("extract_variant_bases returns 0 rows for empty variants", {
  reads = make_reads("r1", 100L, 115L)
  seqs  = c(r1 = "ACGTACGTACGTACGT")
  cigs  = c(r1 = "16M")
  # Build a 0-row variants data.frame directly to avoid data.frame() recycling
  variants = data.frame(
    position = integer(0L),
    ref      = character(0L),
    alt      = character(0L),
    type     = character(0L),
    stringsAsFactors = FALSE
  )

  result = ggmethylation:::extract_variant_bases(reads, seqs, cigs, variants)

  expect_equal(nrow(result), 0L)
  expect_true(all(c("read_name", "position", "base", "variant_class", "lane")
                  %in% names(result)))
})

# --- Test 5: multiple reads, multiple variants -> correct row count ---

test_that("extract_variant_bases handles multiple reads and variants", {
  # Two reads, each spanning 100..119 (20 bp, 20M)
  # Two variants: at positions 105 and 110
  # Each read spans both variants -> 2 reads * 2 variants = 4 rows
  seq1 = "ACGTACGTACGTACGTACGT"  # 20 bp
  seq2 = "TTTTTTTTTTTTTTTTTTTT"  # 20 bp

  reads = rbind(
    make_reads("r1", 100L, 119L, lane = 1L),
    make_reads("r2", 100L, 119L, lane = 2L)
  )
  seqs = c(r1 = seq1, r2 = seq2)
  cigs = c(r1 = "20M", r2 = "20M")

  variants = rbind(
    make_variants(position = 105L, ref = "A", alt = "T"),
    make_variants(position = 110L, ref = "G", alt = "C")
  )

  result = ggmethylation:::extract_variant_bases(reads, seqs, cigs, variants)

  expect_equal(nrow(result), 4L)
  expect_equal(sort(unique(result$read_name)), c("r1", "r2"))
  expect_equal(sort(unique(result$position)), c(105L, 110L))
})

# --- Test 6: clipped start regression (issue #19) ---
# Read's true BAM POS is before the queried region, so reads$start is clipped.
# The CIGAR must be walked from bam_pos, not the clipped start.
#
# Layout: bam_pos=100, clipped start=140, variant at ref pos 143, 50M CIGAR.
#   correct q_pos = 1 + (143-100) = 44  -> seq_str[44] = 'G'
#   wrong   q_pos = 1 + (143-140) =  4  -> seq_str[4]  = 'A'
# ref='G' alt='A': correct gives "ref", wrong gives "alt" -- discriminating.
test_that("extract_variant_bases uses bam_pos for CIGAR walk when start is clipped", {
  seq_str  = "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAA"  # 50 bp
  bam_pos  = 100L
  clipped  = 140L

  # seq_str[44] = 'G' (within GGGG block 41-44), seq_str[4] = 'A' (AAAA block 1-4)
  correct_base = substr(seq_str, 44L, 44L)  # 'G' — base at correct q_pos 44
  wrong_base   = substr(seq_str, 4L,  4L)   # 'A' — base at wrong q_pos 4
  expect_false(correct_base == wrong_base)    # confirm they differ

  reads    = make_reads("r1", start = clipped, end = 149L, lane = 1L,
                         bam_pos = bam_pos)
  seqs     = c(r1 = seq_str)
  cigs     = c(r1 = "50M")
  variants = make_variants(position = 143L, ref = correct_base, alt = wrong_base)

  result = ggmethylation:::extract_variant_bases(reads, seqs, cigs, variants)

  expect_equal(nrow(result), 1L)
  expect_equal(result$base, correct_base)
  expect_equal(result$variant_class, "ref")
})
