# Tests for match_sa_to_vcf_bnd() in R/variant_overlay_bnd.R
# Internal function — call via ggmethylation:::match_sa_to_vcf_bnd()
#
# Signature:
#   match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)
#
# reads:  data.frame with columns read_name, sa_chrom, sa_pos (and optionally
#         others). sa_chrom / sa_pos are NA for non-SA reads.
# bnd_df: data.frame with columns position (int), mate_chrom (chr),
#         mate_pos (int). May be NULL or zero-row.
# tol:    integer position tolerance in bp (default 50L).
#
# Returns: reads with new logical column vcf_validated.
#
# A read's SA tag gives a coordinate on the SA chromosome, i.e. on the BND's
# *mate* side, so sa_pos must be compared with mate_pos. The helpers below put
# the BND's local position (chr1:6000) and its mate (chr7:1010) far apart and on
# different chromosomes so that a comparison against the wrong column fails.

# Helpers --------------------------------------------------------------------

make_reads = function(read_name = "r1",
                       sa_chrom  = "chr7",
                       sa_pos    = 1000L) {
  data.frame(
    read_name = read_name,
    sa_chrom  = sa_chrom,
    sa_pos    = as.integer(sa_pos),
    stringsAsFactors = FALSE
  )
}

make_bnd = function(position   = 6000L,
                     mate_chrom = "chr7",
                     mate_pos   = 1010L) {
  data.frame(
    position   = as.integer(position),
    mate_chrom = mate_chrom,
    mate_pos   = as.integer(mate_pos),
    stringsAsFactors = FALSE
  )
}

# --- NULL / empty bnd_df -------------------------------------------------------

test_that("match_sa_to_vcf_bnd returns vcf_validated FALSE for NULL bnd_df", {
  reads  = make_reads()
  result = ggmethylation:::match_sa_to_vcf_bnd(reads, NULL)

  expect_true("vcf_validated" %in% names(result))
  expect_false(any(result$vcf_validated))
})

test_that("match_sa_to_vcf_bnd returns vcf_validated FALSE for 0-row bnd_df", {
  reads  = make_reads()
  bnd_df = data.frame(
    position   = integer(0L),
    mate_chrom = character(0L),
    mate_pos   = integer(0L),
    stringsAsFactors = FALSE
  )
  result = ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df)

  expect_true("vcf_validated" %in% names(result))
  expect_false(any(result$vcf_validated))
})

# --- Matching against the mate coordinate -------------------------------------

test_that("match_sa_to_vcf_bnd sets vcf_validated TRUE when sa_pos is near mate_pos", {
  reads  = make_reads(sa_chrom = "chr7", sa_pos = 1000L)
  bnd_df = make_bnd(position = 6000L, mate_chrom = "chr7", mate_pos = 1010L)  # |1000 - 1010| = 10

  result = ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)

  expect_true(result$vcf_validated[1L])
})

test_that("match_sa_to_vcf_bnd ignores the BND's local position", {
  # sa_pos equals the local position exactly but is 4990 bp from the mate.
  reads  = make_reads(sa_chrom = "chr7", sa_pos = 6000L)
  bnd_df = make_bnd(position = 6000L, mate_chrom = "chr7", mate_pos = 1010L)

  result = ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)

  expect_false(result$vcf_validated[1L])
})

test_that("match_sa_to_vcf_bnd returns FALSE when chrom does not match", {
  reads  = make_reads(sa_chrom = "chr8", sa_pos = 1000L)
  bnd_df = make_bnd(mate_chrom = "chr7", mate_pos = 1010L)  # delta 10, wrong chrom

  result = ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)

  expect_false(result$vcf_validated[1L])
})

test_that("match_sa_to_vcf_bnd returns FALSE when position is out of tolerance", {
  reads  = make_reads(sa_chrom = "chr7", sa_pos = 1000L)
  bnd_df = make_bnd(mate_chrom = "chr7", mate_pos = 1100L)  # delta = 100 > 50

  result = ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)

  expect_false(result$vcf_validated[1L])
})

test_that("match_sa_to_vcf_bnd returns FALSE for reads with no SA tag", {
  reads = data.frame(
    read_name = c("r1", "r2"),
    sa_chrom  = c(NA_character_, NA_character_),
    sa_pos    = c(NA_integer_,   NA_integer_),
    stringsAsFactors = FALSE
  )
  bnd_df = make_bnd()

  result = ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)

  expect_false(any(result$vcf_validated))
})

test_that("match_sa_to_vcf_bnd handles exact boundary: |delta| == tol -> TRUE", {
  reads  = make_reads(sa_chrom = "chr7", sa_pos = 1000L)
  bnd_df = make_bnd(mate_chrom = "chr7", mate_pos = 1050L)  # delta = 50 == tol

  result = ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)

  expect_true(result$vcf_validated[1L])
})

test_that("match_sa_to_vcf_bnd handles exact boundary: |delta| == tol + 1 -> FALSE", {
  reads  = make_reads(sa_chrom = "chr7", sa_pos = 1000L)
  bnd_df = make_bnd(mate_chrom = "chr7", mate_pos = 1051L)  # delta = 51 > tol

  result = ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)

  expect_false(result$vcf_validated[1L])
})

test_that("match_sa_to_vcf_bnd tolerates BND rows whose mate could not be parsed", {
  reads  = make_reads(sa_chrom = "chr7", sa_pos = 1000L)
  bnd_df = rbind(
    make_bnd(mate_chrom = NA_character_, mate_pos = NA_integer_),
    make_bnd(mate_chrom = "chr7", mate_pos = 1010L)
  )

  result = ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)
  expect_true(result$vcf_validated[1L])

  result_na = ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df[1L, ], tol = 50L)
  expect_false(result_na$vcf_validated[1L])
})

test_that("match_sa_to_vcf_bnd always produces a logical vcf_validated column", {
  reads  = make_reads()
  bnd_df = make_bnd()

  result = ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df)

  expect_type(result$vcf_validated, "logical")
})

test_that("match_sa_to_vcf_bnd handles mixed SA and non-SA reads correctly", {
  reads = data.frame(
    read_name = c("r1", "r2", "r3"),
    sa_chrom  = c("chr7", NA_character_, "chr7"),
    sa_pos    = c(1000L,   NA_integer_,   2000L),
    stringsAsFactors = FALSE
  )
  # Only r1 should match (correct chrom + within tol of mate_pos 1010);
  # r2 has no SA; r3 is 990 bp from the mate.
  bnd_df = make_bnd(mate_chrom = "chr7", mate_pos = 1010L)

  result = ggmethylation:::match_sa_to_vcf_bnd(reads, bnd_df, tol = 50L)

  expect_true(result$vcf_validated[1L])   # r1: match
  expect_false(result$vcf_validated[2L])  # r2: no SA
  expect_false(result$vcf_validated[3L])  # r3: wrong position
})

test_that("SA reads on the fixture are validated against the BND mate", {
  skip_if_not_installed("VariantAnnotation")
  bam = testthat::test_path("fixtures", "synthetic.bam")
  vcf = testthat::test_path("fixtures", "synthetic.vcf.gz")
  md  = read_methylation(bam, "chr1:1000-7000")
  vd  = read_variants(vcf, "chr1:1000-7000")

  # The fixture BND is chr1:6000 -> chr7:5000. Six primary reads carry
  # SA:chr7,5000; five carry SA:chr3,12000 and one supplementary read carries
  # SA:chr1,4000, neither of which matches.
  bnd = vd$variants[vd$variants$type %in% "BND", , drop = FALSE]
  res = ggmethylation:::match_sa_to_vcf_bnd(md$reads, bnd, tol = vd$bnd_match_tol)
  validated = res$read_name[res$vcf_validated]
  expect_setequal(
    validated,
    c("read011", "read022", "read033", "read044", "read055", "read066")
  )
})

# --- Tests for build_bnd_layer() ---------------------------------------------

test_that("build_bnd_layer returns NULL for NULL input", {
  expect_null(ggmethylation:::build_bnd_layer(NULL))
})

test_that("build_bnd_layer returns a single mate-label layer for a one-row bnd_df", {
  bnd_df = data.frame(
    position   = 1000L,
    mate_chrom = "chr7",
    mate_pos   = 5000L,
    stringsAsFactors = FALSE
  )
  result = ggmethylation:::build_bnd_layer(bnd_df)

  expect_type(result, "list")
  expect_length(result, 1L)
})
