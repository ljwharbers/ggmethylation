# read_variants() against the synthetic fixture VCF.

test_that("read_variants stores bnd_match_tol for the SA overlay", {
  skip_if_not_installed("VariantAnnotation")
  vcf = testthat::test_path("fixtures", "synthetic.vcf.gz")
  vd = read_variants(vcf, "chr1:1000-7000")
  expect_identical(vd$bnd_match_tol, 50L)
  expect_setequal(vd$variants$type, c("SNV", "DEL", "BND"))
  expect_identical(
    read_variants(vcf, "chr1:1000-7000", bnd_match_tol = 5)$bnd_match_tol, 5L
  )
})

test_that("read_variants returns the fixture variants exactly", {
  skip_if_not_installed("VariantAnnotation")
  vcf = testthat::test_path("fixtures", "synthetic.vcf.gz")
  v = read_variants(vcf, "chr1:1000-7000")$variants
  expect_identical(v, data.frame(
    position   = c(3000L, 3899L, 6000L),
    ref        = c("A", "A", "A"),
    alt        = c("G", "<DEL>", "A[chr7:5000["),
    type       = c("SNV", "DEL", "BND"),
    end        = c(3000L, 4050L, 6000L),
    mate_chrom = c(NA, NA, "chr7"),
    mate_pos   = c(NA, NA, 5000L),
    stringsAsFactors = FALSE
  ))
})

test_that("read_variants validates bnd_match_tol", {
  skip_if_not_installed("VariantAnnotation")
  vcf = testthat::test_path("fixtures", "synthetic.vcf.gz")
  for (bad in list(-1, c(10, 20), NA_real_, Inf, "50", NULL)) {
    expect_error(read_variants(vcf, "chr1:1000-7000", bnd_match_tol = bad),
                 "bnd_match_tol")
  }
  expect_identical(read_variants(vcf, "chr1:1000-7000", bnd_match_tol = 0)$bnd_match_tol, 0L)
  # Non-whole values are still accepted (stored truncated, as before)
  expect_identical(read_variants(vcf, "chr1:1000-7000", bnd_match_tol = 12.5)$bnd_match_tol, 12L)
  expect_error(read_variants(vcf, "chr1:1000-7000", bnd_match_tol = 3e9), "bnd_match_tol")
})

test_that("read_variants validates vcf", {
  skip_if_not_installed("VariantAnnotation")
  vcf = testthat::test_path("fixtures", "synthetic.vcf.gz")
  for (bad in list(NULL, NA_character_, "", c(vcf, vcf), 1)) {
    expect_error(read_variants(bad, "chr1:1000-7000"),
                 "`vcf` must be a single non-empty character string")
  }
  expect_error(read_variants(file.path(tempdir(), "absent.vcf.gz"), "chr1:1000-7000"),
               "`vcf` file not found")
})

test_that("read_variants warns when the region's chromosome is not in the VCF", {
  skip_if_not_installed("VariantAnnotation")
  vcf = testthat::test_path("fixtures", "synthetic.vcf.gz")
  # The fixture uses "chr1"; "1" is the Ensembl-style name of the same contig.
  expect_warning(vd <- read_variants(vcf, "1:1000-7000"), "chr1' vs '1")
  expect_s3_class(vd, "variant_data")
  expect_equal(nrow(vd$variants), 0L)
  # A declared contig without records in the region is silently empty
  expect_silent(vd7 <- read_variants(vcf, "chr7:1-1000"))
  expect_equal(nrow(vd7$variants), 0L)
})

test_that("read_variants reads a VCF without INFO header lines", {
  skip_if_not_installed("VariantAnnotation")
  path = tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tG\t50\tPASS\t.",
    "chr1\t200\t.\tAC\tGT\t50\tPASS\t."
  ), path)
  gz = Rsamtools::bgzip(path, overwrite = TRUE)
  Rsamtools::indexTabix(gz, format = "vcf")
  expect_silent(vd <- read_variants(gz, "chr1:1-1000"))
  expect_identical(vd$variants$type, c("SNV", "MNV"))
  # No ##contig lines: the tabix index still knows "chr1", but not "1"
  expect_warning(read_variants(gz, "1:1-1000"), "not among the VCF's contigs")
})
