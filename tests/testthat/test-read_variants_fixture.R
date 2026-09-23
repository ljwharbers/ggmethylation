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
