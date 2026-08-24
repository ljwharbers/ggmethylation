# Tests for read_variants() itself.
#
# Only its helpers (parse_bnd_alt, classify_variant_row) were covered before --
# read_variants() was never called in any test, leaving the VariantAnnotation
# VCF-reading path and INFO-field unwrapping untested. Both had real bugs fixed
# in 27a1e37 and 1c6fe29.

vcf_path <- testthat::test_path("fixtures", "hg002_MEG3_snvs.vcf.gz")
region   <- "chr14:95055000-95070000"

test_that("the VCF fixture is present", {
  expect_true(file.exists(vcf_path))
})

test_that("read_variants returns a variant_data object", {
  skip_if_not_installed("VariantAnnotation")

  v <- read_variants(vcf_path, region)

  expect_s3_class(v, "variant_data")
  expect_true(all(c("variants", "region") %in% names(v)))
  expect_true(is.data.frame(v$variants))
  expect_s4_class(v$region, "GRanges")
})

test_that("variants carry the expected columns and types", {
  skip_if_not_installed("VariantAnnotation")

  v <- read_variants(vcf_path, region)

  expect_true(all(c("position", "ref", "alt", "type") %in% names(v$variants)))
  expect_type(v$variants$position, "integer")
  expect_type(v$variants$ref, "character")
  expect_type(v$variants$alt, "character")

  # ALT must be unwrapped to a plain character vector, not left as the
  # DNAStringSetList / CharacterList that VariantAnnotation returns.
  expect_false(inherits(v$variants$alt, "List"))
  expect_length(v$variants$alt, nrow(v$variants))
  expect_false(any(grepl("^c\\(|list\\(", v$variants$alt)))
})

test_that("read_variants classifies the fixture's variants", {
  skip_if_not_installed("VariantAnnotation")

  v <- read_variants(vcf_path, region)

  # The fixture holds 15 phased heterozygous calls in the MEG3 window:
  # 13 biallelic SNVs and 2 short deletions.
  expect_identical(nrow(v$variants), 15L)
  expect_identical(sum(v$variants$type == "SNV"), 13L)
  expect_identical(sum(v$variants$type == "deletion"), 2L)

  snvs <- v$variants[v$variants$type == "SNV", ]
  expect_true(all(nchar(snvs$ref) == 1L))
  expect_true(all(nchar(snvs$alt) == 1L))
})

test_that("variants fall inside the requested region", {
  skip_if_not_installed("VariantAnnotation")

  v <- read_variants(vcf_path, region)

  expect_true(all(v$variants$position >= 95055000L))
  expect_true(all(v$variants$position <= 95070000L))
})

test_that("a narrower region returns a strict subset", {
  skip_if_not_installed("VariantAnnotation")

  wide   <- read_variants(vcf_path, region)
  narrow <- read_variants(vcf_path, "chr14:95066000-95067000")

  expect_lt(nrow(narrow$variants), nrow(wide$variants))
  expect_true(all(narrow$variants$position %in% wide$variants$position))
})

test_that("a region with no variants yields zero rows, not an error", {
  skip_if_not_installed("VariantAnnotation")

  v <- read_variants(vcf_path, "chr14:95055000-95055100")

  expect_s3_class(v, "variant_data")
  expect_identical(nrow(v$variants), 0L)
})

test_that("print.variant_data reports the type breakdown", {
  skip_if_not_installed("VariantAnnotation")

  v <- read_variants(vcf_path, region)

  expect_output(print(v), "variant_data object")
  expect_output(print(v), "Variants: 15 total")
  expect_output(print(v), "SNV: 13")
})

test_that("a missing VCF errors rather than returning empty", {
  expect_error(
    read_variants(file.path(tempdir(), "does_not_exist.vcf.gz"), region)
  )
})
