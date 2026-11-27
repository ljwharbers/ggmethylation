# print() methods for methylation_data, multi_methylation_data and
# variant_data, on the synthetic fixture.

print_fixture_md = function(...) {
  read_methylation(fixture_bam(), "chr1:1000-7000", seed = 1, ...)
}

test_that("print.methylation_data reports region, reads, codes and groups", {
  md = print_fixture_md(group_tag = "HP")
  expect_output(print(md), "^methylation_data object")
  expect_output(print(md), "Region: chr1:1000-7000")
  expect_output(print(md), "Reads: [0-9]+ +\\(\\+ strand: [0-9]+, +- strand: [0-9]+\\)")
  expect_output(print(md), "Median aligned length: [0-9]+ bp")
  expect_output(print(md), "Modification code\\(s\\): m")
  expect_output(print(md), "Group tag: HP \\(2 groups\\)")
  utils::capture.output(out <- withVisible(print(md)))
  expect_false(out$visible)
  expect_identical(out$value, md)
})

test_that("print.methylation_data of an empty object prints no aligned length", {
  md = empty_methylation_data(region_to_granges("chr1:1-100"), "m", NULL)
  expect_output(print(md), "Reads: 0")
  out = utils::capture.output(print(md))
  expect_false(any(grepl("aligned length", out)))
})

test_that("print.methylation_data lists several modification codes", {
  md = print_fixture_md(mod_code = c("m", "h"))
  expect_output(print(md), "Modification code\\(s\\): m, h")
})

test_that("print.multi_methylation_data reports samples and per-sample reads", {
  md = print_fixture_md(group_tag = "HP")
  merged = merge_methylation(A = md, B = md)
  expect_output(print(merged), "^multi_methylation_data object")
  expect_output(print(merged), "Samples: +2")
  expect_output(print(merged), "Mod code: m")
  expect_output(print(merged), "A: [0-9]+ reads \\[1=[0-9]+, 2=[0-9]+")
})

test_that("print.variant_data reports region and counts per type", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("SummarizedExperiment")
  vd = read_variants(fixture_vcf(), "chr1:1000-7000")
  expect_output(print(vd), "^variant_data object")
  expect_output(print(vd), "Region: chr1:1000-7000")
  expect_output(print(vd), "Variants: 3 total")
  expect_output(print(vd), "SNV: 1")
  expect_output(print(vd), "BND: 1")
  expect_output(print(vd), "DEL: 1")
})
