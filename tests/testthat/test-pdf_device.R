# Every label the package draws must work on the default pdf() device, which
# cannot draw characters outside its single-byte encoding (a "conversion
# failure ... in 'mbcsToSbcs'" error or warning). ggsave("x.pdf") uses it.

on_pdf = function(p) {
  f = tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) })
  print(p)
  invisible(TRUE)
}

test_that("the delta panel draws on pdf(), in continuous and binary mode", {
  md = suppressWarnings(read_methylation(fixture_bam(), "chr1:1000-7000",
                                         group_tag = "HP"))
  expect_no_warning(on_pdf(plot_methylation(md, show_delta = TRUE)))
  expect_no_warning(on_pdf(plot_methylation(md, show_delta = TRUE,
                                            call_threshold = 0.5)))
})

test_that("the variant overlay draws on pdf(), with an SV cut at the region edge", {
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("SummarizedExperiment")
  vars = read_variants(fixture_vcf(), "chr1:500-7000")
  expect_true(all(c("SNV", "DEL", "BND") %in% vars$variants$type))
  for (region in c("chr1:1000-7000", "chr1:3950-4050")) {
    md = suppressWarnings(read_methylation(fixture_bam(), region))
    expect_no_warning(on_pdf(plot_methylation(md, variants = vars)))
  }
})

test_that("the insertion locus plot draws on pdf()", {
  md = suppressWarnings(read_methylation(fixture_bam(), "chr1:1000-7000",
                                         mod_code = c("m", "h")))
  loci = list_insertion_loci(md, min_length = 50L)
  expect_gt(nrow(loci), 0L)
  expect_no_warning(on_pdf(plot_insertion_locus(md, loci[1L, ])))
})
