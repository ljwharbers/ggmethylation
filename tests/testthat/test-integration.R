# Integration tests against the committed modBAM fixture.
#
# Fixture: tests/testthat/fixtures/hg002_fiberseq_MEG3.bam
#   HG002 PacBio Fiber-seq, longphase-phased, aligned to CHM13, subset to the
#   MEG3 promoter/DMR.  See data-raw/build_test_fixture.sh for how it is built.
#
# MEG3 is imprinted, so the two haplotypes carry genuinely different
# methylation -- the grouped tests below assert that real difference rather
# than merely that the code ran.
#
# The exact counts below are regression anchors: the fixture is frozen, so a
# change in any of them means parsing behaviour changed and should be reviewed
# deliberately, not silently re-baselined.

bam_path <- testthat::test_path("fixtures", "hg002_fiberseq_MEG3.bam")
region   <- "chr14:95055000-95070000"

# MEG3 TSS in CHM13; the promoter DMR sits around it.
meg3_tss <- 95060996

test_that("the BAM fixture is present", {
  # This must never skip -- the fixture is committed.  If it fails, the
  # checkout is incomplete rather than the environment being unusual.
  expect_true(file.exists(bam_path))
})

test_that("read_methylation returns methylation_data (ungrouped)", {
  result <- read_methylation(bam_path, region)

  expect_s3_class(result, "methylation_data")
  expect_true(is.data.frame(result$reads))
  expect_true(is.data.frame(result$sites))
  expect_true(all(c("read_name", "start", "end", "strand") %in% names(result$reads)))
  expect_true(all(c("position", "mod_prob", "read_name", "mod_code") %in% names(result$sites)))

  expect_identical(nrow(result$reads), 52L)
  expect_true(all(result$sites$mod_prob >= 0 & result$sites$mod_prob <= 1, na.rm = TRUE))
})

test_that("mod_code selects between the two modifications in the fixture", {
  # The fixture carries both C+m (5mC) and A+a (6mA) on every read, so it
  # exercises mod_code filtering -- a path with no coverage before.
  m_5mc <- read_methylation(bam_path, region, mod_code = "m")
  m_6ma <- read_methylation(bam_path, region, mod_code = "a")

  expect_identical(unique(m_5mc$sites$mod_code), "m")
  expect_identical(unique(m_6ma$sites$mod_code), "a")

  expect_identical(nrow(m_5mc$sites), 10205L)
  expect_identical(nrow(m_6ma$sites), 77888L)

  # 6mA is far denser than 5mC in Fiber-seq data.
  expect_gt(nrow(m_6ma$sites), nrow(m_5mc$sites))
})

test_that("HP grouping splits reads by haplotype", {
  result <- read_methylation(bam_path, region, mod_code = "m",
                             group_tag = "HP", drop_na_group = FALSE)

  expect_s3_class(result, "methylation_data")
  expect_true("group" %in% names(result$reads))

  tab <- table(result$reads$group, useNA = "ifany")
  expect_identical(as.integer(tab[["1"]]), 30L)
  expect_identical(as.integer(tab[["2"]]), 21L)
  # Exactly one read in the window is unphased.
  expect_identical(sum(is.na(result$reads$group)), 1L)
})

test_that("drop_na_group removes the unphased read", {
  kept    <- read_methylation(bam_path, region, mod_code = "m",
                              group_tag = "HP", drop_na_group = FALSE)
  dropped <- read_methylation(bam_path, region, mod_code = "m",
                              group_tag = "HP", drop_na_group = TRUE)

  expect_identical(sum(is.na(kept$reads$group)), 1L)
  expect_identical(sum(is.na(dropped$reads$group)), 0L)
  expect_identical(nrow(dropped$reads), nrow(kept$reads) - 1L)
})

test_that("MEG3 promoter shows allele-specific methylation", {
  # The biological check: at an imprinted DMR one haplotype is methylated and
  # the other is not.  A flat result here means grouping or parsing is broken,
  # not that the locus is uninteresting.
  m <- read_methylation(bam_path, region, mod_code = "m",
                        group_tag = "HP", drop_na_group = TRUE)

  near <- m$sites[abs(m$sites$position - meg3_tss) < 2000, ]
  expect_gt(nrow(near), 1000L)

  means <- tapply(near$mod_prob, near$group, mean)
  expect_length(means, 2L)

  # Observed on this fixture: hap1 ~0.71, hap2 ~0.33.
  expect_gt(abs(means[["1"]] - means[["2"]]), 0.2)
})

test_that("plot_methylation returns a ggplot/patchwork for real data", {
  md <- read_methylation(bam_path, region, mod_code = "m")
  p  <- plot_methylation(md)
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})

test_that("plot_methylation renders the grouped delta and CI panels", {
  md <- read_methylation(bam_path, region, mod_code = "m",
                         group_tag = "HP", drop_na_group = TRUE)
  p  <- plot_methylation(md, show_delta = TRUE, show_ci = TRUE)

  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
  # reads + smooth + delta
  expect_gte(length(p$patches$plots) + 1L, 3L)
})

test_that("insertion-aware parsing populates insertion sites", {
  md <- read_methylation(bam_path, region, mod_code = "m")

  expect_true(is.data.frame(md$insertion_sites))
  expect_identical(nrow(md$insertion_sites), 54L)
  expect_true(all(c("read_name", "ref_anchor", "ins_offset", "ins_length",
                    "mod_prob", "mod_code") %in% names(md$insertion_sites)))
  expect_identical(md$insertion_sites, insertion_sites(md))
})
