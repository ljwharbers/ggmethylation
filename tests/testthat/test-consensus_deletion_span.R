# A consensus deletion is judged against the reads that span the whole
# deletion, not every read of the group in the window or every read that
# merely overlaps it.

test_that("carriers are compared with reads spanning the whole deletion", {
  # Deletion 2000-2100 carried by 3 reads. One more read spans it without the
  # deletion; four reads only reach into it or lie elsewhere in the window.
  reads = data.frame(
    read_name = paste0("r", 1:8),
    start = c(1000L, 1000L, 1000L, 1000L, 1000L, 2050L, 3000L, 3000L),
    end   = c(3000L, 3000L, 3000L, 3000L, 2050L, 3000L, 4000L, 4000L),
    group = "1",
    stringsAsFactors = FALSE
  )
  dels = data.frame(
    type = "D", ref_start = 2000L, ref_end = 2100L,
    read_name = c("r1", "r2", "r3"), stringsAsFactors = FALSE
  )
  # 3 carriers of 4 spanning reads = 0.75: a consensus deletion. Of all 8 reads
  # (3/8) or of the 6 overlapping reads (3/6) it would not be.
  out = .consensus_deletion_ranges(dels, reads, "group")
  expect_equal(nrow(out), 1L)
  expect_equal(out$del_start, 2000L)
  expect_equal(out$del_end, 2100L)

  # A second spanning read without the deletion: 3 of 5 = 0.6, below 0.75.
  reads$end[5L] = 3000L
  expect_equal(nrow(.consensus_deletion_ranges(dels, reads, "group")), 0L)
})

test_that("the fixture's shared deletion breaks the group 1 smooth line", {
  md = suppressWarnings(
    read_methylation(fixture_bam(), "chr1:1000-7000", group_tag = "HP")
  )
  out = .consensus_deletion_ranges(
    .large_deletions(md$cigar_features, 50L), md$reads, "group"
  )
  expect_equal(out$del_start, 3900L)
  expect_equal(out$del_end, 4049L)
  expect_equal(out$group, "1")
})
