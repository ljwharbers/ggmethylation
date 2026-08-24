# Tests for write_methylation() -- TSV, BED, gzip, and the no-op guard.
#
# write_methylation() is exported and documented but had no coverage at all
# before this file: `write_methylation(` appeared nowhere under tests/.

make_write_md <- function(grouped = FALSE) {
  reads <- data.frame(
    read_name = c("r1", "r2"),
    start     = c(1000L, 1200L),
    end       = c(1800L, 2000L),
    strand    = c("+", "-"),
    stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = c("r1", "r1", "r2", "r2"),
    position  = c(1100L, 1200L, 1300L, 1400L),
    mod_prob  = c(0.9, 0.1, 0.4, 0.6),
    mod_code  = "m",
    stringsAsFactors = FALSE
  )
  if (grouped) {
    reads$group <- c("1", "2")
    sites$group <- c("1", "1", "2", "2")
  }

  structure(
    list(
      reads           = reads,
      sites           = sites,
      insertion_sites = sites[0, , drop = FALSE],
      region          = GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(1000, 2000)
      ),
      mod_code        = "m",
      group_tag       = if (grouped) "HP" else NULL
    ),
    class = "methylation_data"
  )
}

test_that("write_methylation writes both TSV files and returns data invisibly", {
  md     <- make_write_md()
  prefix <- file.path(withr::local_tempdir(), "out")

  res <- withVisible(write_methylation(md, prefix = prefix))
  expect_false(res$visible)
  expect_s3_class(res$value, "methylation_data")

  expect_true(file.exists(paste0(prefix, "_reads.tsv")))
  expect_true(file.exists(paste0(prefix, "_sites.tsv")))

  reads_out <- read.delim(paste0(prefix, "_reads.tsv"), stringsAsFactors = FALSE)
  expect_identical(nrow(reads_out), 2L)
  expect_identical(names(reads_out),
                   c("read_name", "start", "end", "strand", "mean_mod_prob"))
  # r1 sites are 0.9 and 0.1 -> mean 0.5
  expect_equal(reads_out$mean_mod_prob[reads_out$read_name == "r1"], 0.5)

  sites_out <- read.delim(paste0(prefix, "_sites.tsv"), stringsAsFactors = FALSE)
  expect_identical(nrow(sites_out), 4L)
  expect_identical(names(sites_out),
                   c("position", "mod_prob", "read_name", "mod_code"))
})

test_that("the group column is written only when the object is grouped", {
  prefix_u <- file.path(withr::local_tempdir(), "ungrouped")
  prefix_g <- file.path(withr::local_tempdir(), "grouped")

  write_methylation(make_write_md(grouped = FALSE), prefix = prefix_u)
  write_methylation(make_write_md(grouped = TRUE),  prefix = prefix_g)

  ungrouped <- read.delim(paste0(prefix_u, "_reads.tsv"), stringsAsFactors = FALSE)
  grouped   <- read.delim(paste0(prefix_g, "_reads.tsv"), stringsAsFactors = FALSE)

  expect_false("group" %in% names(ungrouped))
  expect_true("group"  %in% names(grouped))
  expect_identical(as.character(grouped$group), c("1", "2"))
})

test_that("reads = FALSE and sites = FALSE each suppress one file", {
  dir <- withr::local_tempdir()

  p_reads <- file.path(dir, "readsonly")
  write_methylation(make_write_md(), prefix = p_reads, sites = FALSE)
  expect_true(file.exists(paste0(p_reads, "_reads.tsv")))
  expect_false(file.exists(paste0(p_reads, "_sites.tsv")))

  p_sites <- file.path(dir, "sitesonly")
  write_methylation(make_write_md(), prefix = p_sites, reads = FALSE)
  expect_false(file.exists(paste0(p_sites, "_reads.tsv")))
  expect_true(file.exists(paste0(p_sites, "_sites.tsv")))
})

test_that("writing neither reads nor sites warns and writes nothing", {
  dir    <- withr::local_tempdir()
  prefix <- file.path(dir, "nothing")

  expect_warning(
    write_methylation(make_write_md(), prefix = prefix,
                      reads = FALSE, sites = FALSE),
    "Nothing to write"
  )
  expect_length(list.files(dir), 0L)
})

test_that("BED output is 0-based half-open with mod_prob as the score", {
  prefix <- file.path(withr::local_tempdir(), "bed")
  write_methylation(make_write_md(), prefix = prefix, format = "bed")

  bed_path <- paste0(prefix, "_sites.bed")
  expect_true(file.exists(bed_path))
  # BED carries no header.
  expect_false(file.exists(paste0(prefix, "_sites.tsv")))

  bed <- read.delim(bed_path, header = FALSE, stringsAsFactors = FALSE)
  expect_identical(ncol(bed), 6L)
  expect_identical(nrow(bed), 4L)

  names(bed) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")
  expect_true(all(bed$chrom == "chr1"))
  # position 1100 (1-based) -> [1099, 1100) 0-based half-open
  expect_identical(bed$chromStart, c(1099L, 1199L, 1299L, 1399L))
  expect_identical(bed$chromEnd,   c(1100L, 1200L, 1300L, 1400L))
  expect_identical(bed$chromEnd - bed$chromStart, rep(1L, 4L))
  expect_equal(bed$score, c(0.9, 0.1, 0.4, 0.6))
  # strand comes from the read the site belongs to
  expect_identical(bed$strand, c("+", "+", "-", "-"))
})

test_that("gzip appends .gz and the content round-trips", {
  prefix <- file.path(withr::local_tempdir(), "gz")
  write_methylation(make_write_md(), prefix = prefix, gzip = TRUE)

  expect_true(file.exists(paste0(prefix, "_reads.tsv.gz")))
  expect_true(file.exists(paste0(prefix, "_sites.tsv.gz")))
  expect_false(file.exists(paste0(prefix, "_reads.tsv")))

  sites <- read.delim(gzfile(paste0(prefix, "_sites.tsv.gz")),
                      stringsAsFactors = FALSE)
  expect_identical(nrow(sites), 4L)
  expect_equal(sites$mod_prob, c(0.9, 0.1, 0.4, 0.6))
})

test_that("gzip works for BED output too", {
  prefix <- file.path(withr::local_tempdir(), "bedgz")
  write_methylation(make_write_md(), prefix = prefix,
                    format = "bed", gzip = TRUE)

  expect_true(file.exists(paste0(prefix, "_sites.bed.gz")))
  bed <- read.delim(gzfile(paste0(prefix, "_sites.bed.gz")),
                    header = FALSE, stringsAsFactors = FALSE)
  expect_identical(ncol(bed), 6L)
  expect_identical(nrow(bed), 4L)
})

test_that("intermediate directories are created", {
  nested <- file.path(withr::local_tempdir(), "a", "b", "c", "out")
  expect_false(dir.exists(dirname(nested)))

  write_methylation(make_write_md(), prefix = nested)

  expect_true(dir.exists(dirname(nested)))
  expect_true(file.exists(paste0(nested, "_reads.tsv")))
})

test_that("input is validated", {
  prefix <- file.path(withr::local_tempdir(), "bad")

  expect_error(write_methylation(list(), prefix = prefix),
               "methylation_data")
  expect_error(write_methylation(make_write_md(), prefix = ""),
               "non-empty character")
  expect_error(write_methylation(make_write_md(), prefix = prefix,
                                 format = "vcf"),
               "'arg' should be one of")
})
