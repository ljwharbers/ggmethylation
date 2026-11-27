# write_methylation(): TSV and BED export.

# Two reads (r1 on +, r2 on -) with sites given out of coordinate order and
# with two modification codes.
make_write_data = function() {
  md = empty_methylation_data(region_to_granges("chr1:99990-100010"), c("m", "h"), NULL)
  md$reads = data.frame(
    read_name = c("r1", "r2"), start = c(99990L, 99990L), end = c(100010L, 100010L),
    bam_pos = c(99990L, 99990L), strand = c("+", "-"), is_supplementary = FALSE,
    sa_chrom = NA_character_, sa_pos = NA_integer_, clip_side = NA_character_,
    stringsAsFactors = FALSE
  )
  md$sites = data.frame(
    position  = c(100000L, 100005L, 99995L, 99995L),
    mod_prob  = c(0.8567, 0.5, 1, 0.0004),
    read_name = c("r1", "r2", "r2", "r1"),
    mod_code  = c("m", "m", "h", "m"),
    stringsAsFactors = FALSE
  )
  md
}

test_that("BED export has integer scores, mod_code in the name and sorted rows", {
  prefix = file.path(tempfile("bedout"), "x")
  write_methylation(make_write_data(), prefix, format = "bed", reads = FALSE)
  lines = readLines(paste0(prefix, "_sites.bed"))
  expect_identical(lines, c(
    "chr1\t99994\t99995\tr1:m\t0\t+",
    "chr1\t99994\t99995\tr2:h\t1000\t-",
    "chr1\t99999\t100000\tr1:m\t857\t+",
    "chr1\t100004\t100005\tr2:m\t500\t-"
  ))
})

test_that("BED export works gzip-compressed", {
  prefix = file.path(tempfile("bedgz"), "x")
  write_methylation(make_write_data(), prefix, format = "bed", reads = FALSE, gzip = TRUE)
  lines = readLines(gzfile(paste0(prefix, "_sites.bed.gz")))
  expect_length(lines, 4L)
  expect_identical(lines[3], "chr1\t99999\t100000\tr1:m\t857\t+")
})

test_that("an object without sites writes an empty BED and a TSV", {
  md = empty_methylation_data(region_to_granges("chr1:1-100"), "m", NULL)
  prefix = file.path(tempfile("empty"), "x")
  write_methylation(md, prefix, format = "bed")
  expect_length(readLines(paste0(prefix, "_sites.bed")), 0L)
  write_methylation(md, prefix, format = "tsv")
  expect_length(readLines(paste0(prefix, "_sites.tsv")), 1L)  # header only
})

test_that("reads TSV mean_mod_prob pools every modification code of a read", {
  prefix = file.path(tempfile("tsvout"), "x")
  write_methylation(make_write_data(), prefix, sites = FALSE)
  reads = utils::read.delim(paste0(prefix, "_reads.tsv"))
  # r1: m 0.8567 and m 0.0004; r2: m 0.5 and h 1.
  expect_equal(reads$mean_mod_prob, c((0.8567 + 0.0004) / 2, (0.5 + 1) / 2))
})

test_that("write_methylation validates its arguments, naming them", {
  md = make_write_data()
  prefix = file.path(tempfile("bad"), "x")
  expect_error(write_methylation(md, prefix, reads = NA), "`reads` must be TRUE or FALSE")
  expect_error(write_methylation(md, prefix, sites = "yes"), "`sites` must be TRUE or FALSE")
  expect_error(write_methylation(md, prefix, gzip = c(TRUE, FALSE)), "`gzip` must be TRUE or FALSE")
  expect_error(write_methylation(md, NA_character_), "`prefix` must be")
  expect_error(write_methylation(md, ""), "`prefix` must be")
  expect_error(write_methylation(md$sites, prefix), "`data` must be")
  expect_false(dir.exists(dirname(prefix)))

  expect_warning(write_methylation(md, prefix, reads = FALSE, sites = FALSE),
                 "Nothing to write")
  expect_false(dir.exists(dirname(prefix)))
})

# Second sample: only r1 and its sites.
make_multi_write_data = function() {
  md_a = make_write_data()
  md_b = make_write_data()
  md_b$reads = md_b$reads[md_b$reads$read_name == "r1", , drop = FALSE]
  md_b$sites = md_b$sites[md_b$sites$read_name == "r1", , drop = FALSE]
  merge_methylation(A = md_a, B = md_b)
}

test_that("write_methylation writes one file set per sample of a multi_methylation_data", {
  multi = make_multi_write_data()
  dir = tempfile("multi")
  prefix = file.path(dir, "run")

  expect_invisible(out <- write_methylation(multi, prefix))
  expect_identical(out, multi)
  expect_setequal(list.files(dir), c("run_A_reads.tsv", "run_A_sites.tsv",
                                     "run_B_reads.tsv", "run_B_sites.tsv"))
  expect_equal(nrow(utils::read.delim(file.path(dir, "run_A_reads.tsv"))), 2L)
  expect_equal(nrow(utils::read.delim(file.path(dir, "run_B_reads.tsv"))), 1L)
  expect_equal(nrow(utils::read.delim(file.path(dir, "run_A_sites.tsv"))), 4L)
  expect_equal(nrow(utils::read.delim(file.path(dir, "run_B_sites.tsv"))), 2L)

  write_methylation(multi, prefix, format = "bed", reads = FALSE)
  expect_length(readLines(file.path(dir, "run_A_sites.bed")), 4L)
  expect_length(readLines(file.path(dir, "run_B_sites.bed")), 2L)
})

test_that("write_methylation writes gzipped per-sample files", {
  multi = make_multi_write_data()
  dir = tempfile("multigz")
  prefix = file.path(dir, "run")

  write_methylation(multi, prefix, gzip = TRUE)
  write_methylation(multi, prefix, format = "bed", reads = FALSE, gzip = TRUE)
  expect_setequal(list.files(dir), c(
    "run_A_reads.tsv.gz", "run_A_sites.tsv.gz", "run_A_sites.bed.gz",
    "run_B_reads.tsv.gz", "run_B_sites.tsv.gz", "run_B_sites.bed.gz"
  ))
  expect_equal(nrow(utils::read.delim(gzfile(file.path(dir, "run_B_reads.tsv.gz")))), 1L)
  expect_equal(nrow(utils::read.delim(gzfile(file.path(dir, "run_A_sites.tsv.gz")))), 4L)
  expect_length(readLines(gzfile(file.path(dir, "run_A_sites.bed.gz"))), 4L)
  expect_length(readLines(gzfile(file.path(dir, "run_B_sites.bed.gz"))), 2L)
})

test_that("write_methylation rejects sample names containing a path separator", {
  md = make_write_data()
  multi = merge_methylation(.list = list(`a/b` = md, ok = md))
  dir = tempfile("sep")
  expect_error(write_methylation(multi, file.path(dir, "run")),
               "must not contain a path separator: 'a/b'")
  expect_false(dir.exists(dir))
})
