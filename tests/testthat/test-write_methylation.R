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
