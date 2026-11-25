# classify_variant_row() for short variants, and read_variants() handling of
# list-typed INFO fields (Number=.) and spanning-deletion (*) ALT alleles.

# --- MNV / spanning deletion -------------------------------------------------

test_that("equal-length multi-base REF/ALT is an MNV, not a deletion", {
  expect_identical(ggmethylation:::classify_variant_row("AC", "GT", 100L)$type, "MNV")
  expect_identical(ggmethylation:::classify_variant_row("TTA", "GCA", 100L)$type, "MNV")
})

test_that("length-changing alleles are still insertion / deletion", {
  expect_identical(ggmethylation:::classify_variant_row("GTA", "GT", 100L)$type, "insertion")
  expect_identical(ggmethylation:::classify_variant_row("G", "GT", 100L)$type, "deletion")
  expect_identical(ggmethylation:::classify_variant_row("C", "A", 100L)$type, "SNV")
})

test_that("a spanning-deletion '*' ALT is not classified as a variant", {
  expect_identical(ggmethylation:::classify_variant_row("*", "A", 100L)$type, NA_character_)
})

# --- List-typed INFO fields ----------------------------------------------------

test_that(".info_first() flattens S4 List and base list INFO columns", {
  il = IRanges::IntegerList(-300L, integer(0), c(500L, 7L))
  expect_identical(ggmethylation:::.info_first(il, as.integer, NA_integer_),
                   c(-300L, NA_integer_, 500L))
  cl = IRanges::CharacterList("DEL", character(0), "BND")
  expect_identical(ggmethylation:::.info_first(cl, as.character, NA_character_),
                   c("DEL", NA_character_, "BND"))
  expect_identical(ggmethylation:::.info_first(list(1L, NULL), as.integer, NA_integer_),
                   c(1L, NA_integer_))
  expect_identical(ggmethylation:::.info_first(c(4L, NA), as.integer, NA_integer_),
                   c(4L, NA_integer_))
  expect_null(ggmethylation:::.info_first(NULL, as.integer, NA_integer_))
})

# A tiny VCF in the style of pbsv / SVIM / cuteSV: SVTYPE, END and SVLEN are
# declared Number=., so VariantAnnotation returns CharacterList /
# IntegerList columns. It also carries a spanning deletion (*) and an MNV.
.write_list_info_vcf = function() {
  skip_if_not_installed("VariantAnnotation")
  path = tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=100000>",
    "##INFO=<ID=SVTYPE,Number=.,Type=String,Description=\"t\">",
    "##INFO=<ID=END,Number=.,Type=Integer,Description=\"e\">",
    "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"l\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t1000\tsv1\tA\t<DEL>\t50\tPASS\tSVTYPE=DEL;SVLEN=-300\tGT\t0/1",
    "chr1\t2000\tsv2\tA\t<DUP>\t50\tPASS\tSVTYPE=DUP;END=2500;SVLEN=500\tGT\t0/1",
    "chr1\t3000\tsv3\tA\tA[chr7:5000[\t50\tPASS\tSVTYPE=BND\tGT\t0/1",
    "chr1\t4000\tmnv\tGT\tAC\t50\tPASS\t.\tGT\t0/1",
    "chr1\t4100\tspan\tA\t*\t50\tPASS\t.\tGT\t0/1"
  ), path)
  gz = Rsamtools::bgzip(path, overwrite = TRUE)
  Rsamtools::indexTabix(gz, format = "vcf")
  gz
}

test_that("read_variants handles Number=. SVTYPE/END/SVLEN and drops * ALTs", {
  gz = .write_list_info_vcf()
  vd = read_variants(gz, "chr1:1-10000")
  v  = vd$variants

  expect_identical(v$position, c(1000L, 2000L, 3000L, 4000L))
  expect_identical(v$type, c("DEL", "DUP", "BND", "MNV"))
  # DEL has no END: end = 1000 + abs(-300) = 1300
  expect_identical(v$end, c(1300L, 2500L, 3000L, 4000L))
  expect_identical(v$mate_chrom, c(NA, NA, "chr7", NA))
  expect_identical(v$mate_pos, c(NA, NA, 5000L, NA))
})
