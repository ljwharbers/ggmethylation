# Build the vignette cache.
#
# Run this from the package root:
#   Rscript data-raw/build_vignette_cache.R
#
# Unlike the previous version, this builds entirely from the committed test
# fixtures, so anyone with a checkout can regenerate the cache -- no cluster
# paths, no vanished BAMs.  The only external dependency is the CHM13
# ncbiRefSeq GTF that read_annotations() downloads and caches on first use.
#
# The locus is MEG3, an imprinted lncRNA whose promoter DMR carries
# allele-specific methylation.  Every grouped figure in the vignette therefore
# shows a real biological difference rather than noise.

devtools::load_all()

bam    <- "tests/testthat/fixtures/hg002_fiberseq_MEG3.bam"
vcf    <- "tests/testthat/fixtures/hg002_MEG3_snvs.vcf.gz"
region <- "chr14:95055000-95070000"

stopifnot(file.exists(bam), file.exists(vcf))

message("Building meth_basic (ungrouped, 5mC) ...")
meth_basic <- read_methylation(bam, region, mod_code = "m")

message("Building meth_6ma (ungrouped, 6mA) ...")
meth_6ma <- read_methylation(bam, region, mod_code = "a")

message("Building meth_hp (haplotype-grouped) ...")
meth_hp <- read_methylation(bam, region, mod_code = "m",
                            group_tag = "HP", drop_na_group = TRUE)

# A real heterozygous SNV inside the window, taken from the VCF fixture.
# The previous version of this script used `meth_snv <- meth_hp` as a
# placeholder, which made the vignette's SNV section silently untrue.
message("Building meth_snv (grouped by SNV genotype at chr14:95066012 C>T) ...")
meth_snv <- read_methylation(bam, region, mod_code = "m",
                             snv_position = 95066012L,
                             ref_base     = "C",
                             alt_base     = "T")

message("Building variants ...")
variants <- read_variants(vcf, region)

message("Building gene annotations (downloads the CHM13 GTF on first run) ...")
annotations <- read_annotations(genome = "chm13", region = region)

cache <- list(
  meth_basic  = meth_basic,
  meth_6ma    = meth_6ma,
  meth_hp     = meth_hp,
  meth_snv    = meth_snv,
  variants    = variants,
  annotations = annotations,
  region      = region
)

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
saveRDS(cache, "inst/extdata/vignette_cache.rds", compress = "xz")

message("Cache saved to inst/extdata/vignette_cache.rds")
message("File size: ",
        format(file.info("inst/extdata/vignette_cache.rds")$size,
               big.mark = ","),
        " bytes")
