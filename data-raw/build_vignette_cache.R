# Build vignette cache
# Run this script manually on HPC to regenerate inst/extdata/vignette_cache.rds
# Package version: see DESCRIPTION
# BAM path: /staging/leuven/stg_00096/home/lharbers/repositories/ggmethylation/data/PTCL8_PB_tumor_chr21_22_subset.bam
#
# Usage (from the package root):
#   Rscript data-raw/build_vignette_cache.R

devtools::load_all()

bam    <- "/staging/leuven/stg_00096/home/lharbers/repositories/ggmethylation/data/PTCL8_PB_tumor_chr21_22_subset.bam"
region <- "chr21:34500000-34510000"

message("Building meth_basic ...")
meth_basic <- read_methylation(bam, region)

message("Building meth_hp ...")
meth_hp <- read_methylation(bam, region, group_tag = "HP")

# For the SNV example, identify a heterozygous SNV in the region and supply its
# coordinate, reference allele, and alternate allele below.  If no clean SNV is
# available in the region, meth_hp is used as a stand-in so the vignette still
# demonstrates the code path.
#
# To find candidate SNV positions:
#   head(meth_basic$reads)        # inspect the data frame columns
#   table(meth_basic$reads$pos)   # look for positions with bi-allelic coverage
#
# Replace the placeholder below once a suitable position is identified:
#   meth_snv <- read_methylation(
#     bam, region,
#     snv_position = <INTEGER_POSITION>,
#     ref_base     = "C",
#     alt_base     = "T"
#   )

message("Building meth_snv (using HP grouping as placeholder) ...")
meth_snv <- meth_hp   # placeholder; replace with actual SNV call when position is known

cache <- list(
  meth_basic = meth_basic,
  meth_hp    = meth_hp,
  meth_snv   = meth_snv
)

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
saveRDS(cache, "inst/extdata/vignette_cache.rds", compress = "xz")

message("Cache saved to inst/extdata/vignette_cache.rds")
message("File size: ",
        format(file.info("inst/extdata/vignette_cache.rds")$size,
               big.mark = ","),
        " bytes")
