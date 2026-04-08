# =============================================================================
# ggmethylation — interactive showcase script
#
# Demonstrates all major use cases of the package.
# Requires access to the staging BAM/VCF files on the HPC cluster.
#
# Run sections interactively; do not source() the whole file at once.
# =============================================================================

devtools::load_all()  # or: install.packages(".", repos = NULL, type = "source")

library(ggplot2)
library(patchwork)

# -- Shared paths -------------------------------------------------------------

BAM_PTCL8  <- "/staging/leuven/stg_00096/home/averham/LR_SOMATIC_hg38/PTCL8_PB/bamfiles/PTCL8_PB_tumor.bam"
REGION     <- "chr8:127598491-127604951"

# =============================================================================
# 1. Basic single-sample, ungrouped
# =============================================================================

# =============================================================================
# 2. Grouped by haplotype tag (HP)
# =============================================================================

meth_hp <- read_methylation(BAM_PTCL8, REGION,
                             mod_code    = "m",
                             group_tag   = "HP",
                             drop_na_group = TRUE)

plot_methylation(meth_hp, show_cigar = TRUE)
