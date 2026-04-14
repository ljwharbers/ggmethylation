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

BAM_PTCL8  <- "/staging/leuven/stg_00096/home/averham/LR_SOMATIC_T2T/BL1/bamfiles/BL1_tumor.bam"
VCF_SVs = "/staging/leuven/stg_00096/home/averham/LR_SOMATIC_T2T/BL1/variants/severus/somatic_SVs/severus_somatic.vcf.gz"
VCF_SNVs = "/staging/leuven/stg_00096/home/averham/LR_SOMATIC_T2T/BL1/variants/clairsto/somatic.vcf.gz"
# REGION     <- "chr8:127598491-127604951" #hg38
REGION     <- "chr8:128861685-128872405" #CHM13
# REGION     <- "chr8:128801685-128885405" #CHM13
# REGION = "chr12:6542873-6650727"
REGION2 = "chr16:3685304-3885304"

annotation = read_annotations(genome = "chm13", region = REGION)
vars <- read_variants(VCF_SVs, REGION)

meth_hp <- read_methylation(BAM_PTCL8,
                            REGION,
                            mod_code    = "m",
                            group_tag   = "HP",
                            drop_na_group = FALSE)

plt = plot_methylation(meth_hp,
                 annotations = annotation,
                 variants = vars)
ggsave("test.png", plot = plt, width = 8, height = 6)

########

annotation2 = read_annotations(genome = "chm13", region = REGION2)
vars2 <- read_variants(VCF_SNVs, REGION2)

meth_hp2 <- read_methylation(BAM_PTCL8,
                            REGION2,
                            mod_code    = "m",
                            group_tag   = "HP",
                            drop_na_group = FALSE)

plt2 = plot_methylation(meth_hp2,
                 annotations = annotation2,
                 variants = vars2)
ggsave("test2.png", plot = plt2, width = 8, height = 6)

