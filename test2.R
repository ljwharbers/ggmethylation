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
REGION     <- "chr8:127598491-127604951" #hg38
REGION     <- "chr8:128861685-128872405" #CHM13
REGION     <- "chr8:128861685-128865405" #CHM13


annotation = read_annotations(genome = "chm13", region = REGION)


meth_hp <- read_methylation(BAM_PTCL8,
                            REGION,
                            mod_code    = "m",
                            group_tag   = "HP",
                            drop_na_group = FALSE)

plt = plot_methylation(meth_hp,
                 annotations = annotation,
                 show_cigar = TRUE,
                 show_supplementary = TRUE,
                 line_width = 0.2,
                 group_colours = c("1" = "#95babc", "2" = "#efbb76"),
                 panel_heights = c(1, 0.2, 0.04))
ggsave("test.png", plot = plt, width = 8, height = 6)
