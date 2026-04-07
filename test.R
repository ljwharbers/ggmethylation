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

BAM_PTCL8  <- "/staging/leuven/stg_00096/home/averham/LR_SOMATIC_T2T/PTCL8_PB/bamfiles/PTCL8_PB_tumor.bam"
BAM_AITL4  <- "/staging/leuven/stg_00096/home/averham/LR_SOMATIC_T2T/AITL4/bamfiles/AITL4_tumor.bam"
VCF_AITL4  <- "/staging/leuven/stg_00096/home/averham/LR_SOMATIC_T2T/AITL4/variants/clairsto/somatic.vcf.gz"
GTF        <- "/staging/leuven/stg_00096/references/chm13_v2.0_maskedY.rCRS/annotation/chm13v2.0_RefSeq_Liftoff_v5.1.gtf"

REGION     <- "chr8:127730434-127750951"


# =============================================================================
# 1. Basic single-sample, ungrouped
# =============================================================================

meth <- read_methylation(BAM_PTCL8, REGION, mod_code = "m")

print(meth)
summary(meth)

plot_methylation(meth)


# =============================================================================
# 2. Grouped by haplotype tag (HP)
# =============================================================================

meth_hp <- read_methylation(BAM_PTCL8, REGION,
                             mod_code    = "m",
                             group_tag   = "HP",
                             drop_na_group = TRUE)

print(meth_hp)

plot_methylation(meth_hp)


# =============================================================================
# 3. Grouped by SNV genotype
# =============================================================================

meth_snv <- read_methylation(BAM_PTCL8, REGION,
                              mod_code     = "m",
                              snv_position = 127740000L,
                              ref_base     = "A",
                              alt_base     = "T")

plot_methylation(meth_snv)


# =============================================================================
# 4. Read filters — MAPQ, strand, length, downsampling
# =============================================================================

meth_filtered <- read_methylation(BAM_PTCL8, REGION,
                                   mod_code        = "m",
                                   min_mapq        = 20L,
                                   strand_filter   = "+",
                                   min_read_length = 5000L,
                                   max_reads       = 50L)

plot_methylation(meth_filtered)


# =============================================================================
# 5. sort_by options
# =============================================================================

plot_methylation(meth_hp, sort_by = "group")
plot_methylation(meth_hp, sort_by = "start")
plot_methylation(meth_hp, sort_by = "mean_mod_prob")


# =============================================================================
# 6. Colour and aesthetics customisation
# =============================================================================

plot_methylation(meth,
                 colour_low  = "#FFFFFF",
                 colour_high = "#1A237E",
                 dot_size    = 2)

# Colour reads by strand instead of a single colour
plot_methylation(meth,
                 colour_strand  = TRUE,
                 strand_colours = c("+" = "#2196F3", "-" = "#FF5722"))

# Custom group colours
plot_methylation(meth_hp,
                 group_colours = c("1" = "#388E3C", "2" = "#7B1FA2"))

# Adjust panel height ratio (reads : smooth)
plot_methylation(meth_hp, panel_heights = c(4, 1))


# =============================================================================
# 7. Gene annotations — from TxDb and from GTF (with session cache)
# =============================================================================

# From GTF — parsed once, cached for the session
annot <- read_annotations(gtf = GTF, region = REGION)

print(annot)

# Second call reuses cache (no re-parsing message)
annot2 <- read_annotations(gtf = GTF, region = "chr8:127700000-127800000")

# Clear cache to free memory or force re-parse
clear_annotation_cache()

# From a pre-built TxDb (e.g., Bioconductor annotation package)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# annot_txdb <- read_annotations(txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#                                region = REGION)

plot_methylation(meth_hp, annotations = annot)


# =============================================================================
# 8. Variant overlay (VCF)
# =============================================================================

vars <- read_variants(VCF_AITL4, REGION)

print(vars)

meth_aitl4 <- read_methylation(BAM_AITL4, REGION,
                                mod_code  = "m",
                                group_tag = "HP")

plot_methylation(meth_aitl4, variants = vars)

# Combined: annotations + variants
plot_methylation(meth_aitl4,
                 annotations = annot,
                 variants    = vars)


# =============================================================================
# 9. Multi-sample comparison
# =============================================================================

merged <- merge_methylation(PTCL8 = meth_hp,
                             AITL4 = meth_aitl4)

# Equivalent using .list
merged <- merge_methylation(.list = list(PTCL8 = meth_hp, AITL4 = meth_aitl4))

print(merged)
summary(merged)

plot_methylation(merged)

plot_methylation(merged, annotations = annot, variants = vars)


# =============================================================================
# 10. Write methylation data to disk
# =============================================================================

write_methylation(meth, prefix = "/tmp/ptcl8_meth", format = "tsv")
write_methylation(meth, prefix = "/tmp/ptcl8_meth", format = "bed", gzip = TRUE)

# Pipe-friendly (returns data invisibly)
meth |>
  write_methylation(prefix = "/tmp/ptcl8_pipe", format = "tsv") |>
  plot_methylation()


# =============================================================================
# 11. Manipulating patchwork output
# =============================================================================

p <- plot_methylation(meth_hp, annotations = annot)

# Access individual panels
reads_panel  <- p[[1]]
smooth_panel <- p[[2]]
gene_panel   <- p[[3]]

# Modify a panel and reassemble
reads_panel + labs(y = "Read lane") + theme_classic()

# Add a global title
p + plot_annotation(
  title    = REGION,
  subtitle = "CpG methylation — PTCL8 tumour (haplotype-grouped)"
)

# Apply a theme across all panels
p & theme(text = element_text(size = 14))
p & theme_bw()

# Mark a position of interest
reads_panel + geom_vline(xintercept = 127740000, colour = "blue",
                          linetype = "dashed", linewidth = 0.5)
