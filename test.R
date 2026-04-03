install.packages("/staging/leuven/stg_00096/home/lharbers/repositories/ggmethylation/", repos = NULL, type = "source")
require(ggmethylation)
require(ggplot2)
require(patchwork)

region = 'chr8:127730434-127750951'

annotations = read_annotations(gtf = "/staging/leuven/stg_00096/references/chm13_v2.0_maskedY.rCRS/annotation/chm13v2.0_RefSeq_Liftoff_v5.1.gtf",
                               region = region)

bamfile <- '/staging/leuven/stg_00096/home/averham/LR_SOMATIC_T2T/PTCL8_PB/bamfiles/PTCL8_PB_tumor.bam'
bamfile2 <- '/staging/leuven/stg_00096/home/averham/LR_SOMATIC_T2T/AITL4/bamfiles/AITL4_tumor.bam'
vcf_file = "/staging/leuven/stg_00096/home/averham/LR_SOMATIC_T2T/AITL4/variants/clairsto/somatic.vcf.gz"



meth_data <- read_methylation(bamfile, region,
                              mod_code = "m",
                              group_tag = "HP")

meth_data2 <- read_methylation(bamfile2, region,
                               mod_code = "m",
                               group_tag = "HP")

vcf_data = read_variants(vcf_file, region)

merged_data = merge_methylation(.list = list(PTCL8 = meth_data, AITL4 = meth_data2))

p <- plot_methylation(merged_data,
                      dot_size = 0.8,
                      variants = vcf_data,
                      annotations = annotations)

# --- Saving ---
ggsave("test_plot.png", p, width = 14, height = 6, dpi = 150)

# --- Theming ---
# Override theme elements on a patchwork or ggplot object
p & theme(text = element_text(size = 14))
p & theme_bw()

# --- Titles / labels ---
# patchwork::plot_annotation() adds a title/subtitle to the composite
p + plot_annotation(
  title    = "chr8:128830975-130830975",
  subtitle = "CpG methylation — PTCL8 tumour"
)

# --- Accessing individual panels ---
# patchwork stores panels as a list; index with [[
top_panel    <- p[[1]]   # reads panel
bottom_panel <- p[[2]]   # smoothed line panel

# Modify one panel and reassemble
top_modified <- top_panel + labs(y = "Read lane") + theme_classic()
top_modified / bottom_panel   # reassemble with patchwork /

# --- Adding layers to a panel ---
# e.g. mark a position of interest with a vertical line
top_panel + geom_vline(xintercept = 128865000, colour = "blue", linetype = "dashed")

# --- Adjusting panel size ratios ---
plot_methylation(meth_data, panel_heights = c(4, 1))  # taller reads panel

# --- Inspecting the ggplot object ---
# List all layers in the top panel
top_panel$layers

# Extract underlying data used for the smoothed line
bottom_panel$layers[[1]]$data
