require(ggmethylation)

bamfile <- '/staging/leuven/stg_00096/home/averham/LR_SOMATIC_T2T/PTCL8_PB/bamfiles/PTCL8_PB_tumor.bam'
region  <- 'chr8:128862888-128870405'

meth_data <- read_methylation(bamfile, region,
                              mod_code = 'm',
                              group_tag = 'HP')
p <- plot_methylation(meth_data)

# --- Saving ---
ggsave("test_plot.png", p, width = 14, height = 6, dpi = 150)

# --- Theming ---
# Override theme elements on a patchwork or ggplot object
p & theme(text = element_text(size = 14))
p & theme_bw()

# --- Titles / labels ---
# patchwork::plot_annotation() adds a title/subtitle to the composite
p + plot_annotation(
  title    = "chr8:128862888-128870405",
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
