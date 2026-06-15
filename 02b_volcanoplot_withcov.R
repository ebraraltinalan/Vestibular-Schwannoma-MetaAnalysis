#volcano_plot for withcovs

library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)

base_dir <- "withcov_geo/files"

volcano_file <- file.path(
  base_dir,
  "meta_analysis_results_withcov.xlsx"
)

output_dir <- file.path(
  base_dir,
  "volcano_padj_withcov"
)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

meta_analysis_results_withcov <- read_excel(volcano_file)

colnames(meta_analysis_results_withcov) <- trimws(colnames(meta_analysis_results_withcov))

meta_analysis_results_withcov <- meta_analysis_results_withcov %>%
  mutate(
    meta_LFc = as.numeric(meta_LFc),
    meta_pval = as.numeric(meta_pval)
  )

# If padj column does not exist, calculate it from meta_pval
if (!"padj" %in% colnames(meta_analysis_results_withcov)) {
  meta_analysis_results_withcov <- meta_analysis_results_withcov %>%
    mutate(
      padj = p.adjust(meta_pval, method = "BH")
    )
} else {
  meta_analysis_results_withcov <- meta_analysis_results_withcov %>%
    mutate(
      padj = as.numeric(padj)
    )
}


write_xlsx(
  meta_analysis_results_withcov,
  file.path(output_dir, "meta_analysis_results_withcov_with_padj.xlsx")
)


lfc_cutoff <- 1
padj_cutoff <- 0.05
top_n <- 20


#Volcano data
all_genes_volcano_withcov <- meta_analysis_results_withcov %>%
  filter(!is.na(meta_LFc), !is.na(padj)) %>%
  mutate(
    regulation = case_when(
      meta_LFc > lfc_cutoff & padj < padj_cutoff ~ "Up-regulated",
      meta_LFc < -lfc_cutoff & padj < padj_cutoff ~ "Down-regulated",
      TRUE ~ "NS"
    ),
    plot_padj = pmax(padj, 1e-50),
    neg_log10_padj = -log10(plot_padj)
  )

# Check gene counts
table(all_genes_volcano_withcov$regulation)

all_genes_volcano_withcov %>%
  filter(regulation != "NS") %>%
  count(regulation)


# labelled genes
top_genes_withcov <- all_genes_volcano_withcov %>%
  filter(regulation != "NS") %>%
  arrange(padj) %>%
  slice_head(n = top_n)

top_genes_withcov %>%
  select(`Gene Symbol`, meta_LFc, meta_pval, padj, neg_log10_padj, regulation)


#Plot
x_min <- min(all_genes_volcano_withcov$meta_LFc, na.rm = TRUE) - 0.5
x_max <- max(all_genes_volcano_withcov$meta_LFc, na.rm = TRUE) + 0.5
y_max <- max(all_genes_volcano_withcov$neg_log10_padj, na.rm = TRUE) + 1

p_volcano_withcov <- ggplot(
  all_genes_volcano_withcov,
  aes(x = meta_LFc, y = neg_log10_padj, color = regulation)
) +
  geom_point(alpha = 0.6, size = 1.6, shape = 16) +
  scale_color_manual(
    values = c(
      "Down-regulated" = "blue",
      "NS" = "grey",
      "Up-regulated" = "red"
    ),
    breaks = c("Down-regulated", "NS", "Up-regulated")
  ) +
  geom_vline(
    xintercept = c(-lfc_cutoff, lfc_cutoff),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype = "dashed"
  ) +
  geom_text_repel(
    data = top_genes_withcov,
    aes(
      x = meta_LFc,
      y = neg_log10_padj,
      label = `Gene Symbol`
    ),
    inherit.aes = FALSE,
    color = ifelse(top_genes_withcov$regulation == "Up-regulated", "red", "blue"),
    size = 2.9,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.15,
    point.padding = 0.1,
    segment.color = "grey40",
    segment.size = 0.25,
    segment.alpha = 0.7,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  labs(
    title = "Volcano Plot With Covariates Model",
    x = "Meta-analysed log2 fold change (metaLFC)",
    y = "-log10(adjusted p-value)",
    color = NULL
  ) +
  coord_cartesian(
    xlim = c(x_min, x_max),
    ylim = c(0, y_max)
  ) +
  scale_x_continuous(
    breaks = seq(floor(x_min), ceiling(x_max), by = 1)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = 16,
        size = 4,
        alpha = 1
      )
    )
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )



ggsave(
  filename = file.path(output_dir, "Volcano_plot_with_cov_metaLFC1_padj_top20.png"),
  plot = p_volcano_withcov,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(output_dir, "Volcano_plot_with_cov_metaLFC1_padj_top20.pdf"),
  plot = p_volcano_withcov,
  width = 12,
  height = 8,
  bg = "white"
)


# Save labelled genes
write.csv(
  top_genes_withcov,
  file.path(output_dir, "labelled_top20_genes_with_cov_padj.csv"),
  row.names = FALSE
)

p_volcano_withcov


