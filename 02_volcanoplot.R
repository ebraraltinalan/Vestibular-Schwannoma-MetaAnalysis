#volcano plot

library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)


base_dir <- "analysis"

volcano_file <- file.path(
  base_dir,
  "rechecked_without_cov",
  "meta_analysis_results_all_with_padj.xlsx"
)

output_dir <- file.path(
  base_dir,
  "rechecked_without_cov",
  "volcano_padj"
)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


meta_analysis_results_all <- read_excel(volcano_file)

colnames(meta_analysis_results_all) <- trimws(colnames(meta_analysis_results_all))

meta_analysis_results_all <- meta_analysis_results_all %>%
  mutate(
    meta_LFc = as.numeric(meta_LFc),
    meta_pval = as.numeric(meta_pval),
    padj = as.numeric(padj)
  )


lfc_cutoff <- 1
padj_cutoff <- 0.05
top_n <- 20


# Volcano data
all_genes_volcano <- meta_analysis_results_all %>%
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

# Select labelled genes and axis limits


top_genes <- all_genes_volcano %>%
  filter(regulation != "NS") %>%
  arrange(padj) %>%
  slice_head(n = top_n)

top_genes %>%
  select(`Gene Symbol`, meta_LFc, meta_pval, padj, neg_log10_padj, regulation)


x_min <- min(all_genes_volcano$meta_LFc, na.rm = TRUE) - 0.5
x_max <- max(all_genes_volcano$meta_LFc, na.rm = TRUE) + 0.5
y_max <- max(all_genes_volcano$neg_log10_padj, na.rm = TRUE) + 1


#  Plot


p_volcano_nocov <- ggplot(
  all_genes_volcano,
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
    data = top_genes,
    aes(
      x = meta_LFc,
      y = neg_log10_padj,
      label = `Gene Symbol`
    ),
    inherit.aes = FALSE,
    color = ifelse(top_genes$regulation == "Up-regulated", "red", "blue"),
    size = 2.9,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.15,
    point.padding = 0.1,
    segment.color = NA,
    show.legend = FALSE
  ) +
  labs(
    title = "Volcano Plot Without Covariates Model",
    x = "meta log2 fold change (metaLFC)",
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
    legend.title = element_blank()
  )

ggsave(
  filename = file.path(output_dir, "Volcano_plot_without_cov_metaLFC1_padj_top20.png"),
  plot = p_volcano_nocov,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(output_dir, "Volcano_plot_without_cov_metaLFC1_padj_top20.pdf"),
  plot = p_volcano_nocov,
  width = 12,
  height = 8,
  bg = "white"
)

#Save labelled genes
write.csv(
  top_genes,
  file.path(output_dir, "labelled_top20_genes_without_cov_padj2.csv"),
  row.names = FALSE
)
p_volcano_nocov
