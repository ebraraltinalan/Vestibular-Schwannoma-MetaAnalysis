#volcano plot
library(ggrepel)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(patchwork)

# dataframe: meta_analysis_results_all.xlsx
meta_analysis_results_all <- meta_analysis_results_all


lfc_cutoff1 <- 1
pval_cutoff <- 0.05

# Mark the significant / non-significant colours
all_genes_volcano <- meta_analysis_results_all %>%
  mutate(significance1 = ifelse(abs(meta_LFc) > lfc_cutoff1 & meta_pval < pval_cutoff, "Significant", "NS"))

x_min <- min(all_genes_volcano$meta_LFc, na.rm=TRUE) - 0.5
x_max <- max(all_genes_volcano$meta_LFc, na.rm=TRUE) + 0.5


# Only significant 50 genes
top_genes <- all_genes_volcano %>%
  filter(significance1 == "Significant") %>%
  arrange(meta_pval) %>%
  slice(1:50)  # en düşük p-value'lu 50 gen

# Volcano plot
p2 <- ggplot(all_genes_volcano, aes(x = meta_LFc, y = -log10(meta_pval), color = significance1)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Significant" = "red", "NS" = "grey")) +
  geom_vline(xintercept = c(-lfc_cutoff1, lfc_cutoff1), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed") +
  geom_text_repel(
    data = top_genes,
    aes(label = `Gene Symbol`),
    size = 3
  ) +
  labs(
    title = "Volcano Plot |metaLFC| > 1",
    x = "Meta-analysed log2 fold change (metaLFC)",
    y = "-log10(meta p-value)"
  )

# Save as PNG
png("Volcano_plot_LFC1_top50.png", width = 1200, height = 800)
print(p2)
dev.off()



#or classfy by up/down regulated

library(ggplot2)
library(ggrepel)
library(dplyr)

lfc_cutoff <- 1
pval_cutoff <- 0.05

# Volcano plot 
all_genes_volcano2 <- meta_analysis_results_all %>%
  mutate(
    regulation = case_when(
      meta_LFc > lfc_cutoff & meta_pval < pval_cutoff ~ "Up-regulated",
      meta_LFc < -lfc_cutoff & meta_pval < pval_cutoff ~ "Down-regulated",
      TRUE ~ "NS"
    )
  )

x_min2 <- min(all_genes_volcano2$meta_LFc, na.rm=TRUE) - 0.5
x_max2 <- max(all_genes_volcano2$meta_LFc, na.rm=TRUE) + 0.5

# Only signifcant 20 genes
top_genes2 <- all_genes_volcano2 %>%
  filter(regulation != "NS") %>%
  arrange(meta_pval) %>%
  slice(1:20)

# Volcano plot
p5 <- ggplot(all_genes_volcano2, aes(x = meta_LFc, y = -log10(meta_pval), color = regulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "NS" = "grey")) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed") +
  geom_text_repel(
    data = top_genes2,
    aes(label = `Gene Symbol`),
    size = 6,
    fontface = "bold" 
  ) +
  labs(
    title = "Volcano Plot |metaLFC| > 1",
    x = "Meta-analysed log2 fold change (metaLFC)",
    y = "-log10(meta p-value)"
  ) +

  scale_x_continuous(
    limits = c(x_min2, x_max2),  # mevcut aralıkları koru
    breaks = c(-4, -3, -2,-1, 0, 1, 2, 3, 4),      # sadece -1, 0 ve 1 değerlerini göster
    labels = c("-4","-3", "-2", "-1", "0", "1", "2", "3", "4") # istersen metinleri özelleştirebilirsin
  ) +
  theme_minimal() +
theme(
  axis.text.x = element_text(face = "bold"),
  axis.text.y = element_text(face = "bold")
)

# Save as PNG
png("Volcano_plot_LFC1_top20.png", width = 1200, height = 800,bg = "white")
print(p5)
dev.off()



