#volcano_plot for withcovs

library(ggrepel)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(patchwork)

# dataframe: meta_analysis_results_withcov.xlsx
meta_analysis_results_withcov <- meta_analysis_results_withcov

# Cutoff values

lfc_cutoff2 <- 1
pval_cutoff <- 0.05

# Mark the significant / non-significant colours
  all_genes_withcov_volcano <- meta_analysis_results_withcov %>%
  mutate(significance = ifelse(abs(meta_LFc) > lfc_cutoff2 & meta_pval < pval_cutoff, "Significant", "NS"))

x_min <- min(all_genes_withcov_volcano$meta_LFc, na.rm=TRUE) - 0.5
x_max <- max(all_genes_withcov_volcano$meta_LFc, na.rm=TRUE) + 0.5


# Only significant 50 genes
top_genes_withcov <- all_genes_withcov_volcano %>%
  filter(significance2 == "Significant") %>%
  arrange(meta_pval) %>%
  slice(1:50)  

# Volcano plot
p3 <- ggplot(all_genes_withcov_volcano, aes(x = meta_LFc, y = -log10(meta_pval), color = significance2)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Significant" = "red", "NS" = "grey")) +
  geom_vline(xintercept = c(-lfc_cutoff2, lfc_cutoff2), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed") +
  geom_text_repel(
    data = top_genes_withcov,
    aes(label = `Gene Symbol`),
    size = 3
  ) +
  labs(title = "Volcano Plot |Log2FC| > 1", x = "Log2 Fold Change", y = "-log10(p-value)") +
  xlim(x_min, x_max) +
  theme_minimal()

# Save as PNG
png("Volcano_plot_LFC1_top50_withcov.png", width = 1200, height = 800)
print(p3)
dev.off()

colnames(meta_analysis_results_withcov)



#or classfy by up/down regulated

library(ggplot2)
library(ggrepel)
library(dplyr)



lfc_cutoff <- 1
pval_cutoff <- 0.05

all_genes_withcov_volcano2 <- meta_analysis_results_withcov %>%
  mutate(
    regulation = case_when(
      meta_LFc > lfc_cutoff & meta_pval < pval_cutoff ~ "Up-regulated",
      meta_LFc < -lfc_cutoff & meta_pval < pval_cutoff ~ "Down-regulated",
      TRUE ~ "NS"
    )
  )


x_min2 <- min(all_genes_withcov_volcano2$meta_LFc, na.rm=TRUE) - 0.5
x_max2 <- max(all_genes_withcov_volcano2$meta_LFc, na.rm=TRUE) + 0.5

# Only signifcant 20 genes
top_genes3 <- all_genes_withcov_volcano2 %>%
  filter(regulation != "NS") %>%
  arrange(meta_pval) %>%
  slice(1:20)

# Volcano plot
p4 <- ggplot(all_genes_withcov_volcano2, aes(x = meta_LFc, y = -log10(meta_pval), color = regulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "NS" = "grey")) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed") +
  geom_text_repel(
    data = top_genes3,
    aes(label = `Gene Symbol`),
    size = 6,
    fontface = "bold" 
  ) +
  labs(title = "Volcano Plot |Log2FC| > 1", x = "Log2 Fold Change", y = "-log10(p-value)") +
  scale_x_continuous(
    limits = c(x_min2, x_max2),  
    breaks = c(-4, -3, -2,-1, 0, 1, 2, 3, 4),      
    labels = c("-4","-3", "-2", "-1", "0", "1", "2", "3", "4") 
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

# Save as PNG
png("Volcano_plot_LFC1_top20_withcov_new.png", width = 1200, height = 800)
print(p4)
dev.off()
ggsave("Volcano_plot_LFC1_top20_withcov.png", plot = p4, width = 12, height = 8, dpi = 150)


# Count up- and down-regulated genes
gene_counts <- all_genes_withcov_volcano2 %>%
  filter(regulation != "NS") %>%   
  summarise(
    Up = sum(regulation == "Up-regulated"),
    Down = sum(regulation == "Down-regulated"),
    Total = n()
  )

gene_counts 
#  Up  Down Total
#  1   745   479  1224

