#enrichment analysis GO
#withcovs
library(dplyr)
library(tidyverse)
library(writexl)
library(readxl)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)


#added gene name before (gene_info_2)
library(readxl)
meta_analysis_results_withcov = read_xlsx("/Users/ebraraltinalan/Desktop/GEO/final_withcov_geo/files/meta_analysis_results_withcov.xlsx")


meta_analysis_genes_names_withcov <- meta_analysis_results_withcov %>%
  left_join(gene_info2, by = "Gene Symbol")

write_xlsx(meta_analysis_genes_names_withcov, "meta_analysis_genes_names_withcov.xlsx")

# Convert  Gene Symbol to  entrez ID

## —— Gene Symbol → ENTREZ ID  (mapIds) ——
entrez_ids2 <- mapIds(x = org.Hs.eg.db,         #org.Hs.eg.db: Bioconductor package human gene anotation database
                     keys = meta_analysis_genes_names_withcov$`Gene Symbol`,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")

data_mapped_withcov <- meta_analysis_genes_names_withcov %>%
  mutate(ENTREZID = as.character(entrez_ids2)) %>%
  filter(!is.na(ENTREZID))
#21116 entries

data_mapped_clean_withcov <- data_mapped_withcov %>%
  filter(!is.na(ENTREZID)) %>%   # 1) exclude NA
  distinct(ENTREZID, .keep_all = TRUE)


sum(data_mapped_clean_withcov$meta_pval < 0.05, na.rm = TRUE)   
sum(data_mapped_clean_withcov$padj < 0.05, na.rm = TRUE)


## 21,113 entries

## —— 3) Universe and significant list ——

universe_entrez2 <- unique(data_mapped_clean_withcov$ENTREZID)

data_mapped_clean_withcov <- data_mapped_clean_withcov %>%
  mutate(padj = p.adjust(meta_pval, method = "BH")) 

sig_genes_entrez2 <- data_mapped_clean_withcov %>%
  filter(padj < 0.05, abs(meta_LFc) > 0.5) %>%       
  pull(ENTREZID) %>%
  unique()

# Security control
if (length(sig_genes_entrez2) == 0) stop("Significant gen bulunamadı (padj < 0.05). Eşiği gevşetmeyi deneyin.")

## ——  enrichGO (BP) + setReadable ——
ego_bp2 <- enrichGO(gene          = sig_genes_entrez2,
                   universe      = universe_entrez2,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = FALSE)

##  setReadable —— ENTREZ → SYMBOL 
ego_bp2_read <- setReadable(ego_bp2, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")


## save
ego_bp2_df <- as.data.frame(ego_bp2_read)
ego_bp2_df$GeneFirst <- sapply(strsplit(ego_bp2_df$geneID, "/"), `[`, 1)
write.csv(ego_bp2_df, "GO_ORA_BP_results_readable2_withcov.csv", row.names = FALSE)

write_xlsx(ego_bp2_df,"GO_table1_withcov.xlsx")
# logFC filtresi
filtered_withcov <- data_mapped_clean_withcov %>%
  filter(abs(meta_LFc) > 0.5)


# 
filtered_withcov <- filtered_withcov %>%
  arrange(meta_pval)

# top10
top10_withcov <- head(filtered_withcov, 10)

library(enrichplot)
library(ggplot2)
library(stringr)

# barplot
p2 <- barplot(ego_bp2_read,
             showCategory = 10,
             x = "GeneRatio",            # Count yerine oransal gösterim
             title = "GO Biological Process (Top 10)")

# p2
p2 <- p2 +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    plot.title  = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_blank()
  ) +
  labs(x = "Gene Ratio")

print(p2)

# save
ggsave("GO_BP_barplot_top10_withcov.png", plot = p2, width = 7, height = 5, dpi = 300)
ggsave("GO_BP_barplot_top10_withcov.pdf", plot = p2, width = 7, height = 5)


#another option

library(dplyr)
library(ggplot2)
library(stringr)

# enrichment results dataframe
df_5 <- as.data.frame(ego_bp2_read) %>%
  arrange(p.adjust) %>%
  slice(1:15) %>%                             
  mutate(Description_wrapped = str_wrap(Description, width = 45),
         negLog10FDR = -log10(p.adjust))

# barplot + count points
g5 <- ggplot(df_5, aes(x = reorder(Description_wrapped, negLog10FDR),
                      y = negLog10FDR)) +
  geom_col(fill = "steelblue") +
  geom_point(aes(y = negLog10FDR, size = Count), color = "red", shape = 21, stroke = 1.2) +
  coord_flip() +
  labs(title = "GO Biological Process Enrichment (Top 15)",
       x = NULL,
       y = expression(-log[10]~"(adjusted p-value)"),
       size = "Gene Count") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 13, face = "bold"))

print(g5)

ggsave("GO_BP_barplot_top15_custom_withcov.png", g5, width = 8, height = 6, dpi = 300)

#change the colors

df_5 <- df_5 %>%
  mutate(Count_cat = factor(Count)) %>%
  filter(!is.na(negLog10FDR), !is.na(Count))   

g5 <- ggplot(df_5, aes(x = reorder(Description_wrapped, negLog10FDR),
                      y = negLog10FDR)) +
  geom_col(fill = "steelblue") +
  geom_point(aes(color = Count_cat), size = 3) +
  scale_color_viridis_d(option = "magma") +   
  coord_flip() +
  labs(title = "GO Biological Process Enrichment (Top 15)",
       x = NULL,
       y = expression(-log[10]~"(adjusted p-value)"),
       color = "Gene Count") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title.position = "plot",            
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.position = "right",        
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  )
  



print(g5)
ggsave("GO_BP_barplot_top15_colored_withcov.png", g5, width = 8, height = 6, dpi = 300)


#calculation
gene_counts <- filtered_withcov %>%
  mutate(
    regulation = case_when(
      meta_LFc > 0.5 ~ "Up-regulated",
      meta_LFc < -0.5 ~ "Down-regulated",
      TRUE ~ "NS"
    )
  ) %>%
  filter(regulation != "NS") %>%
  summarise(
    Up = sum(regulation == "Up-regulated"),
    Down = sum(regulation == "Down-regulated"),
    Total = n()
  )

gene_counts



#for top20

# enrichment results dataframe
df_6 <- as.data.frame(ego_bp2_read) %>%
  arrange(p.adjust) %>%
  slice(1:20) %>%                             
  mutate(Description_wrapped = str_wrap(Description, width = 45),
         negLog10FDR = -log10(p.adjust))

# barplot + count points
g6 <- ggplot(df_6, aes(x = reorder(Description_wrapped, negLog10FDR),
                       y = negLog10FDR)) +
  geom_col(fill = "steelblue") +
  geom_point(aes(y = negLog10FDR, size = Count), color = "red", shape = 21, stroke = 1.2) +
  coord_flip() +
  labs(title = "GO Biological Process Enrichment (Top 20)",
       x = NULL,
       y = expression(-log[10]~"(adjusted p-value)"),
       size = "Gene Count") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 13, face = "bold"))

print(g6)

# save
ggsave("GO_BP_barplot_top20_custom_withcov.png", g6, width = 8, height = 6, dpi = 300)

#change the colors

df_6 <- df_6 %>%
  mutate(Count_cat = factor(Count)) %>%
  filter(!is.na(negLog10FDR), !is.na(Count))   

g6 <- ggplot(df_6, aes(x = reorder(Description_wrapped, negLog10FDR),
                       y = negLog10FDR)) +
  geom_col(fill = "steelblue") +
  geom_point(aes(color = Count_cat), size = 3) +
  scale_color_viridis_d(option = "magma") +   
  coord_flip() +
  labs(title = "GO Biological Process Enrichment (Top 20)",
       x = NULL,
       y = expression(-log[10]~"(adjusted p-value)"),
       color = "Gene Count") +
  theme(
    plot.title.position = "plot",              
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.position = "right",        
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  )


print(g6)

