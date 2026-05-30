#enrichment analysis GO

#add gene name
gene_info2 <- Homo_sapiens_gene_info %>%      
  select(V3, V9) %>%
  rename(`Gene Symbol` = V3, GeneName = V9)

meta_analysis_all_genes_names <- meta_analysis_results_all %>%
  left_join(gene_info2, by = "Gene Symbol")

write_xlsx(meta_analysis_all_genes_names, "meta_analysis_all_genes_names.xlsx")


meta_analysis_sig_genes_names <- meta_analysis_results_significant_genes %>%
  left_join(gene_info2, by = "Gene Symbol")
write_xlsx(meta_analysis_sig_genes_names, "meta_analysis_sig_genes_names.xlsx")


# Convert  Gene Symbol to  entrez ID
data <- read_excel("/Users/ebraraltinalan/Desktop/meta_analysis_all_genes_names.xlsx")


library(readxl)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot) 

## —— Gene Symbol → ENTREZ ID  (mapIds) ——
entrez_ids <- mapIds(x = org.Hs.eg.db,         #org.Hs.eg.db: Bioconductor package human gene anotation database
                     keys = data$`Gene Symbol`,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")      

dat_mapped <- data %>%
  mutate(ENTREZID = as.character(entrez_ids)) %>%
  filter(!is.na(ENTREZID))

data_mapped_clean <- dat_mapped %>%
  filter(!is.na(ENTREZID)) %>%   # 1) exclude NA
  distinct(ENTREZID, .keep_all = TRUE)  


data_mapped_clean %>%
  filter(padj < 0.05, abs(meta_LFc) > 0.5) %>%
  mutate(direction = ifelse(meta_LFc > 0, "Up", "Down")) %>%
  count(direction)   #Down       1072
                    #Up         2132   



## —— 3) Universe and significant list ——

universe_entrez <- unique(data_mapped_clean$ENTREZID)

data_mapped_clean <- data_mapped_clean %>%
  mutate(padj = p.adjust(meta_pval, method = "BH")) 

sig_genes_entrez <- data_mapped_clean %>%
  filter(padj < 0.05 , abs(meta_LFc) > 0.5) %>%       
  pull(ENTREZID) %>%
  unique()



## ——  enrichGO (BP) + setReadable ——
ego_bp <- enrichGO(gene          = sig_genes_entrez,    
                   universe      = universe_entrez,     
                   OrgDb         = org.Hs.eg.db,        
                   keyType       = "ENTREZID",
                   ont           = "BP",                
                   pAdjustMethod = "BH",               
                   pvalueCutoff  = 0.05,                
                   qvalueCutoff  = 0.2,                 
                   readable      = FALSE)              

##  setReadable 
ego_bp_read <- setReadable(ego_bp, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")


## —— 5) Save the results / visualise 
ego_bp_df <- as.data.frame(ego_bp_read)
ego_bp_df$GeneFirst <- sapply(strsplit(ego_bp_df$geneID, "/"), `[`, 1)

write.csv(ego_bp_df, "GO_ORA_BP_results_readable2.csv", row.names = FALSE)

library(writexl)
write_xlsx(ego_bp_df, "ego_bp_df.xlsx")

names(ego_bp_df)


library(dplyr)
library(ggplot2)
library(stringr)

# enrichment results- dataframe 
df_2 <- as.data.frame(ego_bp_read) %>%
  arrange(p.adjust) %>%
  slice(1:15) %>%          # top15 pathways
  mutate(Description_wrapped = str_wrap(Description, width = 45),
         negLog10FDR = -log10(p.adjust))

df_3 <- as.data.frame(ego_bp_read) %>%
  arrange(p.adjust) %>%
  slice(1:20) %>%          # top15 pathways
  mutate(Description_wrapped = str_wrap(Description, width = 45),
         negLog10FDR = -log10(p.adjust))

# barplot + count points -top15 
g <- ggplot(df_2, aes(x = reorder(Description_wrapped, negLog10FDR),
                    y = negLog10FDR)) +
  geom_col(fill = "steelblue") +
  geom_point(aes(y = negLog10FDR, size = Count), color = "red", shape = 21, stroke = 1.2) +    # shows how many genes are active in that pathway
  coord_flip() +
  labs(title = "GO Biological Process Enrichment (Top 15)",
       x = NULL,
       y = expression(-log[10]~"(adjusted p-value)"),
       size = "Gene Count") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 13, face = "bold"))

print(g)

# save the file
ggsave("GO_BP_barplot_top15_custom.png", g, width = 8, height = 6, dpi = 300)

#color adjestment

df_2 <- df_2 %>%
  mutate(Count_cat = factor(Count)) %>%
  filter(!is.na(negLog10FDR), !is.na(Count))   

g <- ggplot(df_2, aes(x = reorder(Description_wrapped, negLog10FDR),
                      y = negLog10FDR)) +
  geom_col(fill = "steelblue") +
  geom_point(aes(color = Count_cat), size = 3) +
  scale_color_viridis_d(option = "magma") +  
  coord_flip() +
  labs(title = "GO Biological Process Enrichment (Top 15)",
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

print(g)
ggsave("GO_BP_barplot_top15_colored.png", g, width = 8, height = 6, dpi = 300)


#another way
# 1) Filter logFC 
filtered <- data_mapped_clean %>%
  filter(abs(meta_LFc) > 0.5)

# 2) Sort by p-value
filtered <- filtered %>%
  arrange(meta_pval)

# 3) Select Top10 gene
top10 <- head(filtered, 10)

library(enrichplot)
library(ggplot2)
library(stringr)

# Show only the 10 most significant terms + Plot according to GeneRatio
p <- barplot(ego_bp_read,
             showCategory = 10,
             x = "GeneRatio",           
             title = "GO Biological Process (Top 10)")

#  Highlight the labels, reduce the text size, simplify the axis titles
p <- p +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    plot.title  = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_blank()
  ) +
  labs(x = "Gene Ratio")

print(p)

# Save
ggsave("GO_BP_barplot_top10.png", plot = p, width = 7, height = 5, dpi = 300)
ggsave("GO_BP_barplot_top10.pdf", plot = p, width = 7, height = 5)


# barplot + count points -top15 
g2 <- ggplot(df_3, aes(x = reorder(Description_wrapped, negLog10FDR),
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

print(g2)


ggsave("GO_BP_barplot_top20_custom.png", g2, width = 8, height = 6, dpi = 300)



df_3 <- df_3 %>%
  mutate(Count_cat = factor(Count)) %>%
  filter(!is.na(negLog10FDR), !is.na(Count))   

g2 <- ggplot(df_3, aes(x = reorder(Description_wrapped, negLog10FDR),
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

print(g2)

