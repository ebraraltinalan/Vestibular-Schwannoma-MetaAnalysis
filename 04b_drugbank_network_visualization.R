#drugbank_withcovariates

library(dbparser)
library(dplyr)
library(readr)
library(xml2)


# Extract target and polypeptide data
targets_withcov <- dv$cett$targets
polypeptides_withcov <- targets_withcov$polypeptides$general_information
synonyms_withcov <- dv$cett$targets$polypeptides$synonyms

# Replace gene symbols
symbol_map2 <- setNames(Homo_sapiens_gene_info$V3, Homo_sapiens_gene_info$V5)
polypeptides_withcov$gene_name <- ifelse(
  polypeptides_withcov$gene_name %in% names(symbol_map2),
  symbol_map2[polypeptides_withcov$gene_name],
  polypeptides_withcov$gene_name
)

targets_general2 <- targets_withcov$general_information  # tibble
nrow(targets_general2)   #23136

drug_names2 <- drugs_general[, c("drugbank_id", "name")]

nrow(drugs_general)
nrow(synonyms_withcov)
nrow(targets_general2)
nrow(polypeptides_withcov)

#merge tables 
library(dplyr)

# The targets will merge with polypeptides (by target_id)
targets_merged2 <- targets_general2 %>%
  left_join(polypeptides_withcov, by = "target_id") 

# Merge targets + polypeptides with synonym (by target_id)
targets_merged2 <- targets_merged2 %>%
  left_join(synonyms_withcov, by = "target_id")     


# Add drugs informations (by drugbank_id)
final_merged2 <- targets_merged2 %>%
  left_join(drugs_general[, c("drugbank_id", "name")], by = "drugbank_id")

# Filter   the “Known action”  interactions
final_verified2 <- final_merged2 %>%
  filter(known_action == "yes")

# Control
nrow(final_verified2)
head(final_verified2$known_action)

library(readxl)

# 
meta2 <- read_excel("/Users/ebraraltinalan/Desktop/GEO/final_withcov_geo/files/meta_analysis_results_withcov_sig_genes_names.xlsx")

# Control
head(meta2)
colnames(meta2)

#categorize by up-regulated and down regulated
meta2 <- meta2 %>%
  mutate(regulation = case_when(
    meta_LFc > 0  ~ "up-regulated",
    meta_LFc < 0  ~ "down-regulated",
    TRUE          ~ "no_change"
  ))


colnames(final_verified2)



# Convert gene names to uppercase
meta2 <- meta2 %>%
  mutate(`Gene Symbol` = toupper(`Gene Symbol`)) 

final_verified2 <- final_verified2 %>%
  mutate(gene_name = toupper(gene_name))

final_verified_unique2 <- final_verified2 %>%
  distinct(gene_name, drugbank_id, .keep_all = TRUE)


# Merge
drug_gene_interactions2 <- inner_join(
  meta2,
  final_verified_unique2,
  by = c("Gene Symbol" = "gene_name")
)
head(drug_gene_interactions2)

library(writexl)
write_xlsx(drug_gene_interactions2,"drug_gene_interactions_unique_withcov.xlsx")


# 
drug_gene_interactions_unique2 <- drug_gene_interactions2 %>%
  select(
    `Gene Symbol`,
    GeneName,
    regulation,
    meta_LFc,
    meta_pval,
    drugbank_id,
    name,           # Drug primary name
    target_id,
    known_action,
    general_function,
    specific_function,
    synonym
  )

write_xlsx(drug_gene_interactions_unique2, "drug_gene_interactions_unique_filtered2.xlsx")

# Control
head(drug_gene_interactions_unique2)
nrow(drug_gene_interactions_unique2)



#Visualisation

#build drug-gene network (igraph + ggraph)

library(dplyr)
library(igraph)
library(ggraph)

# ---  Top 15 gene ---
top_genes20_withcov <- drug_gene_interactions_unique2 %>%
  arrange(desc(abs(meta_LFc))) %>%
  distinct(`Gene Symbol`, .keep_all = TRUE) %>%
  slice(1:20)

# --
edges_top2 <- drug_gene_interactions_unique2 %>%
  filter(`Gene Symbol` %in% top_genes20_withcov$`Gene Symbol`) %>%
  select(`Gene Symbol`, name, general_function)

colnames(edges_top2) <- c("from", "to", "general_function")

# --- gene nodes
genes_top2 <- unique(top_genes20_withcov[, c('Gene Symbol', "meta_LFc")])
genes_top2$regulation <- ifelse(genes_top2$meta_LFc > 0, "upregulated", "downregulated")
colnames(genes_top2)[1] <- "name"

# --- drug nodes
drugs_top2 <- unique(data.frame(
  name = edges_top2$to,
  regulation = "drug",
  meta_LFc = NA
))

# --- merge the nodes
nodes_top2 <- rbind(genes_top2, drugs_top2)

# --- graph
graph_top2 <- graph_from_data_frame(d = edges_top2, vertices = nodes_top2, directed = TRUE)

# ---
set.seed(42)  
ggraph(graph_top2, layout = "fr") +
  geom_edge_link(alpha = 0.9, color = "gray70") +
  
  # nodes
  geom_node_point(aes(color = regulation), size = 8) +
  
  # labels
  geom_node_text(
    aes(label = name),
    size = 2.3,  
    color = "black",
    fontface = ifelse(nodes_top2$regulation == "drug", "plain", "bold"),
    repel = TRUE,
    max.overlaps = Inf
  ) +
  
  # Legend 
  scale_color_manual(
    values = c("upregulated" = "#D6EAF8",
               "downregulated" = "#FDEBD0",
               "drug" = "#F1948A"),
    breaks = c("upregulated", "downregulated", "drug"),
    labels = c("Upregulated", "Downregulated", "Drugname"),
    name = NULL
  ) +
  
  # theme
  theme_void() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 6)
  )

# --- 8. Save the graph ---
ggsave("DrugBank_network_withcov_Top20.tiff", plot = last_plot(), width = 12, height = 8, dpi = 300, bg = "white")

