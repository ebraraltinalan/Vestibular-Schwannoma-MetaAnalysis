# ============================================================
# 09_DrugBank_network_visualization.R
# Drug–gene interaction network visualization
# ============================================================

library(readxl)
library(dplyr)
library(igraph)
library(ggraph)
library(ggplot2)

# ------------------------------------------------------------
# 1. Read filtered DrugBank interaction table
# ------------------------------------------------------------

drug_gene_interactions_filtered_metalfc1 <- read_excel(
  "data/drug_gene_interactions_filtered_metalfc1.xlsx"
)

# ------------------------------------------------------------
# 2. Select top 10 downregulated and upregulated genes
# ------------------------------------------------------------

down_top10_genes <- drug_gene_interactions_filtered_metalfc1 %>%
  filter(regulation == "downregulated") %>%
  distinct(`Gene Symbol`, meta_LFc) %>%
  arrange(meta_LFc) %>%
  slice(1:10) %>%
  pull(`Gene Symbol`)

df_down_top10 <- drug_gene_interactions_filtered_metalfc1 %>%
  filter(`Gene Symbol` %in% down_top10_genes)

up_top10_genes <- drug_gene_interactions_filtered_metalfc1 %>%
  filter(regulation == "upregulated") %>%
  distinct(`Gene Symbol`, meta_LFc) %>%
  arrange(desc(meta_LFc)) %>%
  slice(1:10) %>%
  pull(`Gene Symbol`)

df_up_top10 <- drug_gene_interactions_filtered_metalfc1 %>%
  filter(`Gene Symbol` %in% up_top10_genes)

# ------------------------------------------------------------
# 3. Function for network plotting
# ------------------------------------------------------------

plot_gene_drug_network <- function(df_sub) {
  
  edges <- df_sub %>%
    select(from = `Gene Symbol`, to = name)
  
  genes <- df_sub %>%
    distinct(name = `Gene Symbol`, meta_LFc, regulation)
  
  drugs <- df_sub %>%
    distinct(name = name) %>%
    mutate(
      meta_LFc = NA,
      regulation = "drug"
    )
  
  nodes <- bind_rows(genes, drugs) %>%
    mutate(
      regulation = factor(
        regulation,
        levels = c("upregulated", "downregulated", "drug")
      )
    )
  
  graph <- graph_from_data_frame(
    d = edges,
    vertices = nodes,
    directed = FALSE
  )
  
  ggraph(graph, layout = "fr") +
    geom_edge_link(color = "gray70", alpha = 0.8) +
    geom_node_point(aes(color = regulation), size = 7) +
    geom_node_text(
      aes(
        label = name,
        fontface = ifelse(regulation == "drug", "plain", "bold")
      ),
      repel = TRUE,
      max.overlaps = Inf,
      size = 3
    ) +
    scale_color_manual(
      values = c(
        "upregulated" = "#D6EAF8",
        "downregulated" = "#FDEBD0",
        "drug" = "#F1948A"
      ),
      breaks = c("upregulated", "downregulated", "drug"),
      labels = c("Upregulated", "Downregulated", "Drugname"),
      name = NULL
    ) +
    theme_void() +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 6)
    )
}

# ------------------------------------------------------------
# 4. Generate and save figures
# ------------------------------------------------------------

fig3_downregulated <- plot_gene_drug_network(df_down_top10)
fig4_upregulated <- plot_gene_drug_network(df_up_top10)

ggsave(
  filename = "figures/Fig3_DrugBank_downregulated_network.pdf",
  plot = fig3_downregulated,
  width = 12,
  height = 6
)

ggsave(
  filename = "figures/Fig4_DrugBank_upregulated_network.tiff",
  plot = fig4_upregulated,
  width = 12,
  height = 6,
  dpi = 300
)