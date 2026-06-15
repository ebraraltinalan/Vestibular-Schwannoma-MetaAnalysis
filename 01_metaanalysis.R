#meta analyses
install.packages("metafor", type = "binary", repos = "https://cloud.r-project.org")
library(metafor)
packageVersion("metafor")

library(dplyr)
library(metafor)
library(readxl)
library(writexl)

metaanalyses2 <- read_excel("/Users/ebraraltinalan/Desktop/merged_cleaned_data.xlsx")

library(dplyr)

#only genes that appear in at least three datasets were taken into analysis

merged_cleaned_filtered <- metaanalyses2 %>%
  group_by(`Gene Symbol`) %>%                       
  filter(n_distinct(Study) >= 3) %>%                  
  ungroup()                                       


# meta analysis: Random effects for each gene  Sidik-Jonkman estimator/ heterogeneity
library(dplyr)
library(metafor)
library(readxl)
library(tibble)
library(writexl)

# create a unique_genes list #It does not filter the rows. It only creates a list based on the gene symbol for the loop.
unique_genes <- unique(merged_cleaned_filtered$`Gene Symbol`)

# create a list for save the results
results_list <- list()

# x is iterated in a loop/function
for (x in unique_genes) {
  tmp_df_genes <- merged_cleaned_filtered[merged_cleaned_filtered$`Gene Symbol` == x, ]
  
  # at least three study
  if (nrow(tmp_df_genes) >= 3) {
    tmp_meta_model <- rma.uni(      
      yi = tmp_df_genes$meanLogFC,
      vi = tmp_df_genes$maxSE^2,
      method = "SJ"     #Sidik–Jonkman
    )
    
    # add the results to list
    results_list[[x]] <- list(
      meta_LFc = tmp_meta_model$b,
      meta_se = tmp_meta_model$se,
      meta_pval = tmp_meta_model$pval,
      tau2 = tmp_meta_model$tau2,
      I2 = tmp_meta_model$I2,
      H2 = tmp_meta_model$H2,
      Q = tmp_meta_model$QE,
      Q.p = tmp_meta_model$QEp
    )
  }
}

# convert list to data frame
meta_analysis_results <- bind_rows(
  lapply(names(results_list), function(gene) {
    cbind(`Gene Symbol` = gene, as.data.frame(results_list[[gene]]))
  })
)

# Ensure numeric columns are correctly formatted
meta_analysis_results <- meta_analysis_results %>%
  mutate(
    meta_LFc = as.numeric(meta_LFc),
    meta_se = as.numeric(meta_se),
    meta_pval = as.numeric(meta_pval),
    tau2 = as.numeric(tau2),
    I2 = as.numeric(I2),
    H2 = as.numeric(H2),
    Q = as.numeric(Q),
    Q.p = as.numeric(Q.p)
  )

# Calculate Benjamini-Hochberg adjusted p-values
# Multiple-testing correction was performed across all genes included in the meta-analysis,
# not only among nominally significant genes.
meta_analysis_results <- meta_analysis_results %>%
  mutate(
    padj = p.adjust(meta_pval, method = "BH")
  )

# Save full meta-analysis results including nominal and adjusted p-values
write_xlsx(meta_analysis_results, "meta_analysis_results_with_padj.xlsx")

# Nominally significant genes: meta_pval < 0.05
meta_analysis_results_nominal_significant <- meta_analysis_results %>%
  filter(meta_pval < 0.05)

write_xlsx(
  meta_analysis_results_nominal_significant,
  "meta_analysis_results_nominal_pval_0.05.xlsx"
)

# FDR-significant genes: padj < 0.05
meta_analysis_results_FDR_significant <- meta_analysis_results %>%
  filter(padj < 0.05)

write_xlsx(
  meta_analysis_results_FDR_significant,
  "meta_analysis_results_padj_0.05.xlsx"
)

# DEGs using the primary threshold: |metaLFC| > 0.5 and padj < 0.05
DEGs_metaLFC_0.5_padj <- meta_analysis_results %>%
  filter(abs(meta_LFc) > 0.5, padj < 0.05) %>%
  mutate(
    regulation = case_when(
      meta_LFc > 0 ~ "Upregulated",
      meta_LFc < 0 ~ "Downregulated",
      TRUE ~ "No direction"
    )
  )

write_xlsx(
  DEGs_metaLFC_0.5_padj,
  "DEGs_abs_metaLFC_0.5_padj_0.05.xlsx"
)

# More stringent subset: |metaLFC| > 1 and padj < 0.05
DEGs_metaLFC_1_padj <- meta_analysis_results %>%
  filter(abs(meta_LFc) > 1, padj < 0.05) %>%
  mutate(
    regulation = case_when(
      meta_LFc > 0 ~ "Upregulated",
      meta_LFc < 0 ~ "Downregulated",
      TRUE ~ "No direction"
    )
  )

write_xlsx(
  DEGs_metaLFC_1_padj,
  "DEGs_abs_metaLFC_1_padj_0.05.xlsx"
)

# Optional: genes used for nominal volcano-style threshold
genes_metaLFC_1_nominal <- meta_analysis_results %>%
  filter(abs(meta_LFc) > 1, meta_pval < 0.05) %>%
  mutate(
    regulation = case_when(
      meta_LFc > 0 ~ "Upregulated",
      meta_LFc < 0 ~ "Downregulated",
      TRUE ~ "No direction"
    )
  )

write_xlsx(
  genes_metaLFC_1_nominal,
  "genes_abs_metaLFC_1_nominal_meta_pval_0.05.xlsx"
)

# Summary counts
summary_counts <- tibble::tibble(
  Criterion = c(
    "Total genes",
    "meta_pval < 0.05",
    "padj < 0.05",
    "|metaLFC| > 0.5 and padj < 0.05",
    "|metaLFC| > 1 and padj < 0.05",
    "|metaLFC| > 1 and meta_pval < 0.05"
  ),
  Total = c(
    nrow(meta_analysis_results),
    nrow(meta_analysis_results_nominal_significant),
    nrow(meta_analysis_results_FDR_significant),
    nrow(DEGs_metaLFC_0.5_padj),
    nrow(DEGs_metaLFC_1_padj),
    nrow(genes_metaLFC_1_nominal)
  )
)

write_xlsx(summary_counts, "summary_counts_meta_analysis.xlsx")

# Up/down counts for reporting
DEGs_metaLFC_0.5_padj %>%
  count(regulation)

DEGs_metaLFC_1_padj %>%
  count(regulation)

genes_metaLFC_1_nominal %>%
  count(regulation)

summary_counts
write_xlsx(meta_analysis_results, "meta_analysis_results.xlsx") 


# Nominal meta-analysis p-values were generated using the random-effects model.
# BH-adjusted p-values (padj) were calculated across all genes included in the meta-analysis.
# DEG filtering was performed using |metaLFC| thresholds together with padj < 0.05.
