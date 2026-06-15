#metanalysis_withcov
library(dplyr)
library(metafor)
library(readxl)
library(writexl)

metaanalyses_withcov <- read_excel("/Users/ebraraltinalan/Desktop/GEO/final_withcov_geo/files/merged_cleaned_data_withcov.xlsx")

#only genes that appear in at least three datasets were taken into analysis

merged_cleaned_filtered_withcov <- metaanalyses_withcov %>%
  group_by(`Gene Symbol`) %>%                        
  filter(n_distinct(Study) >= 3) %>%               
  ungroup() 

# meta analysis: Random effects for each gene  Sidik-Jonkman estimator/ heterogeneity
library(dplyr)
library(metafor)
library(readxl)
library(tibble)
library(writexl)

# create a unique_genes list 
unique_genes_withcov <- unique(merged_cleaned_filtered_withcov$`Gene Symbol`)

# create a list for save the results
results_list_withcov <- list()

# x is iterated in a loop/function
for (x in unique_genes_withcov) {
  tmp_df_genes_with <- merged_cleaned_filtered_withcov[merged_cleaned_filtered_withcov$`Gene Symbol` == x, ]
  
  # at least three study
  if (nrow(tmp_df_genes_with) >= 3) {
    tmp_meta_model <- rma.uni(
      yi = tmp_df_genes_with$meanLogFC,
      vi = tmp_df_genes_with$maxSE^2,
      method = "SJ"
    )
    
    
 # add the results to list
    results_list_withcov[[x]] <- list(
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


# Convert list to data frame
meta_analysis_results_withcov <- bind_rows(
  lapply(names(results_list_withcov), function(gene) {
    cbind(`Gene Symbol` = gene, as.data.frame(results_list_withcov[[gene]]))
  })
)

# Ensure numeric columns are correctly formatted
meta_analysis_results_withcov <- meta_analysis_results_withcov %>%
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
meta_analysis_results_withcov <- meta_analysis_results_withcov %>%
  mutate(
    padj = p.adjust(meta_pval, method = "BH")
  )

# Save full meta-analysis results including nominal and adjusted p-values
write_xlsx(
  meta_analysis_results_withcov,
  "meta_analysis_results_withcov.xlsx"
)

# Nominally significant genes: meta_pval < 0.05
meta_analysis_results_withcov_nominal_significant <- meta_analysis_results_withcov %>%
  filter(meta_pval < 0.05)

write_xlsx(
  meta_analysis_results_withcov_nominal_significant,
  "meta_analysis_results_withcov_nominal_pval_0.05.xlsx"
)

# FDR-significant genes: padj < 0.05
meta_analysis_results_withcov_FDR_significant <- meta_analysis_results_withcov %>%
  filter(padj < 0.05)

write_xlsx(
  meta_analysis_results_withcov_FDR_significant,
  "meta_analysis_results_withcov_padj_0.05.xlsx"
)

# DEGs using the primary threshold: |metaLFC| > 0.5 and padj < 0.05
DEGs_withcov_metaLFC_0.5_padj <- meta_analysis_results_withcov %>%
  filter(abs(meta_LFc) > 0.5, padj < 0.05) %>%
  mutate(
    regulation = case_when(
      meta_LFc > 0 ~ "Upregulated",
      meta_LFc < 0 ~ "Downregulated",
      TRUE ~ "No direction"
    )
  )

write_xlsx(
  DEGs_withcov_metaLFC_0.5_padj,
  "DEGs_withcov_abs_metaLFC_0.5_padj_0.05.xlsx"
)

# More stringent subset: |metaLFC| > 1 and padj < 0.05
DEGs_withcov_metaLFC_1_padj <- meta_analysis_results_withcov %>%
  filter(abs(meta_LFc) > 1, padj < 0.05) %>%
  mutate(
    regulation = case_when(
      meta_LFc > 0 ~ "Upregulated",
      meta_LFc < 0 ~ "Downregulated",
      TRUE ~ "No direction"
    )
  )

write_xlsx(
  DEGs_withcov_metaLFC_1_padj,
  "DEGs_withcov_abs_metaLFC_1_padj_0.05.xlsx"
)


# Summary counts
summary_counts_withcov <- tibble::tibble(
  Criterion = c(
    "Total genes",
    "meta_pval < 0.05",
    "padj < 0.05",
    "|metaLFC| > 0.5 and padj < 0.05",
    "|metaLFC| > 1 and padj < 0.05",
    "|metaLFC| > 1 and meta_pval < 0.05"
  ),
  Total = c(
    nrow(meta_analysis_results_withcov),
    nrow(meta_analysis_results_withcov_nominal_significant),
    nrow(meta_analysis_results_withcov_FDR_significant),
    nrow(DEGs_withcov_metaLFC_0.5_padj),
    nrow(DEGs_withcov_metaLFC_1_padj),
    nrow(genes_withcov_metaLFC_1_nominal)
  )
)

write_xlsx(
  summary_counts_withcov,
  "summary_counts_withcov_meta_analysis.xlsx"
)

# Up/down counts for reporting
updown_counts_withcov <- list(
  DEG_metaLFC_0.5_padj = DEGs_withcov_metaLFC_0.5_padj %>% count(regulation),
  DEG_metaLFC_1_padj = DEGs_withcov_metaLFC_1_padj %>% count(regulation),
  metaLFC_1_nominal = genes_withcov_metaLFC_1_nominal %>% count(regulation)
)

write_xlsx(
  updown_counts_withcov,
  "updown_counts_withcov_meta_analysis.xlsx"
)

summary_counts_withcov
updown_counts_withcov
