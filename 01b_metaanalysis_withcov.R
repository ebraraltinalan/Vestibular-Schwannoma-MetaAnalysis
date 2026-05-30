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


# conver list to data frame
meta_analysis_results_withcov <- bind_rows(
  lapply(names(results_list_withcov), function(gene) {
    cbind(`Gene Symbol` = gene, as.data.frame(results_list_withcov[[gene]]))
  })
)  
#21113entires

meta_analysis_results_withcov <- meta_analysis_results_withcov %>%
  mutate(padj = p.adjust(meta_pval, method = "BH"))



write_xlsx(meta_analysis_results_withcov, "meta_analysis_results_withcov.xlsx") 


meta_analysis_results_withcov_signif_genes <- meta_analysis_results_withcov %>% 
filter(padj < 0.05, abs(meta_LFc) > 0.5)

write_xlsx(meta_analysis_results_withcov_signif_genes, "meta_analysis_results_withcov_signif_genes.xlsx")
