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

# conver list to data frame
meta_analysis_results <- bind_rows(
  lapply(names(results_list), function(gene) {
    cbind(`Gene Symbol` = gene, as.data.frame(results_list[[gene]]))
  })
)

write_xlsx(meta_analysis_results, "meta_analysis_results.xlsx") 

# p-value < 0.05 
meta_analysis_results_significant_genes <- meta_analysis_results %>% filter(meta_pval < 0.05) 
write_xlsx(meta_analysis_results_significant_genes, "meta_analysis_results_significant_genes.xlsx")

# In the meta-analysis script, nominal p-values were generated.
# BH-adjusted p-values (padj) and DEG filtering were applied in the enrichment/annotation script before downstream analyses.
