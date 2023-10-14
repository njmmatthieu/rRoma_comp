library(GSVA)
library(ggfortify)
library(ggplot2)
library(tidyr)
library(tidyverse)

load(file="/Users/matthieu/ownCloud/Thèse/Systems Biology/project rROMA/Reviews/Rehman_exp_RNA_matrix_norm.Rda")
load(file ="/Users/matthieu/ownCloud/Thèse/Systems Biology/project rROMA/Reviews/Rehman_exp_sample_labels.Rda")

control_samples <- Rehman_exp_coldata %>% 
  filter(Treatment=="control") %>%
  pull(Simple_id)

inf_samples <- Rehman_exp_coldata %>% 
  filter(Treatment=="IL-17+TNFa") %>%
  pull(Simple_id) 
inf_samples <- gsub(x = inf_samples,
                    pattern = "[-|+]",
                    replacement = ".")

# expr.matrix <- as.matrix(MatData)
# 
# colnames(Rehman_exp_RNA_matrix_norm) <- sapply(colnames(Rehman_exp_RNA_matrix_norm), function(id){
#   Treatment_id <- Rehman_exp_coldata %>%
#     filter(sample_id==id) %>%
#     pull(Simple_id)
#   
#   return(Treatment_id)
# })

# Hallmark Gene set
hallmark_genesets <- GSEABase::getGmt("./h.all.v2023.1.Hs.symbols.gmt")
genesets_list <- lapply(genesets, function(geneset) {
  return(geneset@geneIds)
})
names(genesets_list) <- lapply(genesets, function(geneset) {
  return(geneset@setName)
})

# PLAGE

plage_es <- gsva(Rehman_exp_RNA_matrix_norm, 
                 genesets_list, 
                 verbose=FALSE,
                 method="plage")

plage_output_df <- data.frame(plage_es)

# Perform t-tests for each gene
plage_t_test_p_values <- apply(X = plage_output_df,
                               MARGIN = 1,
                               FUN = function(x) {
                                 t_test_result <- t.test(x[control_samples], x[inf_samples])
                                 return(t_test_result$p.value)
                               })

plage_adjusted_p_values <- p.adjust(plage_t_test_p_values, 
                                    method = "fdr")

names(plage_adjusted_p_values[plage_adjusted_p_values < 0.05])
length(plage_adjusted_p_values < 0.05)

plage_adjusted_p_values[order(plage_adjusted_p_values, decreasing = F)]

pheatmap(plage_output_df)

# GSVA

gsva_es <- gsva(Rehman_exp_RNA_matrix_norm, 
                genesets_list, 
                verbose=FALSE,
                method="gsva",
                kcdf="Poisson")

gsva_output_df <- data.frame(gsva_es)

# Perform t-tests for each gene
gsva_t_test_p_values <- apply(X = gsva_output_df,
                              MARGIN = 1,
                              FUN = function(x) {
                                t_test_result <- t.test(x[control_samples], x[inf_samples])
                                return(t_test_result$p.value)
                              })

gsva_adjusted_p_values <- p.adjust(gsva_t_test_p_values, 
                                   method = "fdr")

names(gsva_adjusted_p_values[gsva_adjusted_p_values < 0.05])
table(gsva_adjusted_p_values < 0.05)

gsva_adjusted_p_values[order(gsva_adjusted_p_values, decreasing = F)]

pheatmap(gsva_output_df)

## SSGSEA

ssgsea_es <- gsva(Rehman_exp_RNA_matrix_norm, 
                  genesets_list, 
                  verbose=FALSE,
                  method="ssgsea")

ssgsea_output_df <- data.frame(ssgsea_es)

# Perform t-tests for each gene
ssgsea_t_test_p_values <- apply(X = ssgsea_output_df,
                                MARGIN = 1,
                                FUN = function(x) {
                                  # print(x[CF_samples])
                                  # print(x[non_CF_samples])
                                  t_test_result <- t.test(x[control_samples], x[inf_samples])
                                  return(t_test_result$p.value)
                                })

ssgsea_adjusted_p_values <- p.adjust(ssgsea_t_test_p_values, 
                                     method = "fdr")

names(ssgsea_adjusted_p_values[ssgsea_adjusted_p_values < 0.05])
table(ssgsea_adjusted_p_values < 0.05)

ssgsea_adjusted_p_values[order(ssgsea_adjusted_p_values, decreasing = F)]

pheatmap(ssgsea_output_df)

## Zscore

zscore_es <- gsva(Rehman_exp_RNA_matrix_norm, 
                  genesets_list, 
                  verbose=FALSE,
                  method="zscore")

zscore_output_df <- data.frame(zscore_es)

# Perform t-tests for each gene
zscore_t_test_p_values <- apply(X = zscore_output_df,
                                MARGIN = 1,
                                FUN = function(x) {
                                  # print(x[CF_samples])
                                  # print(x[non_CF_samples])
                                  t_test_result <- t.test(x[control_samples], 
                                                          x[inf_samples])
                                  return(t_test_result$p.value)
                                })

zscore_adjusted_p_values <- p.adjust(zscore_t_test_p_values, 
                                     method = "fdr")

names(zscore_adjusted_p_values[zscore_adjusted_p_values < 0.05])
table(zscore_adjusted_p_values < 0.05)

zscore_adjusted_p_values[order(zscore_adjusted_p_values, decreasing = F)]

pheatmap(zscore_output_df)

