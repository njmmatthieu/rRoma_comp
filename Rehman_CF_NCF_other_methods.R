library(GSVA)
library(ggfortify)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(tidyverse)

load(file="./Rehman_CF_NCF_RNA_exp_matrix_norm.rda")
load(file = "./Rehman_CF_NCF_sample_labels.rda")

# Sample labels

sample_labels$Type <- gsub(x = sample_labels$Type,
                            pattern = "-",
                           replacement = "_")
sample_labels$CF_sample_id <- paste(sample_labels$Type, 
                                     rep(1:6, 2),
                                     sep = " ")

CF_samples <- sample_labels %>% 
  filter(Type=="CF") %>%
  pull(CF_sample_id)

non_CF_samples <- sample_labels %>% 
  filter(Type=="non_CF") %>%
  pull(CF_sample_id)
non_CF_samples <- gsub(x = non_CF_samples,
                    pattern = "[_]",
                    replacement = "-")

# Pheatmap colors

mycolors <- colorRampPalette(brewer.pal(3, "RdBu"))(50)
pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100)
myBreaks <- seq(-1, 1, length.out = 100)

# annotation column to pheatmap
condition <- data.frame(group = sample_labels$Type)
condition$group <- gsub(pattern = "_", 
                        replacement = "-",
                        x = condition$group)
rownames(condition) <- gsub(pattern = "_", 
                            replacement = "-",
                            x = sample_labels$CF_sample_id)

# List with colors for each annotation.
condition_colors <- list(group = brewer.pal(3, "Set2"))
names(condition_colors$group) <- unique(condition$group)


expr.matrix <- as.matrix(MatData)

colnames(expr.matrix) <- sapply(colnames(expr.matrix), function(id){
  CF_id <- sample_labels %>%
    filter(sample_id==id) %>%
    pull(CF_sample_id)
  
  return(CF_id)
})

# gene.pca <- prcomp(t(expr.matrix),
#        center = T)
# 
# gene.pca.plot <- autoplot(gene.pca,
#                           data = t(expr.matrix),
#                           label = T)



# Hallmark GMT
genesets <- GSEABase::getGmt("./h.all.v2023.1.Hs.symbols.gmt")
genesets_list <- lapply(genesets, function(geneset) {
  return(geneset@geneIds)
})
names(genesets_list) <- lapply(genesets, function(geneset) {
  geneset_name <- geneset@setName
  geneset_name <- gsub(pattern = "HALLMARK_", 
                       replacement = "", 
                       x = geneset_name)
  geneset_name <- gsub(pattern = "_", 
                       replacement = " ", 
                       x = geneset_name)
  return(geneset_name)
})

# PLAGE

plage_es <- gsva(expr.matrix, 
                genesets_list, 
                verbose=FALSE,
                method="plage")

plage_output_df <- data.frame(plage_es)
colnames(plage_output_df) <- gsub(pattern = "_",
                                  replacement = "-",
                                  x = colnames(plage_output_df))
colnames(plage_output_df) <- gsub(pattern = "[.]",
                                  replacement = " ",
                                  x = colnames(plage_output_df))

plage_heatmap <- pheatmap(plage_output_df,
         annotation_col = condition,
         annotation_colors = condition_colors,
         color = pal,
         breaks = myBreaks,
         border_color = "NA")

# png(filename = "./Rehman_CF_NCF_plage_heatmap.png", 
#     width = 1800, 
#     height = 1500, 
#     res = 150)
# plage_heatmap
# dev.off()



# Perform t-tests for each gene
plage_t_test_p_values <- apply(X = plage_output_df,
                               MARGIN = 1,
                               FUN = function(x) {
                                 t_test_result <- t.test(x[CF_samples], x[non_CF_samples])
                                 return(t_test_result$p.value)
                                 })

plage_adjusted_p_values <- p.adjust(plage_t_test_p_values, 
                              method = "fdr")

length(plage_adjusted_p_values[plage_adjusted_p_values < 0.05])
length(plage_adjusted_p_values < 0.05)

# GSVA

gsva_es <- gsva(expr.matrix, 
                 genesets_list, 
                 verbose=FALSE,
                 method="gsva",
                 kcdf="Poisson")

gsva_output_df <- data.frame(gsva_es)
colnames(gsva_output_df) <- gsub(pattern = "_",
                                  replacement = "-",
                                  x = colnames(gsva_output_df))
colnames(gsva_output_df) <- gsub(pattern = "[.]",
                                  replacement = " ",
                                  x = colnames(gsva_output_df))

gsva_heatmap <- pheatmap(gsva_output_df,
         annotation_col = condition,
         annotation_colors = condition_colors,
         color = pal,
         breaks = myBreaks,
         border_color = "NA")

# png(filename = "./Rehman_CF_NCF_gsva_heatmap.png", 
#     width = 1800, 
#     height = 1500, 
#     res = 150)
# gsva_heatmap
# dev.off()

# Perform t-tests for each gene
gsva_t_test_p_values <- apply(X = gsva_output_df,
                         MARGIN = 1,
                         FUN = function(x) {
                           # print(x[CF_samples])
                           # print(x[non_CF_samples])
                           t_test_result <- t.test(x[CF_samples], x[non_CF_samples])
                           return(t_test_result$p.value)
                         })

gsva_adjusted_p_values <- p.adjust(gsva_t_test_p_values, 
                              method = "fdr")

length(gsva_adjusted_p_values[gsva_adjusted_p_values < 0.05])
table(gsva_adjusted_p_values < 0.05)

pheatmap(gsva_output_df)

## SSGSEA

ssgsea_es <- gsva(expr.matrix, 
                genesets_list, 
                verbose=FALSE,
                method="ssgsea")

ssgsea_output_df <- data.frame(ssgsea_es)
colnames(ssgsea_output_df) <- gsub(pattern = "_",
                                 replacement = "-",
                                 x = colnames(ssgsea_output_df))
colnames(ssgsea_output_df) <- gsub(pattern = "[.]",
                                 replacement = " ",
                                 x = colnames(ssgsea_output_df))

ssgsea_heatmap <- pheatmap(ssgsea_output_df,
         annotation_col = condition,
         annotation_colors = condition_colors,
         color = pal,
         breaks = myBreaks,
         border_color = "NA")

# png(filename = "./Rehman_CF_NCF_ssgsea_heatmap.png", 
#     width = 1800, 
#     height = 1500, 
#     res = 150)
# ssgsea_heatmap
# dev.off()

# Perform t-tests for each gene
ssgsea_t_test_p_values <- apply(X = ssgsea_output_df,
                              MARGIN = 1,
                              FUN = function(x) {
                                # print(x[CF_samples])
                                # print(x[non_CF_samples])
                                t_test_result <- t.test(x[CF_samples], x[non_CF_samples])
                                return(t_test_result$p.value)
                              })

ssgsea_adjusted_p_values <- p.adjust(ssgsea_t_test_p_values, 
                                   method = "fdr")

names(ssgsea_adjusted_p_values[ssgsea_adjusted_p_values < 0.05])
table(ssgsea_adjusted_p_values < 0.05)

pheatmap(ssgsea_output_df)

## Zscore

zscore_es <- gsva(expr.matrix, 
                  genesets_list, 
                  verbose=FALSE,
                  method="zscore")

zscore_output_df <- data.frame(zscore_es)
colnames(zscore_output_df) <- gsub(pattern = "_",
                                   replacement = "-",
                                   x = colnames(zscore_output_df))
colnames(zscore_output_df) <- gsub(pattern = "[.]",
                                   replacement = " ",
                                   x = colnames(zscore_output_df))

zscore_mycolors <- colorRampPalette(brewer.pal(3, "RdBu"))(50)
zscore_pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "PiYG")))(100)
zscore_myBreaks <- seq(min(zscore_output_df), max(zscore_output_df), length.out = 100)

zscore_heatmap <- pheatmap(zscore_output_df,
         annotation_col = condition,
         annotation_colors = condition_colors,
         color = zscore_pal,
         breaks = zscore_myBreaks,
         border_color = "NA")

# png(filename = "./Rehman_CF_NCF_zscore_heatmap.png", 
#     width = 1800, 
#     height = 1500, 
#     res = 150)
# zscore_heatmap
# dev.off()

# Perform t-tests for each gene
zscore_t_test_p_values <- apply(X = zscore_output_df,
                                MARGIN = 1,
                                FUN = function(x) {
                                  # print(x[CF_samples])
                                  # print(x[non_CF_samples])
                                  t_test_result <- t.test(x[CF_samples], x[non_CF_samples])
                                  return(t_test_result$p.value)
                                })

zscore_adjusted_p_values <- p.adjust(zscore_t_test_p_values, 
                                     method = "fdr")

length(zscore_adjusted_p_values[zscore_adjusted_p_values < 0.05])
table(zscore_adjusted_p_values < 0.05)

pheatmap(zscore_output_df)

## rROMA

rROMA_output_df <- read.table(file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/project rROMA/Reviews/rROMA_output_SampleMatrix_2023_09_19.csv",
                              sep = "\t")

colnames(rROMA_output_df) <- sapply(colnames(rROMA_output_df), function(id){
  CF_id <- sample_labels %>%
    filter(sample_id==id) %>%
    pull(CF_sample_id)
  
  return(CF_id)
})
colnames(rROMA_output_df) <- gsub(pattern = "_",
                                   replacement = "-",
                                   x = colnames(rROMA_output_df))
colnames(rROMA_output_df) <- gsub(pattern = "[.]",
                                   replacement = " ",
                                   x = colnames(rROMA_output_df))

rRoma_heatmap <- pheatmap(rROMA_output_df,
                           annotation_col = condition,
                           annotation_colors = condition_colors,
                           color = pal,
                           breaks = myBreaks,
                           border_color = "NA")

# Perform t-tests for each gene
rROMA_t_test_p_values <- apply(X = rROMA_output_df,
                                MARGIN = 1,
                                FUN = function(x) {
                                  # print(x[CF_samples])
                                  # print(x[non_CF_samples])
                                  t_test_result <- t.test(x[CF_samples], x[non_CF_samples])
                                  return(t_test_result$p.value)
                                })

rROMA_adjusted_p_values <- p.adjust(rROMA_t_test_p_values, 
                                     method = "fdr")

length(rROMA_adjusted_p_values[rROMA_adjusted_p_values < 0.05])

# png(filename = "./Rehman_CF_NCF_rRoma_heatmap.png", 
#     width = 1800, 
#     height = 1500, 
#     res = 150)
# rRoma_heatmap
# dev.off()

