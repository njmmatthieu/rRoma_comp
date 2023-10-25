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
plage_t_test <- data.frame(t(apply(X = plage_output_df,
                                   MARGIN = 1,
                                   FUN = function(x) {
                                     t_test_result <- t.test(x[CF_samples], x[non_CF_samples])
                                     print(class(t_test_result))
                                     return(c(t_test_result$statistic, t_test_result$p.value))
                                   })))
colnames(plage_t_test) <- c("statistic", "p.value")

plage_t_test$statistic.abs <- abs(plage_t_test$statistic)
plage_t_test$adj.p.value <- p.adjust(plage_t_test$p.value, 
                                     method = "fdr")

plage_t_test <- plage_t_test[order(plage_t_test$statistic.abs, decreasing = T),]

write.table(plage_t_test,
          file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/project rROMA/rRoma_comp/Rehman_CF_NCF_plage_t_test_results.tsv",
          sep = "\t",
          row.names = T,
          row)



length(which(plage_t_test$adj.p.value < 0.05))
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
gsva_t_test <- data.frame(t(apply(X = gsva_output_df,
                                  MARGIN = 1,
                                  FUN = function(x) {
                                    t_test_result <- t.test(x[CF_samples], x[non_CF_samples])
                                    print(class(t_test_result))
                                    return(c(t_test_result$statistic, t_test_result$p.value))
                                  })))
colnames(gsva_t_test) <- c("statistic", "p.value")

gsva_t_test$statistic.abs <- abs(gsva_t_test$statistic)
gsva_t_test$adj.p.value <- p.adjust(gsva_t_test$p.value, 
                                    method = "fdr")

gsva_t_test <- gsva_t_test[order(gsva_t_test$statistic.abs, decreasing = T),]

write.table(gsva_t_test[,c("statistic","adj.p.value")],
            file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/project rROMA/rRoma_comp/Rehman_CF_NCF_gsva_t_test_results.tsv",
            sep = "\t",
            row.names = T)

length(which(gsva_t_test$adj.p.value < 0.05))
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
ssgsea_t_test <- data.frame(t(apply(X = ssgsea_output_df,
                                    MARGIN = 1,
                                    FUN = function(x) {
                                      t_test_result <- t.test(x[CF_samples], x[non_CF_samples])
                                      return(c(t_test_result$statistic, t_test_result$p.value))
                                    })))
colnames(ssgsea_t_test) <- c("statistic", "p.value")

ssgsea_t_test$statistic.abs <- abs(ssgsea_t_test$statistic)
ssgsea_t_test$adj.p.value <- p.adjust(ssgsea_t_test$p.value, 
                                      method = "fdr")

ssgsea_t_test <- ssgsea_t_test[order(ssgsea_t_test$statistic.abs, decreasing = T),]

write.table(ssgsea_t_test[,c("statistic","adj.p.value")],
            file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/project rROMA/rRoma_comp/Rehman_CF_NCF_ssgsea_t_test_results.tsv",
            sep = "\t",
            row.names = T)


length(which(ssgsea_t_test$adj.p.value < 0.05))
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
zscore_t_test <- data.frame(t(apply(X = zscore_output_df,
                                    MARGIN = 1,
                                    FUN = function(x) {
                                      t_test_result <- t.test(x[CF_samples], x[non_CF_samples])
                                      return(c(t_test_result$statistic, t_test_result$p.value))
                                    })))
colnames(zscore_t_test) <- c("statistic", "p.value")

zscore_t_test$statistic.abs <- abs(zscore_t_test$statistic)
zscore_t_test$adj.p.value <- p.adjust(zscore_t_test$p.value, 
                                      method = "fdr")

zscore_t_test <- zscore_t_test[order(zscore_t_test$statistic.abs, decreasing = T),]

write.table(zscore_t_test[,c("statistic","adj.p.value")],
            file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/project rROMA/rRoma_comp/Rehman_CF_NCF_zscore_t_test_results.tsv",
            sep = "\t",
            row.names = T)

length(which(zscore_t_test$adj.p.value < 0.05))
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

rROMA_t_test <- data.frame(t(apply(X = rROMA_output_df,
                                   MARGIN = 1,
                                   FUN = function(x) {
                                     t_test_result <- t.test(x[CF_samples], x[non_CF_samples])
                                     return(c(t_test_result$statistic, t_test_result$p.value))
                                   })))
colnames(rROMA_t_test) <- c("statistic", "p.value")

rROMA_t_test$statistic.abs <- abs(rROMA_t_test$statistic)
rROMA_t_test$adj.p.value <- p.adjust(rROMA_t_test$p.value, 
                                     method = "fdr")

rROMA_t_test <- rROMA_t_test[order(rROMA_t_test$statistic.abs, decreasing = T),]

write.table(rROMA_t_test[,c("statistic","adj.p.value")],
            file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/project rROMA/rRoma_comp/Rehman_CF_NCF_rROMA_t_test_results.tsv",
            sep = "\t",
            row.names = T)

length(which(rROMA_t_test$adj.p.value<0.05))
which(rownames(rROMA_t_test)=="TNFA SIGNALING VIA NFKB")

# png(filename = "./Rehman_CF_NCF_rRoma_heatmap.png", 
#     width = 1800, 
#     height = 1500, 
#     res = 150)
# rRoma_heatmap
# dev.off()

