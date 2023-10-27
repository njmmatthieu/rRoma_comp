library(GSVA)
library(ggfortify)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(tidyverse)

load(file="./Rehman_IL17_TNFa_RNA_exp_matrix_norm.Rda")
load(file ="./Rehman_IL17_TNFa_sample_labels.Rda")

Rehman_exp_coldata$Simple_id <- gsub(pattern = "_",
                                     replacement = " ",
                                     x = Rehman_exp_coldata$Simple_id)
colnames(Rehman_exp_RNA_matrix_norm) <- Rehman_exp_coldata$Simple_id

control_samples <- Rehman_exp_coldata %>% 
  filter(Treatment=="control") %>%
  pull(Simple_id)

inf_samples <- Rehman_exp_coldata %>% 
  filter(Treatment=="IL-17+TNFa") %>%
  pull(Simple_id) 
# inf_samples <- gsub(x = inf_samples,
#                     pattern = "[_]",
#                     replacement = " ")
# inf_samples <- gsub(x = inf_samples,
#                     pattern = "[-|+]",
#                     replacement = ".")

# expr.matrix <- as.matrix(MatData)
# 
# colnames(Rehman_exp_RNA_matrix_norm) <- sapply(colnames(Rehman_exp_RNA_matrix_norm), function(id){
#   Treatment_id <- Rehman_exp_coldata %>%
#     filter(sample_id==id) %>%
#     pull(Simple_id)
#   
#   return(Treatment_id)
# })

# Pheatmap colors

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100)
myBreaks <- seq(-1, 1, length.out = 100)

# annotation column to pheatmap
condition <- data.frame(group = Rehman_exp_coldata$Treatment)
# condition$group <- gsub(pattern = "_", 
#                         replacement = "-",
#                         x = condition$group)
rownames(condition) <- gsub(pattern = "_", 
                            replacement = " ",
                            x = Rehman_exp_coldata$Simple_id)

# List with colors for each annotation.
condition_colors <- list(group = brewer.pal(3, "Set3"))
names(condition_colors$group) <- unique(condition$group)

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

plage_es <- gsva(Rehman_exp_RNA_matrix_norm, 
                 genesets_list, 
                 verbose=FALSE,
                 method="plage")

plage_output_df <- data.frame(plage_es)
colnames(plage_output_df) <- Rehman_exp_coldata$Simple_id

plage_heatmap <- pheatmap(plage_output_df,
                          annotation_col = condition,
                          annotation_colors = condition_colors,
                          color = pal,
                          breaks = myBreaks,
                          border_color = "NA")

# png(filename = "./Rehman_IL17_TNFa_plage_heatmap.png", 
#     width = 1800, 
#     height = 1500, 
#     res = 150)
# plage_heatmap
# dev.off()


# Perform t-tests for each gene
plage_t_test <- data.frame(t(apply(X = plage_output_df,
                               MARGIN = 1,
                               FUN = function(x) {
                                 t_test_result <- t.test(x[control_samples], x[inf_samples])
                                 print(class(t_test_result))
                                 return(c(t_test_result$statistic, t_test_result$p.value))
                               })))
colnames(plage_t_test) <- c("statistic", "p.value")

plage_t_test$statistic.abs <- abs(plage_t_test$statistic)
plage_t_test$adj.p.value <- p.adjust(plage_t_test$p.value, 
                                    method = "fdr")

plage_t_test <- plage_t_test[order(plage_t_test$statistic.abs, decreasing = T),]

write.table(plage_t_test[,c("statistic","adj.p.value")],
            file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/project rROMA/rRoma_comp/Rehman_TNF_IL17_plage_t_test_results.tsv",
            sep = "\t",
            row.names = T)


length(which(plage_t_test$adj.p.value<0.05))


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
colnames(gsva_output_df) <- Rehman_exp_coldata$Simple_id

# gsva_heatmap <- pheatmap(gsva_output_df,
#                           annotation_col = condition,
#                           annotation_colors = condition_colors,
#                           color = pal,
#                           breaks = myBreaks,
#                           border_color = "NA")

# png(filename = "./Rehman_IL17_TNFa_gsva_heatmap.png", 
#     width = 1800, 
#     height = 1500, 
#     res = 150)
# gsva_heatmap
# dev.off()

# Perform t-tests for each gene
gsva_t_test <- data.frame(t(apply(X = gsva_output_df,
                                   MARGIN = 1,
                                   FUN = function(x) {
                                     t_test_result <- t.test(x[control_samples], x[inf_samples])
                                     print(class(t_test_result))
                                     return(c(t_test_result$statistic, t_test_result$p.value))
                                   })))
colnames(gsva_t_test) <- c("statistic", "p.value")

gsva_t_test$statistic.abs <- abs(gsva_t_test$statistic)
gsva_t_test$adj.p.value <- p.adjust(gsva_t_test$p.value, 
                                     method = "fdr")

gsva_t_test <- gsva_t_test[order(gsva_t_test$statistic.abs, decreasing = T),]

write.table(gsva_t_test[,c("statistic","adj.p.value")],
            file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/project rROMA/rRoma_comp/Rehman_TNF_IL17_gsva_t_test_results.tsv",
            sep = "\t",
            row.names = T)

length(which(gsva_t_test$adj.p.value<0.05))
which(rownames(gsva_t_test)=="TNFA SIGNALING VIA NFKB")

names(gsva_adjusted_p_values[gsva_adjusted_p_values < 0.05])
table(gsva_adjusted_p_values < 0.05)

gsva_adjusted_p_values[order(gsva_adjusted_p_values, decreasing = F)]

# pheatmap(gsva_output_df)

## SSGSEA

ssgsea_es <- gsva(Rehman_exp_RNA_matrix_norm, 
                  genesets_list, 
                  verbose=FALSE,
                  method="ssgsea")

ssgsea_output_df <- data.frame(ssgsea_es)
colnames(ssgsea_output_df) <- Rehman_exp_coldata$Simple_id

ssgsea_heatmap <- pheatmap(ssgsea_output_df,
                         annotation_col = condition,
                         annotation_colors = condition_colors,
                         color = pal,
                         breaks = myBreaks,
                         border_color = "NA")

# png(filename = "./Rehman_IL17_TNFa_ssgsea_heatmap.png", 
#     width = 1800, 
#     height = 1500, 
#     res = 150)
# ssgsea_heatmap
# dev.off()

# Perform t-tests for each gene
ssgsea_t_test <- data.frame(t(apply(X = ssgsea_output_df,
                                  MARGIN = 1,
                                  FUN = function(x) {
                                    t_test_result <- t.test(x[control_samples], x[inf_samples])
                                    return(c(t_test_result$statistic, t_test_result$p.value))
                                  })))
colnames(ssgsea_t_test) <- c("statistic", "p.value")

ssgsea_t_test$statistic.abs <- abs(ssgsea_t_test$statistic)
ssgsea_t_test$adj.p.value <- p.adjust(ssgsea_t_test$p.value, 
                                    method = "fdr")

ssgsea_t_test <- ssgsea_t_test[order(ssgsea_t_test$statistic.abs, decreasing = T),]

write.table(ssgsea_t_test[,c("statistic","adj.p.value")],
            file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/project rROMA/rRoma_comp/Rehman_TNF_IL17_ssgsea_t_test_results.tsv",
            sep = "\t",
            row.names = T)

length(which(ssgsea_t_test$adj.p.value<0.05))
which(rownames(ssgsea_t_test)=="TNFA SIGNALING VIA NFKB")

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
colnames(zscore_output_df) <- Rehman_exp_coldata$Simple_id

zscore_t_test <- data.frame(t(apply(X = zscore_output_df,
                                    MARGIN = 1,
                                    FUN = function(x) {
                                      t_test_result <- t.test(x[control_samples], x[inf_samples])
                                      return(c(t_test_result$statistic, t_test_result$p.value))
                                    })))
colnames(zscore_t_test) <- c("statistic", "p.value")

zscore_t_test$statistic.abs <- abs(zscore_t_test$statistic)
zscore_t_test$adj.p.value <- p.adjust(zscore_t_test$p.value, 
                                      method = "fdr")

zscore_t_test <- zscore_t_test[order(zscore_t_test$statistic.abs, decreasing = T),]

write.table(zscore_t_test[,c("statistic","adj.p.value")],
            file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/project rROMA/rRoma_comp/Rehman_TNF_IL17_zscore_t_test_results.tsv",
            sep = "\t",
            row.names = T)

length(which(zscore_t_test$adj.p.value<0.05))
which(rownames(zscore_t_test)=="TNFA SIGNALING VIA NFKB")

zscore_pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "PiYG")))(100)
zscore_myBreaks <- seq(min(zscore_output_df), max(zscore_output_df), length.out = 100)

zscore_heatmap <- pheatmap(zscore_output_df,
                           annotation_col = condition,
                           annotation_colors = condition_colors,
                           color = zscore_pal,
                           breaks = zscore_myBreaks,
                           border_color = "NA")

# png(filename = "./Rehman_IL17_TNFa_zscore_heatmap.png", 
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
                                  t_test_result <- t.test(x[control_samples], 
                                                          x[inf_samples])
                                  return(t_test_result$p.value)
                                })

zscore_adjusted_p_values <- p.adjust(zscore_t_test_p_values, 
                                     method = "fdr")



names(zscore_adjusted_p_values[zscore_adjusted_p_values < 0.05])
table(zscore_adjusted_p_values < 0.05)

zscore_adjusted_p_values[order(zscore_adjusted_p_values, decreasing = F)]

# pheatmap(zscore_output_df)

## rROMA

load(file="./Rehman_IL17_TNFa_rRoma_obj.Rdata")

rROMA_stat_outpout_df <- rRoma.output$ModuleMatrix

rownames(rROMA_stat_outpout_df) <- sapply(rownames(rROMA_stat_outpout_df), function(geneset) {
  geneset_name <- gsub(pattern = "HALLMARK_", 
                       replacement = "", 
                       x = geneset)
  geneset_name <- gsub(pattern = "_", 
                       replacement = " ", 
                       x = geneset_name)
  return(geneset_name)
})

write.table(rROMA_stat_outpout_df,
            file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/project rROMA/rRoma_comp/Rehman_TNF_IL17_rROMA_module_matrix.tsv",
            sep = "\t",
            row.names = T)

rROMA_output_df <- rRoma.output$SampleMatrix
colnames(rROMA_output_df) <- gsub(pattern = "_",
                                  replacement = " ",
                                  x = colnames(rROMA_output_df))

rownames(rROMA_output_df) <- sapply(rownames(rROMA_output_df), function(geneset) {
  geneset_name <- gsub(pattern = "HALLMARK_", 
                       replacement = "", 
                       x = geneset)
  geneset_name <- gsub(pattern = "_", 
                       replacement = " ", 
                       x = geneset_name)
  return(geneset_name)
})

rROMA_t_test <- data.frame(t(apply(X = rROMA_output_df,
                                    MARGIN = 1,
                                    FUN = function(x) {
                                      t_test_result <- t.test(x[control_samples], x[inf_samples])
                                      return(c(t_test_result$statistic, t_test_result$p.value))
                                    })))
colnames(rROMA_t_test) <- c("statistic", "p.value")

rROMA_t_test$statistic.abs <- abs(rROMA_t_test$statistic)
rROMA_t_test$adj.p.value <- p.adjust(rROMA_t_test$p.value, 
                                      method = "fdr")

rROMA_t_test <- rROMA_t_test[order(rROMA_t_test$statistic.abs, decreasing = T),]

length(which(rROMA_t_test$adj.p.value<0.05))
which(rownames(rROMA_t_test)=="TNFA SIGNALING VIA NFKB")


rRoma_heatmap <- pheatmap(rROMA_output_df,
                          annotation_col = condition,
                          annotation_colors = condition_colors,
                          color = pal,
                          breaks = myBreaks,
                          border_color = "NA")

# png(filename = "./Rehman_IL17_TNFa_rRoma_heatmap.png", 
#     width = 1800, 
#     height = 1500, 
#     res = 150)
# rRoma_heatmap
# dev.off()

