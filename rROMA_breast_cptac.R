# Load libraries

library("DESeq2")
library("limma")
library("FactoMineR")
library("factoextra")
library("ggplot2")
library("pheatmap")
library("rRoma")
library("clusterProfiler")
library("org.Hs.eg.db")
library("ggrepel")
library("dendextend")
library("rRoma")
library("GSEABase")
library("ggdendro")
library("dendextend")
library("missMDA")


# Load dataset
setwd("~/draft_roma_2023/submission_npj/revision/proteomics/breast_cptac")
prot <- read.csv("./prot_breast_cptac.txt", header = TRUE, row.names = 1, sep='\t')
dim(prot)
design <- read.csv("./metadata_breast_cptac.txt", header = TRUE, row.names = 1,sep='\t' )

design

# Convert ProtID to Gene_Name
# Remove duplicated gene names in prot matrix (keep the protein quantification with most peptides)

df <- prot %>% 
  group_by(GeneSymbol) %>% 
  dplyr::slice(which.max(numSpectraProteinObserved)) 

df <- as.data.frame(df)
rownames(df) <- df[,1]
prot <- as.matrix(df[,3:ncol(df)])

# Run PCA, check for outlier samples
pca.prot <- PCA(t(prot), graph = FALSE)
fviz_eig(pca.prot)
fviz_pca_ind(pca.prot, habillage = as.factor(design$PAM50))

# Create complete matrix with imputed data
imputed.prot <- missMDA::imputePCA(prot, ncp = 6)
imputed.prot <- imputed.prot$completeObs
  
# Run hierarchical clustering
d <- dist(t(imputed.prot))

hc <- hclust(d, method = "ward.D")
dend <- as.dendrogram(hc)
plot(dend)

#ggdendrogram(hc)
#colored_bars(colors = design$PAM50, dend = dend)

#plot(hc, hang = -1)
#dend <- as.dendrogram(hc)
#colored_bars(as.numeric(as.factor(design$PAM50)), dend = dend)

# Run rROMA analysis

## Set group labels

Group <- design$PAM50
names(Group) <- rownames(design)
table(Group)

## Selecting the module list

Hallmarks <- SelectFromMSIGdb("HALLMARK")

Hallmarks <- lapply(Hallmarks, function(x){
  x$Name <- sub("HALLMARK_", "", x$Name)
  x
})


## Run rROMA

rRoma.prot <- rRoma.R(imputed.prot, Hallmarks)

### Overdispersed modules
overdispersed.modules <- which(rRoma.prot$ModuleMatrix[, "ppv L1"] <= 0.05 & rRoma.prot$ModuleMatrix[, "ppv Median Exp"] > 0.05)
Plot.Genesets.Samples(rRoma.prot, Selected = overdispersed.modules, GroupInfo = Group, cluster_cols = TRUE)

### Shifted modules

shifted.modules <- which(rRoma.prot$ModuleMatrix[, "ppv Median Exp"] <= 0.05)
Plot.Genesets.Samples(rRoma.prot, Selected = shifted.modules, GroupInfo = Group, cluster_cols = TRUE)

df <- data.frame(Subtype = design$PAM50, row.names = colnames(rRoma.prot$SampleMatrix))

pheatmap(rRoma.prot$SampleMatrix[shifted.modules,], annotation_col = df, show_colnames = FALSE)

# Perform t-tests for each pathway

roma_output_df <- data.frame(rRoma.output$SampleMatrix)

roma_t_test_p_values <- apply(X = roma_output_df,
                              MARGIN = 1,
                              FUN = function(x) {
                                t_test_result <- t.test(x[which(design$Group == "CF")], x[which(design$Group == "NCF")])
                                return(t_test_result$p.value)
                              })

roma_adjusted_p_values <- p.adjust(roma_t_test_p_values, 
                                   method = "fdr")

names(roma_t_test_p_values[roma_t_test_p_values < 0.05])
names(roma_adjusted_p_values[roma_adjusted_p_values < 0.05])
table(roma_t_test_p_values < 0.05)
table(roma_adjusted_p_values < 0.05)

## Run PLAGE

genesets <- GSEABase::getGmt("./h.all.v2023.1.Hs.symbols.gmt")
genesets_list <- lapply(genesets, function(geneset) {
  return(geneset@geneIds)
})
names(genesets_list) <- lapply(genesets, function(geneset) {
  return(geneset@setName)
})

plage_es <- gsva(log.prot, 
                 genesets_list, 
                 verbose=FALSE,
                 method="plage")

plage_output_df <- data.frame(plage_es)

# Perform t-tests for each pathway
plage_t_test_p_values <- apply(X = plage_output_df,
                               MARGIN = 1,
                               FUN = function(x) {
                                 t_test_result <- t.test(x[which(design$Group == "CF")], x[which(design$Group == "NCF")])
                                 return(t_test_result$p.value)
                               })

plage_adjusted_p_values <- p.adjust(plage_t_test_p_values, 
                                    method = "fdr")

names(plage_t_test_p_values[plage_t_test_p_values < 0.05])
names(plage_adjusted_p_values[plage_adjusted_p_values < 0.05])
table(plage_t_test_p_values < 0.05)
table(plage_adjusted_p_values < 0.05)

pheatmap(plage_output_df)

# GSVA

gsva_es <- gsva(log.prot, 
                genesets_list, 
                verbose=FALSE,
                method="gsva",
                kcdf="Poisson")

gsva_output_df <- data.frame(gsva_es)

# Perform t-tests for each gene
gsva_t_test_p_values <- apply(X = gsva_output_df,
                              MARGIN = 1,
                              FUN = function(x) {
                                # print(x[CF_samples])
                                # print(x[non_CF_samples])
                                t_test_result <- t.test(x[which(design$Group == "CF")], x[which(design$Group == "NCF")])
                                return(t_test_result$p.value)
                              })

gsva_adjusted_p_values <- p.adjust(gsva_t_test_p_values, 
                                   method = "fdr")

names(gsva_t_test_p_values[gsva_t_test_p_values < 0.05])
names(gsva_adjusted_p_values[gsva_adjusted_p_values < 0.05])

table(gsva_t_test_p_values < 0.05)
table(gsva_adjusted_p_values < 0.05)

#pheatmap(gsva_output_df)
pheatmap(gsva_output_df)

## SSGSEA

ssgsea_es <- gsva(log.prot, 
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
                                  t_test_result <- t.test(x[which(design$Group == "CF")], x[which(design$Group == "NCF")])
                                  return(t_test_result$p.value)
                                })

ssgsea_adjusted_p_values <- p.adjust(ssgsea_t_test_p_values, 
                                     method = "fdr")

names(ssgsea_t_test_p_values[ssgsea_t_test_p_values < 0.05])
names(ssgsea_adjusted_p_values[ssgsea_adjusted_p_values < 0.05])
table(ssgsea_t_test_p_values < 0.05)
table(ssgsea_adjusted_p_values < 0.05)

pheatmap(ssgsea_output_df)

## Zscore

zscore_es <- gsva(log.prot, 
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
                                  t_test_result <- t.test(x[which(design$Group == "CF")], x[which(design$Group == "NCF")])
                                  return(t_test_result$p.value)
                                })

zscore_adjusted_p_values <- p.adjust(zscore_t_test_p_values, 
                                     method = "fdr")

names(zscore_t_test_p_values[zscore_t_test_p_values < 0.05])
names(zscore_adjusted_p_values[zscore_adjusted_p_values < 0.05])
table(zscore_t_test_p_values < 0.05)
table(zscore_adjusted_p_values < 0.05)

pheatmap(zscore_output_df)
