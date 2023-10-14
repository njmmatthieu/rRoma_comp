library("GEOquery")
library(tidyr)
library(tidyverse)
library(ggfortify)
library(ggplot2)
library(edgeR)

# Phenotypes

gse <- getGEO(filename="~/ownCloud/Thèse/Systems Biology/Transcriptomic studies/Transcriptomic Analyses/Analyse Rehman/Data/GSE176121_series_matrix.txt",GSEMatrix=TRUE,getGPL = FALSE)

Rehman_phenotypes <- gse@phenoData@data
Rehman_phenotypes <- Rehman_phenotypes[,c("batch:ch1", 
                                          "disease:ch1", 
                                          "donor:ch1", 
                                          "treatment:ch1")]
colnames(Rehman_phenotypes) <- gsub(pattern = ":ch1", 
                                    replacement = "", 
                                    x = colnames(Rehman_phenotypes))
Rehman_phenotypes$sample_id <- rownames(Rehman_phenotypes)

# Scenario 2 - control vs IL-17+TNFa

exp_samples <- rownames(Rehman_phenotypes[which(Rehman_phenotypes$disease=="non-CF"),])
Rehman_exp_phenotypes <- Rehman_phenotypes[exp_samples, ]

Rehman_exp_coldata <- data.frame(Sample_id = rownames(Rehman_exp_phenotypes),
                            Batch=Rehman_exp_phenotypes$batch, 
                            CF_State=Rehman_exp_phenotypes$disease,
                            Treatment=Rehman_exp_phenotypes$treatment)

Rehman_exp_coldata$Simple_id <- paste(Rehman_exp_coldata$Treatment, 
                                    rep(1:6, 2),
                                    sep = "_")

# # Removing GSM5356221
# Rehman_coldata <- Rehman_coldata[which(Rehman_coldata$Sample_id!="GSM5356221"),]

# Matrix

Rehman_raw_counts <- read.table("~/ownCloud/Thèse/Systems Biology/Transcriptomic studies/Transcriptomic Analyses/Analyse Rehman/Data/GSE176121_raw_counts_aggregated.csv",
                                sep = "\t", 
                                header = T)

colnames(Rehman_raw_counts) <- gsub(pattern = "_Welsh[[:digit:]]+", 
                                    replacement = "", 
                                    colnames(Rehman_raw_counts))

# keep only control samples
rownames(Rehman_raw_counts) <- Rehman_raw_counts$gene_name
Rehman_raw_counts$gene_name <- NULL
Rehman_exp_raw_counts <- Rehman_raw_counts[,colnames(Rehman_raw_counts) %in% exp_samples]

# # Removing GSM5356221
# Rehman_raw_counts$GSM5356221 <- NULL

Rehman_exp_RNA_matrix <- as.matrix(Rehman_exp_raw_counts)

colnames(Rehman_exp_RNA_matrix) <- sapply(colnames(Rehman_exp_RNA_matrix), function(id){
  Treatment_id <- Rehman_exp_coldata %>%
    filter(Sample_id==id) %>%
    pull(Simple_id)
  
  return(Treatment_id)
})

Rehman_exp_RNA_matrix_norm <- cpm(Rehman_exp_RNA_matrix, 
            log=TRUE)

# save(Rehman_exp_RNA_matrix_norm,
#      file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/project rROMA/Reviews/Rehman_exp_RNA_matrix_norm.Rda")

# save(Rehman_exp_coldata,
#      file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/project rROMA/Reviews/Rehman_exp_sample_labels.Rda")

# PCA visualisation

gene.pca <- prcomp(t(lcpm),
                   center = T)

gene.pca.plot <- autoplot(gene.pca,
                          data = t(lcpm),
                          label = T)
