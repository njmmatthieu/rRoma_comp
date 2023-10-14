library(rRoma)

load("./Rehman_IL17_TNFa_RNA_exp_matrix_norm.Rda")
load("./Rehman_IL17_TNFa_sample_labels.Rda")

Treatment <- Rehman_exp_coldata$Treatment
names(Treatment) <- Rehman_exp_coldata$Simple_id

table(Treatment)

# Selecting the modules

AllHall <- SelectFromMSIGdb("HALLMARK")

AllHall <- lapply(AllHall, function(x){
  x$Name <- sub("HALLMARK_", "", x$Name)
  x
})

set.seed(123)
rRoma.output <- rRoma.R(Rehman_exp_RNA_matrix_norm, 
                        AllHall)

save(rRoma.output,
     file = "./Rehman_IL17_TNFa_rRoma_obj.Rdata")

shifted.modules <- which(rRoma.output$ModuleMatrix[, "ppv Median Exp"] <= 0.05)

rROMA_samples_mat_shifted <- rRoma.output$SampleMatrix[shifted.modules,]
pheatmap(rROMA_samples_mat_shifted,
         )

# RColorBrewer::

Plot.Genesets.Samples(rRoma.output, 
                      Selected = shifted.modules, 
                      ColorGradient = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(50),
                      # Fixed = 50,
                      GroupInfo = Treatment, 
                      cluster_cols = TRUE)




overdispersed.modules <- which(rRoma.output$ModuleMatrix[, "ppv L1"] <= 0.05 & rRoma.output$ModuleMatrix[, "ppv Median Exp"] > 0.05)

Plot.Genesets.Samples(rRoma.output, 
                      Selected = overdispersed.modules, 
                      ColorGradient = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(50),,
                      # Fixed = 50,
                      GroupInfo = Treatment, 
                      cluster_cols = TRUE)
