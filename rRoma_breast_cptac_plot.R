# Load rRoma output
load(file="./rRoma_obj_breast_cptac.Rdata")

rRoma.prot.sample_mat <- rRoma.prot$SampleMatrix

rownames(rRoma.prot.sample_mat) <- lapply(rownames(rRoma.prot.sample_mat), function(geneset) {
  geneset_name <- gsub(pattern = "_", 
                       replacement = " ", 
                       x = geneset)
  return(geneset_name)
})

df <- data.frame(Subtype = design$PAM50, row.names = colnames(rRoma.prot$SampleMatrix))
# List with colors for each annotation.
df_colors <- list(Subtype = brewer.pal(5, "Set2"))
names(df_colors$Subtype) <- unique(df$Subtype)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100)
myBreaks <- seq(-0.4, 0.4, length.out = 100)

### Overdispersed modules
overdispersed.modules <- which(rRoma.prot$ModuleMatrix[, "ppv L1"] <= 0.05 & rRoma.prot$ModuleMatrix[, "ppv Median Exp"] > 0.05)

overdispersed_heatmap <- pheatmap(rRoma.prot.sample_mat[overdispersed.modules,], 
                          annotation_col = df,
                          annotation_colors = df_colors,
                          color = pal,
                          breaks = myBreaks,
                          border_color = "NA",
                          show_colnames = F)

png(filename = "./rRoma_breast_cptac_activities_overdispersed_pathways_heatmap.png", 
    width = 1800, 
    height = 600, 
    res = 150)
overdispersed_heatmap
dev.off()

### Shifted modules

shifted.modules <- which(rRoma.prot$ModuleMatrix[, "ppv Median Exp"] <= 0.05)

shifted_heatmap <- pheatmap(rRoma.prot.sample_mat[shifted.modules,], 
                                  annotation_col = df,
                                  annotation_colors = df_colors,
                                  color = pal,
                                  breaks = myBreaks,
                                  border_color = "NA",
                                  show_colnames = F)

png(filename = "./rRoma_breast_cptac_activities_shifted_pathways_heatmap.png", 
    width = 1800, 
    height = 1500, 
    res = 150)
shifted_heatmap
dev.off()
