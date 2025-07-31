#Heat map construction
row.names(pdata(my_eset_norm))
library(stringr)
library(pheatmap)
pData(my_eset_norm)
disease_names <- ifelse(str_detect(pData(my_eset_norm)$Factor.Value..DISEASE.STATE. ,
                                   "normal"), "normal" , "lung cancer")
annotation_for_heatmap <- data.frame(disease=disease_names)
annotation_for_heatmap
row.names(annotation_for_heatmap) <- row.names(pData(my_eset_norm))
row.names(pData(my_eset_norm))
row.names(annotation_for_heatmap)
dists <- as.matrix(dist(t(my_eset_norm),method = "manhattan"))
dists
rownames(dists) <- row.names(pData(my_eset_norm))
rownames(dists)
#till there was the part for data collection for heatmap

hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9,"YlorRd"))(255))
colnames(dists) <- NULL
diag(dists) <- NA

ann_colors <- list(
  disease = c(normal= "green" , lung cancer = "red"))
pheatmap(dists,col=(hmcol),
         annotation_row = annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE,
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = TRUE),
                           max(dists, na.rm = TRUE)),
         legend_labels = (c("small distance) , 
)
  
)
