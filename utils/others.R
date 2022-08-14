library(umap)
library(Seurat)
library(ggplot2)
library(ggsci)

sig_gene_umap = function(gcmat, col_group){
  X_pseudo = gcmat
  X_pseudo[X_pseudo == 0] = 0.1
  X_normalize = log(X_pseudo)
  X_pc = prcomp(t(X_normalize), scale = TRUE)$x[,1:50]
  umap_df = umap(d = X_pc, random_state = 123, min_dist = 0.3,  spread = 1, n_neighbors = 10)
  df <- data.frame(x = umap_df$layout[,1],
                   cluster = factor(col_group, levels = c(1:10)),
                   y = umap_df$layout[,2])
  p1 = ggplot(df, aes(x, y, colour = cluster)) + geom_point(size = 1, alpha = 0.5) +  
    theme_bw()  + scale_color_d3("category10") + 
    theme(axis.title=element_blank(),
          legend.title = element_text(size=15), legend.text = element_text(size=15),
          axis.text = element_blank())
  p2 = p1 + guides(color = guide_legend(override.aes = list(size = 5)))
  return(p2)
}