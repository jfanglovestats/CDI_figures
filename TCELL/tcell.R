library(Matrix)
library(matrixStats)
library(MASS)
library(Rcpp)
library(truncnorm)
library(data.table)
library(mclust)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(ggsci)
library(umap)
library(cluster)
library(BiocParallel)
library(dplyr)




X = readRDS("Tcell_5type_filtered.rds")
cell_info = readRDS("Tcell_5type_filtered_labels.rds")

para.set = list(K = length(unique(cell_info$Group)), 
                ngene = nrow(X),
                celltype = cell_info$Group,
                nfeature = 500,
                ncore = 20)
set.indx = "tcell"
saveRDS(para.set, paste0(set.indx, "_paraset.rds"))


# ------------------------------------------------------------------------------
#                       Generate Candidate Labels
# ------------------------------------------------------------------------------
source("utils/clustering.R")
ncl_vec = c(2,   3,    4,    5,    6,    7,    8,    10,   12,   14,   16)


## Seurat
seu_res = c(0.02, 0.05, 0.1,  0.3,  0.5,  0.6,  0.9,  1.3,  1.75, 2.1,  2.5)
df = data.frame(seu_res, ncluster = seurat_lab(X, seu_res))


## SC3
start_time = Sys.time()
sc3_lab(X, ncl_vec)
end_time = Sys.time()
end_time - start_time

## CIDR
start_time = Sys.time()
cidr_lab(X, ncl_vec)
end_time = Sys.time()
end_time - start_time



X_pc = gcmat_pc(X, npc = 200)
# saveRDS(X_pc, paste0(set.indx, "_bc_pc200_countmat.rds"))
# X_pc = readRDS(paste0(set.indx, "_bc_pc200_countmat.rds"))

## spectral
start_time = Sys.time()
spec_lab_pc(X_pc, ncl_vec)
end_time = Sys.time()
end_time - start_time


## kmeans
start_time = Sys.time()
kmeans_lab_pc(X_pc, ncl_vec)
end_time = Sys.time()
end_time - start_time


## HC
start_time = Sys.time()
hc_lab_pc(X_pc, ncl_vec)
end_time = Sys.time()
end_time - start_time




method_name = c("CIDR", "KMeans", "HC", "SC3", "Seurat", "Spectral")
lab_df = data.frame(tmp = rep(NA, ncol(X)))
for(k in seq_len(length(method_name))){
  for(i in seq_len(length(ncl_vec))){
    cur_df = readRDS(paste0(set.indx, "_", method_name[k], "_k", ncl_vec[i], "_labs.rds"))
    colnames(cur_df) = paste0(method_name[k], "_k", ncl_vec[i])
    lab_df = cbind(lab_df, cur_df)
  }
}
lab_df = lab_df[,-1]
# saveRDS(lab_df, paste0(set.indx, "_clustering_labels.rds"))

# ------------------------------------------------------------------------------
#                           Calculate CDI 
# ------------------------------------------------------------------------------

library(CDI)
## feature selection
selected_feature = feature_gene_selection(
  gcmat = X, 
  Seurat_obj = NULL,
  method = "wds",
  nfeature = 500,
  batch_label = NULL,
  zp_threshold = 0.95)

X_sub = X[selected_feature, ]
saveRDS(X_sub, paste0(set.indx, "_X_sub.rds"))


X_sc = size_factor(X)


cdi_df = calculate_CDI(
  sub_gcmat = X_sub,
  cand_lab_df = lab_df, 
  cell_size_factor = X_sc, 
  batch_label = NULL,
  BPPARAM = MulticoreParam(para.set$ncore))


saveRDS(cdi_df, paste0(set.indx, "_cdi_df.rds"))

benchmark_return = calculate_CDI(
  sub_gcmat = X_sub,
  cand_lab_df = para.set$celltype, 
  cell_size_factor = X_sc, 
  batch_label = NULL,
  BPPARAM = MulticoreParam(para.set$ncore))

saveRDS(benchmark_return, paste0(set.indx, "_benchmark_return.rds"))




# ------------------------------------------------------------------------------
#                            Figures (lineplot & Contingency) 
# ------------------------------------------------------------------------------

pdf(paste0(set.indx, "_bic_lineplot.pdf"), width = 3, height = 3.5)
CDI_lineplot(cdi_dataframe = cdi_df, 
             cdi_type = "CDI_BIC",
             benchmark_maintype_ncluster =length(unique(cell_info$Group)), 
             benchmark_maintype_cdi = lapply(benchmark_return, function(x) x/10^6)) +
  scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
                                rgb(136, 146, 179, max = 255), #light purpleblue -- HC
                                rgb(112, 186, 211, max = 255), #lightblue -- K-means
                                rgb(229, 160, 133, max = 255),#pink -- SC3
                                "#BCBD22B2",#-- Seurat
                                rgb(76, 158, 137, max = 255))) # Spectral
dev.off()




ctname = data.frame(ct = para.set$celltype)
ctname = ctname %>% mutate(ct = recode(ct, 
                                       'Active EM-like Treg' = 'Active EM-like Treg', 
                                       'CD8 Tcm, IL17RA+ & CD28+' = 'CD8 Tcm',
                                       'CD8 Trm cells' = 'CD8 Trm',
                                       'Classical CD4 Tem' = 'Classical CD4 Tem',
                                       'regulatory Trm cells' = 'Regulatory Trm'))

pdf(paste0(set.indx, "_bicopt_contingency.pdf"), width = 4, height = 2.5)
contingency_heatmap(as.character(ctname$ct), 
                    as.character(readRDS(paste0("./clustering_labels/", set.indx, "_SC3_k5_labs.rds"))[,1]))
dev.off()


# ------------------------------------------------------------------------------
#                      UMAPS from selected Features -- Fig.3
# ------------------------------------------------------------------------------

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




## Seurat signature genes
Seurat.obj = CreateSeuratObject(counts = X, project = "sc20a", min.cells = 0, min.features = 0)
Seurat.obj = NormalizeData(Seurat.obj, verbose = FALSE)
Seurat.obj = FindVariableFeatures(Seurat.obj, selection.method = "vst", nfeatures = para.set$nfeature, verbose = FALSE)
seu_sig_gname = Seurat.obj@assays[["RNA"]]@var.features
seu_sig_indx = as.numeric(sub("^.", "", Seurat.obj@assays[["RNA"]]@var.features))

## CDI signature genes
selected_feature = feature_gene_selection(
  gcmat = X, 
  method = "wds",
  nfeature = 500,
  batch_label = NULL,
  zp_threshold = 0.95)

pdf(paste0(dir, set.indx, "_seurat_umap.pdf"), width = 4, height = 4)
print(sig_gene_umap(X[seu_sig_indx, ], para.set$clust.label.real))
dev.off()

pdf(paste0(dir, set.indx, "_cdi_umap.pdf"), width = 4, height = 4)
print(sig_gene_umap(X[selected_feature, ], para.set$clust.label.real))
dev.off()


# ------------------------------------------------------------------------------------
#                  VST-selected feature genes CDI
# ------------------------------------------------------------------------------------

seu_X_sub = X[seu_sig_indx, ]
X_sc = size_factor(X)

seu_cdi_df = calculate_CDI(
  sub_gcmat = seu_X_sub,
  cand_lab_df = lab_df, 
  cell_size_factor = X_sc, 
  batch_label = NULL,
  BPPARAM = MulticoreParam(para.set$ncore))

benchmark_return = calculate_CDI(
  sub_gcmat = seu_X_sub,
  cand_lab_df = para.set$celltype, 
  cell_size_factor = X_sc,
  BPPARAM = MulticoreParam(para.set$ncore))


pdf(paste0(set.indx, "_seu_bic_lineplot.pdf"), width = 3, height = 3.5)
CDI_lineplot(cdi_dataframe = seu_cdi_df, 
             cdi_type = "CDI_BIC",
             benchmark_maintype_ncluster =length(unique(cell_info$Group)), 
             benchmark_maintype_cdi = lapply(benchmark_return, function(x) x/10^6)) +
  scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
                                rgb(136, 146, 179, max = 255), #light purpleblue -- HC
                                rgb(112, 186, 211, max = 255), #lightblue -- K-means
                                rgb(229, 160, 133, max = 255),#pink -- SC3
                                "#BCBD22B2",#-- Seurat
                                rgb(76, 158, 137, max = 255))) # Spectral
dev.off()


# ------------------------------------------------------------------------------------
#        CDI performance on T-CELL using different numbers of feature genes
# ------------------------------------------------------------------------------------

ng_vec = c(200, 300, 400)

i = 1
# i = 2
# i = 3
selected_feature = feature_gene_selection(
  gcmat = X, 
  Seurat_obj = NULL,
  method = "wds",
  nfeature = ng_vec[i],
  batch_label = NULL,
  zp_threshold = 0.95)

X_sub = X[selected_feature, ]
X_sc = size_factor(X)


cdi_df = calculate_CDI(
  sub_gcmat = X_sub,
  cand_lab_df = lab_df, 
  cell_size_factor = X_sc, 
  batch_label = NULL,
  BPPARAM = MulticoreParam(para.set$ncore))

benchmark_return = calculate_CDI(
  sub_gcmat = X_sub,
  cand_lab_df = para.set$celltype, 
  cell_size_factor = X_sc, 
  batch_label = NULL,
  BPPARAM = MulticoreParam(para.set$ncore))

pdf(paste0(set.indx, "_bic_lineplot_",ng_vec[1],"gene.pdf"), width = 3, height = 3.5)
CDI_lineplot(cdi_dataframe = cdi_df, 
             cdi_type = "CDI_BIC",
             benchmark_maintype_ncluster =length(unique(cell_info$Group)), 
             benchmark_maintype_cdi = lapply(benchmark_return, function(x) x/10^6)) +
  scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
                                rgb(136, 146, 179, max = 255), #light purpleblue -- HC
                                rgb(112, 186, 211, max = 255), #lightblue -- K-means
                                rgb(229, 160, 133, max = 255),#pink -- SC3
                                "#BCBD22B2",#-- Seurat
                                rgb(76, 158, 137, max = 255))) # Spectral
dev.off()







# ------------------------------------------------------------------------------------
#                  Heatmaps for six five-cluster candidate label sets
# ------------------------------------------------------------------------------------

cur_k = 5
mname = c("CIDR", "KMeans", "HC", "SC3", "Seurat", "Spectral")
for(met in mname){
  labs = readRDS(paste0("./clustering_labels/", set.indx, "_",met,"_k5_labs.rds"))[,1]
  pdf(paste0(set.indx, "_bicopt_",met,"_k5_contingency.pdf"), width = 5, height = 5)
  contingency_heatmap(as.character(ctname$ct), labs)
  dev.off()
}

