# ------------------------------------------------------------------------------
#                            Libraries 
# ------------------------------------------------------------------------------

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(dplyr)
library(matrixStats)
library(data.table)
library(mclust)
library(ggsci)
library(cluster)
library(BiocParallel)


# SeuratDisk in installed via conda
# virtual env name "r4"

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)


# ------------------------------------------------------------------------------
#                   Read data and check batch effect
# ------------------------------------------------------------------------------

su <- LoadH5Seurat("su_2020_processed.HDF5")
head(su@meta.data, 2)

Idents(su) <- su@meta.data["predicted.celltype.l1"]
pdf("su_maintype_umap.pdf", width = 7, height = 6)
DimPlot(su, reduction = "ref.umap") + labs(x = "UMAP_1", y = "UMAP_2")
dev.off()



Idents(su) <- su@meta.data["predicted.celltype.l2"]
pdf("su_subtype_umap.pdf", width = 10, height = 6)
DimPlot(su, reduction = "ref.umap") + labs(x = "UMAP_1", y = "UMAP_2")
dev.off()



Start = Sys.time()
stephen <- LoadH5Seurat("stephenson_2021_processed.HDF5")
End = Sys.time()
End - Start





head(stephen@meta.data, 2)
Idents(stephen) <- stephen@meta.data["predicted.celltype.l1"]
pdf("stephen_maintype_umap.pdf", width = 7, height = 6)
DimPlot(stephen, reduction = "ref.umap") + labs(x = "UMAP_1", y = "UMAP_2")
dev.off()



Idents(stephen) <- stephen@meta.data["predicted.celltype.l2"]
pdf("stephen_subtype_umap.pdf", width = 10, height = 6)
DimPlot(stephen, reduction = "ref.umap") + labs(x = "UMAP_1", y = "UMAP_2")
dev.off()



# ------------------------------------------------------------------------------
#                   Merge by hand
# ------------------------------------------------------------------------------

su_feature <- rownames(su@assays$RNA@counts)
ste_feature <- rownames(stephen@assays$RNA@counts)
length(intersect(su_feature, ste_feature))

# > length(su_feature)
# [1] 33538
# > length(ste_feature)
# [1] 23491
# > sum(ste_feature %in% intersect(su_feature, ste_feature))
# [1] 23491

su_sub <- subset(su, features = ste_feature)
# su_sub_feature <- rownames(su_sub@assays$RNA@counts)
# sum(su_sub_feature != ste_feature)
# [1] 0
# SaveH5Seurat(su_sub, file.name = "su_sub_2020_processed.HDF5")
feature_order = readRDS("covid_gene_names.rds")
saveRDS(stephen@assays$SCT@counts[feature_order, ], "stephen_count.rds")


com_meta_data <- data.frame(
  dataset = c(rep("su", ncol(su@assays$RNA@counts)), rep("stephen", ncol(stephen@assays$RNA@counts))),
  predicted.celltype.l1 = c(su@meta.data$predicted.celltype.l1, stephen@meta.data$predicted.celltype.l1),
  predicted.celltype.l2 = c(su@meta.data$predicted.celltype.l2, stephen@meta.data$predicted.celltype.l2),
  predicted.celltype.l3 = c(su@meta.data$predicted.celltype.l3, stephen@meta.data$predicted.celltype.l3),
  disease_status_standard = c(su@meta.data$disease_status_standard, stephen@meta.data$disease_status_standard))
# > sum(colnames(su@assays$RNA@counts) %in% colnames(stephen@assays$RNA@counts))
# [1] 0

rownames(com_meta_data) <- c(colnames(su_sub@assays$RNA@counts), colnames(stephen@assays$RNA@counts))
seu_ste <- CreateSeuratObject(
  counts = cbind(su_sub@assays$RNA@counts, stephen@assays$RNA@counts),
  assay = "RNA",
  meta.data = com_meta_data)
seu_ste[["REFUMAP"]] <- CreateDimReducObject(embeddings = rbind(su_sub@reductions$ref.umap@cell.embeddings, 
                                                                stephen@reductions$ref.umap@cell.embeddings), 
                                             key = "REFUMAP_", 
                                             assay = DefaultAssay(seu_ste))
seu_ste[["REFPC"]] <- CreateDimReducObject(embeddings = rbind(su_sub@reductions$ref.pca@cell.embeddings, 
                                                              stephen@reductions$ref.pca@cell.embeddings), 
                                           key = "REFPC_", 
                                           assay = DefaultAssay(seu_ste))


SaveH5Seurat(seu_ste, file.name = "su_ste.HDF5")





saveRDS(seu_ste@meta.data, "covid_cellinfo.rds")
saveRDS(su_sub@assays$RNA@counts, "su_count.rds")
saveRDS(stephen@assays$RNA@counts, "stephen_count.rds")




# ----------------------------------------------
#            Write MM files
# ----------------------------------------------
setwd("/hpc/group/xielab/jf243/covid_immune")
su_count = readRDS("su_count.rds")
ste_count = readRDS("stephen_count.rds")
ng = nrow(su_count)
su_nc = ncol(su_count); ste_nc = ncol(ste_count); nc = as.integer(su_nc + ste_nc)
nc_size = 2000; nc_file = as.integer((su_nc + ste_nc) %/% nc_size + 1)
ng_size = 200; ng_file = as.integer(ng %/% ng_size + 1)
com_count = cbind(su_count, ste_count)

# split by genes 
for(i in seq_len(ng_file - 1)){
  if(i %% 5 == 0){cat("i=", i,"\n")}
  writeMM(com_count[((i-1)*ng_size+1):(i*ng_size), ], paste0("./covid_split_gene/covid_gene_p", i, ".txt"))
}
writeMM(com_count[((ng_file-1)*ng_size):ng, ], paste0("./covid_split_gene/covid_gene_p", ng_file, ".txt"))

# split by cells
for(j in seq_len(nc_file - 1)){
  writeMM(com_count[,((j-1)*nc_size+1):(j*nc_size)], paste0("./covid_split_cell/covid_cell_p", j, ".txt"))
}
writeMM(com_count[,((nc_file-1)*nc_size):nc], paste0("./covid_split_cell/covid_cell_p", nc_file, ".txt"))



# ----------------------------------------------
#            Size factor
# ----------------------------------------------

gene_ref = c()
for(i in seq_len(ng_file)){
  X = readMM(paste0("./covid_split_gene/covid_gene_p", i, ".txt"))
  gene_ref = c(gene_ref, exp(rowMeans(log(pmax(as.matrix(X), 0.5)))))
}
names(gene_ref) = rownames(com_count)
saveRDS(gene_ref, "covid_gene_ref_size.rds")

sc_vec = c()
for(j in seq_len(nc_file)){
  if(j %% 5 == 0){cat("j = ", j, "\n")}
  X = readMM(paste0("./covid_split_cell/covid_cell_p", j, ".txt"))
  sc_vec = c(sc_vec, colMedians(sweep(pmax(as.matrix(X), 0.5), 1, gene_ref, "/")))
}
saveRDS(sc_vec, "covid_size_factor.rds")


# ----------------------------------------------
#            Feature Gene Selection
# ----------------------------------------------

su_nc = 559517; ste_nc = 691683
feature_df = data.frame()
for(i in seq_len(ng_file)){
  if(i %% 5 == 0){cat("i = ", i, "\n")}
  X = as.matrix(readMM(paste0("./covid_split_gene/covid_gene_p", i, ".txt")))
  X1 = X[, c(1:su_nc)]; X2 = X[, (su_nc+1):(su_nc + ste_nc)]
  mu1 = rowMeans(X1); mu2 = rowMeans(X2)
  var1 = rowVars(X1); var2 = rowVars(X2)
  gene_zp1 = rowMeans(X1 == 0); gene_zp2 = rowMeans(X2 == 0)
  gene_phi1 = ifelse(gene_zp1 > 0.95, NA, (var1 - mu1)/mu1^2)
  gene_phi2 = ifelse(gene_zp2 > 0.95, NA, (var2 - mu2)/mu2^2)
  tmp_feature = data.frame(zp1 = gene_zp1, phi1 = gene_phi1, zp2 = gene_zp2, phi2 = gene_phi2)
  feature_df = rbind(feature_df, tmp_feature)
}







tmp_phi1 = ifelse(is.na(feature_df$phi1), -10000, feature_df$phi1)
tmp_phi2 = ifelse(is.na(feature_df$phi2), -10000, feature_df$phi2)
phi_rank1 = rank(-tmp_phi1); phi_rank2 = rank(-tmp_phi2)
gene_min_rank <- pmin(phi_rank1, phi_rank2)
feature_df["rank1"] = phi_rank1
feature_df["rank2"] = phi_rank2
feature_df["rank"] = gene_min_rank
saveRDS(feature_df, "covid_feature_df.rds")




feature_df = readRDS("covid_feature_df.rds")
feature_gene_indx = c(980, which(feature_df$rank < 259))
# indx 980, 3144 -> 259


X_sub = com_count[feature_gene_indx, ]
writeMM(X_sub, "covid_X_sub.txt")



# ----------------------------------------------
#            Generate labels
# ----------------------------------------------


# setwd("/hpc/group/xielab/jf243/covid_immune")
seu_obj <- LoadH5Seurat("su_ste.h5Seurat")


su <- LoadH5Seurat("su_2020_processed.HDF5")
ste = LoadH5Seurat("stephenson_2021_processed.HDF5")

all_su_spca = su@reductions$ref.spca@cell.embeddings
all_ste_spca = ste@reductions$ref.spca@cell.embeddings

seu_obj[["REFSPC"]] <- CreateDimReducObject(embeddings = rbind(all_su_spca, all_ste_spca), 
                                            key = "REFSPC_", 
                                            assay = DefaultAssay(seu_obj))

dim(seu_obj@reductions$REFSPC)


SaveH5Seurat(seu_obj, file.name = "su_ste.h5Seurat")

Start = Sys.time()
seu_obj <- FindNeighbors(seu_obj, dims = 1:10, reduction = "REFSPC")
End = Sys.time()
End - Start
# Time difference of 9.507046 mins

BS_resolution2 = function(ideal_K, Seurat.obj, res_seq = seq(0.1, 2, by = 0.1)){
  nres = length(res_seq)
  res_indx = round(nres / 2)
  max_iter = 20
  iter = 1
  while((iter <= max_iter) & (res_indx >= 1) & (res_indx <= nres)){
    cur_res = res_seq[res_indx]
    Seurat.obj = FindClusters(Seurat.obj, resolution = cur_res, verbose = F)
    cur_ncl = length(unique(Seurat.obj@meta.data$seurat_clusters))
    if(cur_ncl == ideal_K){
      return(cur_res)
    } else if(cur_ncl > ideal_K){
      cat("ncl = ", cur_ncl, "res=", cur_res,"\n")
      res_indx = round(res_indx / 2)
    } else if(cur_ncl < ideal_K){
      cat("ncl = ", cur_ncl, "res=", cur_res,"\n")
      res_indx = round(res_indx * 1.5)
    }
    iter = iter + 1
  }
  if(iter == max_iter){
    cat("last resolution = ", cur_res, ". Exceed maximum iteration!")
  }
  return("Not found!")
}


cur_res = c(0.1, seq(0.5, 8, by = 0.5))
for(res in cur_res){
  seu_obj <- FindClusters(seu_obj, resolution = res, verbose = FALSE)
  nk <- length(unique(seu_obj@meta.data$seurat_clusters))
  df <- data.frame(Seurat = seu_obj@meta.data$seurat_clusters)
  saveRDS(df, paste0("covid_Seurat_k", nk, "_labs.rds"))
  cat("ncl = ", nk, "res=", res,"\n")
}



ncl_vec = c(12, 21, 35, 47, 55, 66, 73, 83, 94, 101, 116, 120, 126, 137, 145, 157, 162)
method_name = c("Seurat")
lab_df = readRDS(paste0("covid_Seurat_k12_labs.rds"))
for(k in seq_len(length(method_name))){
  for(i in seq_len(length(ncl_vec))){
    cur_df = readRDS(paste0("covid_", method_name[k], "_k", ncl_vec[i], "_labs.rds"))
    colnames(cur_df) = paste0(method_name[k], "_k", ncl_vec[i])
    lab_df = cbind(lab_df, cur_df)
  }
}
lab_df = lab_df[,-1]
saveRDS(lab_df, "covid_lab_df.rds")







# ----------------------------------------------
#            Calculate CDI (put in cluster)
# ----------------------------------------------

# code: covid_cluster_l3.R, covid_cluster_s1.R
# suggest to calculate 1 label for each bash script to save some space

# set.indx = "covid"
# info = readRDS("covid_cellinfo.rds")
# X_sub = as.matrix(readMM("covid_X_sub.txt"))
# X_sc = readRDS("covid_size_factor.rds")



# ------------------------------------------------------------------------------
#                             Line plot 
# ------------------------------------------------------------------------------


library(ggplot2)
library(ggsci)
library(CDI)
# setwd("/Users/jiyuanfang/Documents/CDI_review/results_for_review/cdi_covid/")
# l1 = readRDS("cdi_covid_benchmark_l1.rds")
# l2 = readRDS("cdi_covid_benchmark_l2.rds")
l3 = readRDS("cdi_covid_benchmark_l3.rds")

ncl_vec = c(12, 21, 35, 47, 55, 66, 73, 83, 94, 101, 116, 120, 126, 137, 145, 157, 162)
cdi_df = data.frame(N_cluster = ncl_vec, CDI_AIC = numeric(length(ncl_vec)), 
                    CDI_BIC = numeric(length(ncl_vec)), Cluster_method = rep("Seurat", length(ncl_vec)))
for(i in seq_len(length(ncl_vec))){
  cdi_ret = readRDS(paste0("cdi_covid_seurat_k", ncl_vec[i], ".rds"))
  cdi_df[i, "N_cluster"] = ncl_vec[i]
  cdi_df[i, "CDI_AIC"] = cdi_ret$CDI_AIC / 10^6
  cdi_df[i, "CDI_BIC"] = cdi_ret$CDI_BIC / 10^6
}

# setwd("/Users/jiyuanfang/Documents/CDI_review/results_for_review/figs/")

pdf("covid_bic_lineplot.pdf", width = 3, height = 3.5)
CDI_lineplot(cdi_df, "CDI_BIC", benchmark_maintype_cdi = lapply(l2, "/", 10^6), 
             benchmark_maintype_ncluster = 31, show_axis_names = FALSE, show_method_legend = FALSE) +
  scale_color_manual(values = c("#BCBD22B2")) 
dev.off()

# ------------------------------------------------------------------------------
#                           heatmap of reference labels 
# ------------------------------------------------------------------------------
setwd("/Users/jiyuanfang/Documents/CDI_review/results_for_review/")
tab = readRDS("covid_cell_label_table.rds")



longData <- reshape2::melt(tab)
longData <- longData[longData$value != 0, ]

pdf("covid_lab_table.pdf", height =  10, width = 8)
ggplot(longData, aes(x = factor(Var1), y = factor(Var2), fill = value)) + 
  geom_tile() + scale_fill_distiller(palette = "Spectral", 
                                     values = rescale(c(1000, 10000, 20000, 50000, 100000, 300000)), 
                                     breaks = c(1000, 10000, 20000, 50000, 100000, 300000)) + 
  theme_classic() + 
  labs(x = "Layer 2 annotation", y = "Layer 3 annotation", fill = "ncells") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=10), 
        axis.text.y = element_text(size=10, angle=0, vjust=0.3), 
        axis.line = element_line(size = 1),
        legend.key.height = unit(2.5, "cm")) 
dev.off()





## order benchmark cell types according to the composition of candidate cell types
order_level <- NULL
for(cell_type in levels(as.factor(longData$benchmark_label))){
  large_prop_indx <- names(which(prop_mat[,cell_type] >= 0.2))
  large_prop_sort <- large_prop_indx[order(prop_mat[large_prop_indx, cell_type], decreasing = TRUE)]
  distinct_indx <- setdiff(large_prop_sort, order_level)
  order_level <- c(order_level, distinct_indx)
}
p1 <- ggplot(longData, 
             aes(x = factor(candidate_label, levels = order_level), y = factor(benchmark_label), fill = PropValue)) + 
  scale_fill_gradientn(limits = c(0,1), 
                       colours = c(rgb(204,204,204, maxColorValue = 255), rgb(25,150,125, maxColorValue = 255))) + 
  labs(x = "Candidate label", y = "Benchmark label", fill = "Proportion") + 
  theme_classic() + geom_tile() + 
  theme(axis.text = element_text(size=10, angle=0, vjust=0.3), axis.line = element_line(size = 1))

