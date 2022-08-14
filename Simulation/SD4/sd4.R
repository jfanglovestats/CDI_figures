library(Matrix)
library(matrixStats)
library(MASS)
library(doParallel) 
library(Rcpp)
library(splatter)
library(data.table)
library(mclust)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(ggsci)
library(umap)
library(cluster)
library(BiocParallel)

SelectCells = function(scdata, nzg_prop){
  nzg_of_each_cell = colSums(scdata>0)
  ng = nrow(scdata)
  select_cell = c(1:ncol(scdata))[(nzg_of_each_cell/ng > nzg_prop)]
  return(select_cell)
}
SelectGenes = function(scdata, nzc_prop){
  nzc_of_each_gene = rowSums(scdata>0)
  nc = ncol(scdata)
  select_gene = c(1:nrow(scdata))[(nzc_of_each_gene/nc > nzc_prop)]
  return(select_gene)
}

# ------------------------------------------------------------------------------
#                          Simulation setting
# ------------------------------------------------------------------------------


rm(list = setdiff(ls(), lsf.str()))
set.indx = "sd4"

ng = 5000; nc = 3000; ncelltype = 5; nfeature = 500
params = newSplatParams(batchCells = nc, nGenes = ng, lib.scale = 0.1, bcv.common = 0.1, 
                        group.prob = rep(0.2, ncelltype), de.prob = rep(0.01, ncelltype),
                        de.facLoc = 0.4, de.facScale = 0.1, de.downProb = 0, seed = 1)

splat_object = splatSimulateGroups(params,
                                   verbose = FALSE,
                                   dropout.type = "none")


sig_gene_indicator = matrix(NA, nrow = ng, ncol = ncelltype)
colnames(sig_gene_indicator) = paste0("Group", 1:ncelltype)
for(j in 1:ncelltype){
  group_name = paste0("DEFacGroup", j)
  sig_gene_indicator[,j] = (rowData(splat_object)[group_name][,1] > 1)
}

splat_count = counts(splat_object)
process_cell_indx = SelectCells(splat_count, 0.01)
cell_labels = colData(splat_object)$Group[process_cell_indx]
splat_count_process_cell = splat_count[,process_cell_indx]

## remove genes with more than 99% 0
process_gene_indx = SelectGenes(splat_count_process_cell, 0.01)
X = splat_count_process_cell[process_gene_indx, ]
# > dim(X)
# [1] 4887 3000


## signature genes in splat_count 
sig_gene_indx = c(1:nrow(splat_count))[which(rowSums(sig_gene_indicator) >= 1)]
## signature gene name in X (signature genes, after process)
X_sig_gene_name = rowData(splat_object)$Gene[intersect(sig_gene_indx, process_gene_indx)]
X_sig_gene_indx = c(1:nrow(X))[rownames(X) %in% X_sig_gene_name]



para.set = list(K = 5, 
                ngene = nrow(X), 
                nsig = 50,
                ncell = ncol(X),
                total_sig = length(X_sig_gene_indx),
                nfeature = nfeature,
                ncore = 25,
                clust.label.real = cell_labels,
                process.sig.gene = X_sig_gene_indx)

# saveRDS(para.set, paste0(set.indx, "_paraset.rds"))



# ------------------------------------------------------------------------------
#                       Generate Candidate Labels
# ------------------------------------------------------------------------------


source("utils/clustering.R")
ncl_vec = c(2,       3,         4,      5,       6,      7,    8)

## Seurat
seu_res = c(0.0080,  0.0135,  0.0155, 2.0000, 2.4015, 2.5310,  2.847)
seurat_lab(X, seu_res)


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


## spectral
start_time = Sys.time()
spec_lab_pc(X_pc, ncl_vec)
end_time = Sys.time()
end_time - start_time





method_name = c("CIDR", "KMeans", "HC", "SC3", "Seurat", "Spectral")
lab_df = data.frame(tmp = rep(NA, sum(para.set$ncell)))
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
#                           Feature Selection 
# ------------------------------------------------------------------------------

library(CDI)


## via Seurat - VST
Seurat.obj = CreateSeuratObject(counts = X, project = "hh", min.cells = 0, min.features = 0)
Seurat.obj = NormalizeData(Seurat.obj, verbose = FALSE)
Seurat.obj = FindVariableFeatures(Seurat.obj, selection.method = "vst", nfeatures = para.set$nfeature, verbose = FALSE)
seu_sig_gname = Seurat.obj@assays[["RNA"]]@var.features
# seu_sig_indx = as.numeric(sub("^Gene", "", Seurat.obj@assays[["RNA"]]@var.features))
cat("vst = ", sum(seu_sig_gname %in% X_sig_gene_name),"\n")
length(X_sig_gene_name)

## via CDI - WDS

# selected_feature = feature_gene_selection(
#   gcmat = X, 
#   method = "wds",
#   nfeature = 500,
#   batch_label = NULL,
#   zp_threshold = 0.95)

X_sub = X[seu_sig_gname, ]
# saveRDS(X_sub, paste0(set.indx, "_X_sub.rds"))



# ------------------------------------------------------------------------------
#                          Calculate CDI 
# ------------------------------------------------------------------------------


## size factor
X_sc = size_factor(X)
# saveRDS(aicbic_df, paste0(set.indx, "_cdi_df.rds"))

## calculate CDI
aicbic_df = calculate_CDI(
  sub_gcmat = X_sub,
  cand_lab_df = lab_df, 
  cell_size_factor = X_sc, 
  batch_label = NULL,
  BPPARAM = MulticoreParam(para.set$ncore))

# saveRDS(aicbic_df, paste0(set.indx, "_cdi_df.rds"))

cdi_df = aicbic_df %>% 
  mutate(AIC = AIC/(1e+06), BIC = BIC/(1e+06))


start_time = Sys.time()	
benchmark_return = calculate_CDI(
  sub_gcmat = X_sub, 
  cand_lab_df = para.set$clust.label.real, 
  batch_label = NULL, 
  cell_size_factor = X_sc, 
  BPPARAM = MulticoreParam(para.set$ncore))
end_time = Sys.time()
difftime(end_time, start_time)
# saveRDS(benchmark_return, paste0(set.indx, "_benchmark_return.rds"))



# ------------------------------------------------------------------------------
#                     Lineplot & Heatmap -- Fig.4
# ------------------------------------------------------------------------------

pdf(paste0(set.indx, "_bic_lineplot.pdf"), width = 3, height = 3.5)
CDI_lineplot(cdi_dataframe = cdi_df, 
             cdi_type = "CDI_BIC",
             benchmark_maintype_ncluster = length(unique(para.set$clust.label.real)), 
             benchmark_maintype_cdi = lapply(benchmark_return, function(x) x/10^6)) +
  scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
                                rgb(136, 146, 179, max = 255), #light purpleblue -- HC
                                rgb(112, 186, 211, max = 255), #lightblue -- K-means
                                rgb(229, 160, 133, max = 255),#pink -- SC3
                                "#BCBD22B2",#-- Seurat
                                rgb(76, 158, 137, max = 255))) # Spectral
dev.off()


pdf(paste0(set.indx, "_bicopt_contingency.pdf"), width = 2.5, height = 2.5)
contingency_heatmap(as.character(para.set$clust.label.real), 
                    as.character(readRDS(paste0(set.indx, "_Seurat_k5_labs.rds"))[,1]))

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





