library(Matrix)
library(matrixStats)
library(MASS)
library(doParallel) 
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

# ------------------------------------------------------------------------------
#                          Simulation setting
# ------------------------------------------------------------------------------


set_mainsub = list(K = 4, 
                   ngene = 10000, 
                   main_sig = 30, 
                   sub_sig = 15,
                   nfeature = 500,
                   log2fc = 4,
                   sub_log2fc = 1.8,
                   ncore = 20)

set_mainsub$ncell = c(1000, 600, 600, 600)
set_mainsub$subtype = rep(1:set_mainsub$K, set_mainsub$ncell)
set_mainsub$celltype = c(rep(1, 1000), rep(2, 1800))
para.set = set_mainsub
set.indx = "sd3"
para.set$N0 = sum(para.set$ncell)
# saveRDS(para.set, paste0(set.indx, "_paraset.rds"))






set.seed(2)

## baseline
generic_mu = 0.4; generic_r = 1
std_dev = 0.1
mu = matrix(rep(rtruncnorm(para.set$ngene, a = generic_mu - 2*std_dev, b =  generic_mu + 2*std_dev, mean = generic_mu, sd = std_dev), para.set$K), 
            nrow=para.set$ngene, byrow = F)
r = matrix(rep(rtruncnorm(para.set$ngene, a = generic_r - 2*std_dev, b = generic_r + 2*std_dev, mean = generic_r, sd = 0.1), para.set$K), 
           nrow=para.set$ngene, byrow = F)



## signature
mu[1:para.set$main_sig, 2:4] = mu[1:para.set$main_sig, 1] *(2^para.set$log2fc)
for(j in 1:3){
  start_indx = para.set$main_sig + (j-1)*para.set$sub_sig + 1
  end_indx = para.set$main_sig + j*para.set$sub_sig
  mu[start_indx:end_indx, j+1] = mu[start_indx:end_indx, j+1] *(2^para.set$sub_log2fc)
  r[start_indx:end_indx, j+1] = r[start_indx:end_indx, j+1]+ rnorm(length(start_indx:end_indx), 0, sd = 0.05) 
}

X = matrix(NA, nrow = para.set$ngene, ncol = para.set$N0)
for(g in 1:para.set$ngene){
  tmp = NULL
  for(cl in 1:para.set$K){
    cur.mu = mu[g,cl]
    if(cur.mu == 0){ tmp = c(tmp, rep(0, para.set$ncell[cl]))}
    else{ tmp = c(tmp, rnbinom(para.set$ncell[cl], mu = cur.mu, size = r[g,cl]))}
  }
  X[g,] = tmp
}
colnames(X) = paste0("c",1:para.set$N0)
rownames(X) = paste0("g",1:para.set$ngene)

# ------------------------------------------------------------------------------
#                       Generate Candidate Labels
# ------------------------------------------------------------------------------

source("utils/clustering.R")
ncl_vec = c(2,       3,    4,       5,      6,  7,  8)


## Seurat
seu_res = c(0.449,  0.467, 1,   1.347,  1.374, 1.549, 1.641)
## The resoultion for a specific number of clusters is found by binary search
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
## saveRDS(X_pc, paste0(set.indx, "_bc_pc200_countmat.rds"))
## X_pc = readRDS(paste0(set.indx, "_bc_pc200_countmat.rds"))

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
#                          Calculate CDI 
# ------------------------------------------------------------------------------

library(CDI)
## feature selection
selected_feature = feature_gene_selection(
  gcmat = X, 
  method = "wds",
  nfeature = 500,
  batch_label = NULL,
  zp_threshold = 0.95)

X_sub = X[selected_feature, ]
# saveRDS(X_sub, paste0(set.indx, "_X_sub.rds"))


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
maintype_return = calculate_CDI(
  sub_gcmat = X_sub, 
  cand_lab_df = para.set$maintype, 
  batch_label = NULL, 
  cell_size_factor = X_sc, 
  BPPARAM = MulticoreParam(para.set$ncore))

celltype_return = calculate_CDI(
  sub_gcmat = X_sub, 
  cand_lab_df = para.set$clust.label.real, 
  batch_label = NULL, 
  cell_size_factor = X_sc, 
  BPPARAM = MulticoreParam(para.set$ncore))
end_time = Sys.time()
difftime(end_time, start_time)

# saveRDS(list(benchmark_main = maintype_return, 
#              benchmark_cell = celltype_return), 
#         paste0(set.indx, "_benchmark_return.rds"))




# ------------------------------------------------------------------------------
#                     Lineplot & Heatmap -- Fig.4
# ------------------------------------------------------------------------------

pdf(paste0(set.indx, "_bic_lineplot.pdf"), width = 4, height = 4)
CDI_lineplot(cdi_dataframe = cdi_df, 
             cdi_type = "CDI_BIC",
             benchmark_maintype_ncluster = length(unique(para.set$celltype)), 
             benchmark_maintype_cdi = lapply(benchmark_return$benchmark_main, function(x) x/10^6)) +
  scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
                                rgb(136, 146, 179, max = 255), #light purpleblue -- HC
                                rgb(112, 186, 211, max = 255), #lightblue -- K-means
                                rgb(229, 160, 133, max = 255),#pink -- SC3
                                "#BCBD22B2",#-- Seurat
                                rgb(76, 158, 137, max = 255))) # Spectral
dev.off()


pdf(paste0(set.indx, "_aic_lineplot.pdf"), width = 4, height = 4)
CDI_lineplot(cdi_dataframe = cdi_df, 
             cdi_type = "CDI_AIC",
             benchmark_celltype_ncluster = length(unique(para.set$subtype)), 
             benchmark_celltype_cdi = lapply(benchmark_return$benchmark_cell, function(x) x/10^6)) +
  scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
                                rgb(136, 146, 179, max = 255), #light purpleblue -- HC
                                rgb(112, 186, 211, max = 255), #lightblue -- K-means
                                rgb(229, 160, 133, max = 255),#pink -- SC3
                                "#BCBD22B2",#-- Seurat
                                rgb(76, 158, 137, max = 255))) # Spectral
dev.off()





pdf(paste0(set.indx, "_bicopt_contingency.pdf"), width = 4, height = 4)
contingency_heatmap(as.character(para.set$celltype), 
                    as.character(readRDS(paste0(set.indx, "_CIDR_k2_labs.rds"))[,1]))
dev.off()

pdf(paste0(set.indx, "_aicopt_contingency.pdf"), width = 4, height = 4)
contingency_heatmap(as.character(para.set$subtype), 
                    as.character(readRDS(paste0(set.indx, "_Seurat_k4_labs.rds"))[,1]))
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
    theme_bw()  + 
    scale_color_manual(values = c(rgb(1, 0.2, 0.7), rgb(0.2, 0.7, 1), rgb(0, 0.6, 0.6), rgb(0, 1, 1))) +
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
print(sig_gene_umap(X[seu_sig_indx, ], para.set$subtype))
dev.off()

pdf(paste0(dir, set.indx, "_cdi_umap.pdf"), width = 4, height = 4)
print(sig_gene_umap2(X[selected_feature, ], para.set$subtype))
dev.off()



