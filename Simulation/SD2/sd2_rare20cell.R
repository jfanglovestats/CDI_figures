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
set_rare = list(K = 4, 
                ngene = 10000, 
                sig_indx = list(c(1:40), c(41:80), c(81:130), c(131:180)), 
                nfeature = 500,
                log2fc = c(1.5, 1.5, 2.8, 3.2),
                total_sig = 180,
                ncore = 20)

set_rare$ncell = c(2000, 2000, 20, 100)
set_rare$clust.label.real = rep(1:set_rare$K, set_rare$ncell)
para.set = set_rare
set.indx = "sd2_rare20cell"
para.set$N0 = sum(para.set$ncell)
# saveRDS(para.set, paste0(set.indx, "_paraset.rds"))

set.seed(1000)

## baseline
generic_mu = 0.2
mu = matrix(rep(round(rtruncnorm(para.set$ngene, a=0.05, b=1, mean = generic_mu, sd=0.1), 4), para.set$K), 
            nrow=para.set$ngene)
r = matrix(rep(round(rtruncnorm(para.set$ngene, a =0.2, b=5, mean = 0.5, sd=0.1), 4), para.set$K), 
           nrow=para.set$ngene)



## signature
ls.sig = para.set$sig_indx


for(j in 1:para.set$K){
  sig_g_indx = ls.sig[[j]]
  non_sig_celltype = setdiff(1:para.set$K, j)
  zero_col = sample(non_sig_celltype, round(length(non_sig_celltype)/2))
  zero_g_indx = sample(ls.sig[[j]], round(length(ls.sig[[j]])/3))
  mu[zero_g_indx, zero_col] = 0
}

for(j in 1:para.set$K){
  mu[ls.sig[[j]],j] = pmax(mu[ls.sig[[j]],j] *(2^para.set$log2fc[j]), (generic_mu/1.5) * (2^para.set$log2fc[j])) 
  r[ls.sig[[j]],j] = r[ls.sig[[j]],j] + rnorm(length(ls.sig[[j]]), 0, sd = 0.05) 
}


X = matrix(NA, nrow = para.set$ngene, ncol = para.set$N0)
set.seed(111135)
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
ncl_vec = c(2,       3,       4,        5,       6,       7,     8)

## Seurat
seu_res = c(0.0125,    0.7495,    1.067,    1.124,    1.2343,   1.3711,   1.3437)
hh = data.frame(res = seu_res, ncl = seurat_lab(X, seu_res))

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

## kmeans
start_time = Sys.time()
kmeans_lab_pc(X_pc, ncl_vec)
end_time = Sys.time()
end_time - start_time

## spectral
start_time = Sys.time()
spec_lab_pc(X_pc, ncl_vec)
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
saveRDS(lab_df, paste0(set.indx, "_clustering_labels.rds"))


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


pdf(paste0(set.indx, "_bicaicopt_contingency.pdf"), width = 2.5, height = 2.5)
contingency_heatmap(as.character(para.set$clust.label.real), 
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
    theme_bw()  + scale_color_d3("category10") + 
    theme(axis.title=element_blank(),
          legend.title = element_text(size=15), legend.text = element_text(size=15),
          axis.text = element_blank())
  p2 = p1 + guides(color = guide_legend(override.aes = list(size = 5)))
  return(p2)
}




## CDI signature genes
selected_feature = feature_gene_selection(
  gcmat = X, 
  method = "wds",
  nfeature = 500,
  batch_label = NULL,
  zp_threshold = 0.95)



pdf(paste0(dir, set.indx, "_cdi_umap_rare85cell.pdf"), width = 4, height = 4)
print(sig_gene_umap(X[selected_feature, ], para.set$clust.label.real))
dev.off()





