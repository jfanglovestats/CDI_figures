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

# ------------------------------------------------------------------------------
#                          Simulation setting
# ------------------------------------------------------------------------------


rm(list = setdiff(ls(), lsf.str()))
## parameters
set_10cl = list(K = 10, 
                ngene = 10000, 
                nsig = 25,
                total_sig = 300,
                nfeature = 500,
                log2fc = 2.4,
                ncore = 20)
set_10cl$ncell = rep(400, set_10cl$K)
set_10cl$clust.label.real = rep(1:set_10cl$K, set_10cl$ncell)
para.set = set_10cl
set.indx = "sim_10cl3"
para.set$N0 = sum(para.set$ncell)

set.seed(1)


## mean and dispersion parameter
generic_mu = 0.2
mu = matrix(rep(round(rtruncnorm(para.set$ngene, a=0.001, b=1, mean=generic_mu, sd=0.1), 4), para.set$K), 
            nrow=para.set$ngene)
r = matrix(rep(round(rtruncnorm(para.set$ngene, a =0.2, b=5, mean=0.5, sd=0.1), 4), para.set$K), 
           nrow=para.set$ngene)


## signature genes
ns = para.set$nsig
nz = round(para.set$nsig/3)
ls.sig = split(1:(para.set$nsig * para.set$K), 
               f = factor(rep(1:para.set$K, rep(para.set$nsig,para.set$K))))

for(j in 1:para.set$K){
  sig_g_indx = ls.sig[[j]]
  non_sig_celltype = setdiff(1:para.set$K, j)
  zero_col = sample(non_sig_celltype, round(length(non_sig_celltype)/2))
  zero_g_indx = sample(ls.sig[[j]], nz)
  mu[zero_g_indx, zero_col] = 0
}

for(j in 1:para.set$K){
  mu[ls.sig[[j]],j] = pmax(mu[ls.sig[[j]],j] *(2^para.set$log2fc), (generic_mu/2) * (2^para.set$log2fc)) 
}

for(j in 1:para.set$K){
  r[ls.sig[[j]],j] = r[ls.sig[[j]],j] + rnorm(ns, 0, sd = 0.05) 
}


## count matrix

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

## Seurat
ncl_vec = c(4,          6,       7,        8,         9,        10,    11,    12,   14,    16)
seu_res = c(0.024,   0.051,   0.062,   0.078,   0.093,   2.000,   3.000,      4,     5.106,  5.366)
## The resoultion for a specific number of clusters is found by binary search
## unlist(lapply(tmp, BS_resolution, gcmat = X, res_seq = seq(3.001, 8, by = 0.005)))
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


## combine all labels as a data frame

method_name = c("CIDR", "KMeans", "HC", "SC3", "Seurat", "Spectral")
lab_df = data.frame(tmp = rep(NA, sum(set_10cl$ncell)))
for(k in seq_len(length(method_name))){
  for(i in seq_len(length(ncl_vec))){
    cur_df = readRDS(paste0(set.indx, "_", method_name[k], "_k", ncl_vec[i], "_labs.rds"))
    colnames(cur_df) = paste0(method_name[k], "_k", ncl_vec[i])
    lab_df = cbind(lab_df, cur_df)
  }
}
lab_df = lab_df[,-1]

# ------------------------------------------------------------------------------
#                           Feature Selection 
# ------------------------------------------------------------------------------



selected_feature = feature_gene_selection(
  gcmat = X, 
  Seurat_obj = NULL,
  method = "wds",
  nfeature = 500,
  batch_label = NULL,
  zp_threshold = 0.95)

X_sub = X[selected_feature, ]


# ------------------------------------------------------------------------------
#                           Calculate CDI data frame
# ------------------------------------------------------------------------------

# library(BiocParallel)
# library(stringr)
X_sc = size_factor(X)


aicbic_df = calculate_CDI(
  sub_gcmat = X_sub,
  Seurat_obj = NULL,
  cand_lab_df = lab_df, 
  cell_size_factor = X_sc, 
  batch_label = NULL,
  lrt_pval_threshold = 0.01,
  clustering_method = NULL,
  ncore = 20)

saveRDS(aicbic_df, "~/to_save/sd1/sd1_cdi_df.rds")



cdi_df = aicbic_df %>% 
	mutate(AIC = AIC/(1e+06), BIC = BIC/(1e+06))


start_time = Sys.time()	
benchmark_return2 = calculate_CDI(
  sub_gcmat = X_sub, 
  cand_lab_df = para.set$clust.label.real, 
  batch_label = NULL, 
  cell_size_factor = X_sc, 
  ncore = 20)
end_time = Sys.time()
difftime(end_time, start_time)


saveRDS(benchmark_return2, "~/to_save/sd1/sd1_benchmark_return.rds")
# benchmark_return = readRDS("~/to_save/sd1/sd1_benchmark_return.rds")




# ------------------------------------------------------------------------------
#                            Figures (lineplot & Contingency) 
# ------------------------------------------------------------------------------









pdf("sd1_bic_lineplot.pdf", width = 3, height = 3.5)
CDILineplot(cdi_df, "CDI_BIC",
	benchmark_maintype_ncluster = length(unique(para.set$clust.label.real)), 
	benchmark_maintype_cdi = lapply(benchmark_return, function(x) x/10^6), 
	show_axis_names = FALSE, 
	show_method_legend = FALSE) +
	scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
		rgb(136, 146, 179, max = 255), #light purpleblue -- HC
    rgb(112, 186, 211, max = 255), #lightblue -- K-means
		rgb(229, 160, 133, max = 255),#pink -- SC3
    "#BCBD22B2",#-- Seurat
		rgb(76, 158, 137, max = 255))) # Spectral
dev.off()






LabComp = function(trlab, givlab){
  if(any(givlab == "0")){
    labs_new = as.character(as.numeric(as.character(givlab)) +1)
  } else {
    labs_new = givlab
  }
  ct.mat = as.matrix(table(data.frame(trlab, givlab = labs_new)))
  prop = t(round(ct.mat/matrix(rep(colSums(ct.mat), nrow(ct.mat)), 
                               nrow = nrow(ct.mat), byrow = T),2))
  prop_indx = cbind(prop, as.numeric(rownames(prop)))
  prop_sort = prop_indx[order(prop_indx[,ncol(prop_indx)]), c(1:ncol(prop))]
  longData<-reshape2::melt(prop_sort)
  longData<-longData[longData$value!=0,]
  colnames(longData)[1:2] = c("givlab", "trlab")
  order_level = NULL
  for(cell_type in levels(as.factor(longData$trlab))){
    large_prop_indx = which(prop_sort[,cell_type] >= 0.2)
    large_prop_sort = large_prop_indx[order(prop[large_prop_indx, cell_type], 
                                            decreasing = TRUE)]
    distinct_indx = setdiff(large_prop_sort, order_level)
    order_level = c(order_level, distinct_indx)
  }
  p<- ggplot(longData, aes(x = factor(givlab, levels = order_level), 
                           y = factor(trlab), fill=value)) + 
    geom_tile() + 
    scale_fill_gradientn(limits = c(0,1),
                         colours = c(rgb(204,204,204, maxColorValue = 255), 
                                     rgb(25,150,125, maxColorValue = 255))) + 
  	scale_x_discrete(labels=c(1:length(unique(longData$givlab)))) +
    labs(fill = "Proportion") +  
    theme_classic() + theme(axis.text = element_text(size=15, angle=0, vjust=0.3), 
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank(),
                            axis.line = element_line(size = 1), 
                            #axis.title=element_text(size=17), 
                            legend.position = "none")
  return(list(p, order_level))
}


pdf("sd1_bicaicopt_contingency.pdf", width = 2.5, height = 2.5)
LabComp(as.character(para.set$clust.label.real), 
	as.character(readRDS(paste0(set.indx, "_Seurat_k10_labs.rds"))[,1]))

dev.off()


# ------------------------------------------------------------------------------
#                           			Other metrics
# ------------------------------------------------------------------------------








# ------------------------------------------------------------------------------
#                      UMAPS from Selected Features
# ------------------------------------------------------------------------------




