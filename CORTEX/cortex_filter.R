library(Matrix)
library(stats)
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
library(RColorBrewer)
library(BiocParallel)

#### Data source:

## https://github.com/mohuangx/SAVER-paper
## under the SAVER-data folder, it provides the google drive link to datasets:
## https://drive.google.com/drive/folders/1uos3DDymyrh3ZyxveDhqQIcOyD1brxeL?usp=sharing


hrvatin_info = read.csv("hrvatin_celltypes.csv")
hrvatin = readRDS("hrvatin.rds")
# dim(hrvatin) = 19155, 10000
hrvatin_cell_name = data.frame(X = colnames(hrvatin))
merged_info = merge(hrvatin_cell_name, hrvatin_info, sort = F)


labeled_maintype_indx = c(1:nrow(merged_info))[!is.na(merged_info$maintype)]
labeled_maintype_info = merged_info[labeled_maintype_indx,]
rownames(labeled_maintype_info) = as.character(labeled_maintype_info$X)
hrvatin_process = hrvatin[, rownames(labeled_maintype_info)]
# > dim(hrvatin_process)
# [1] 19155  7390


qc1ht = function(scdata){
  nzg_of_each_cell = colSums(scdata > 0)
  ng = nrow(scdata)
  return(scdata[, (nzg_of_each_cell > 300) | (nzg_of_each_cell/ng > 0.03)])
}


qc2ht = function(scdata){
  nzc_of_each_gene = rowSums(scdata>0)
  nc = ncol(scdata)
  select_gene = c(1:nrow(scdata))[(nzc_of_each_gene > 50) | (nzc_of_each_gene/nc > 0.01)]
  return(select_gene)
}


hrvatin_process = qc1ht(hrvatin_process)
hrvatin_process = hrvatin_process[qc2ht(hrvatin_process),]
# gene remain: 12887
# cell remain: 7376
maintype_labels = labeled_maintype_info[colnames(hrvatin_process), "maintype"]
celltype_labels = labeled_maintype_info[colnames(hrvatin_process), "celltype"]


X = hrvatin_process
para.set = list(K = length(unique(maintype_labels)),
                ncell = as.numeric(table(maintype_labels)),
                subK = length(unique(celltype_labels)),
                subncell = as.numeric(table(celltype_labels)),
                nfeature = 500,
                ngene = nrow(X),
                celltype = maintype_labels,
                subtype = celltype_labels)

set.indx = "cortex"

saveRDS(X, paste0(set.indx, "_filtered.rds"))
saveRDS(para.set, paste0(set.indx, "_paraset.rds"))



