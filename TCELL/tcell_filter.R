library(Matrix)
library(matrixStats)
library(MASS)
library(dplyr)

# ------------------------------------------------------------------------------
#                   Combine three LSC datasets
# ------------------------------------------------------------------------------


rm(list = ls())

## pick 5 distinct clusters
lsc_celltype = list(
  lsc1 = c("CD8 Trm cells", "regulatory Trm cells"),
  lsc2 = c("CD8 Tcm, IL17RA+ & CD28+"),
  lsc3 = c("Active EM-like Treg","Classical CD4 Tem")
)

## dataset lsc1
ds_indx = 1
data_path = paste0("~/lsc12356/lsc_datasets/lsc", ds_indx)
mat_raw = as.matrix(readMM(paste0(data_path, "/matrix.mtx")))
barcodes = read.table(paste0(data_path, "/barcodes.tsv"), sep = "\t")
cell_labels = read.csv(paste0(data_path, "/LSC",ds_indx,"_meta_group.csv"))
cell_labels$Cell_name = paste0(cell_labels$Cell_name, "-1")

barcode_label1 = barcodes %>% 
	mutate(dataset = rep("lsc1", ncol(mat_raw))) %>%
	rename("Cell_name" = "V1") %>%
	mutate(Indx = c(1:ncol(mat_raw))) %>%
	left_join(cell_labels, by = "Cell_name") %>%
	filter(!is.na(Group)) %>%
	filter(Group %in% lsc_celltype[[ds_indx]])

gene_indx = read.table(paste0(data_path, "/genes.tsv"), sep = "\t")
colnames(gene_indx) = c("EnsemblID", "GeneName")
lsc1_mat = mat_raw[, barcode_label1$Indx]
rownames(lsc1_mat) = gene_indx$GeneName
colnames(lsc1_mat) = barcode_label1$Cell_name



## dataset lsc2
ds_indx = 2
data_path = paste0("~/lsc12356/lsc_datasets/lsc", ds_indx)
mat_raw = as.matrix(readMM(paste0(data_path, "/matrix.mtx")))
barcodes = read.table(paste0(data_path, "/barcodes.tsv"), sep = "\t")
cell_labels = read.csv(paste0(data_path, "/LSC",ds_indx,"_meta_group.csv"))
cell_labels$Cell_name = paste0(cell_labels$Cell_name, "-1")

barcode_label2 = barcodes %>% 
	mutate(dataset = rep("lsc2", ncol(mat_raw))) %>%
	rename("Cell_name" = "V1") %>%
	mutate(Indx = c(1:ncol(mat_raw))) %>%
	left_join(cell_labels, by = "Cell_name") %>%
	filter(!is.na(Group)) %>%
	filter(Group %in% lsc_celltype[[ds_indx]])

gene_indx = read.table(paste0(data_path, "/genes.tsv"), sep = "\t")
colnames(gene_indx) = c("EnsemblID", "GeneName")
lsc2_mat = mat_raw[, barcode_label2$Indx]
rownames(lsc2_mat) = gene_indx$GeneName
colnames(lsc2_mat) = barcode_label2$Cell_name


## dataset lsc3
ds_indx = 3
data_path = paste0("~/lsc12356/lsc_datasets/lsc", ds_indx)
mat_raw = as.matrix(readMM(paste0(data_path, "/matrix.mtx")))
barcodes = read.table(paste0(data_path, "/barcodes.tsv"), sep = "\t")
cell_labels = read.csv(paste0(data_path, "/LSC",ds_indx,"_meta_group.csv"))
cell_labels$Cell_name = paste0(cell_labels$Cell_name, "-1")

barcode_label3 = barcodes %>% 
	mutate(dataset = rep("lsc3", ncol(mat_raw))) %>%
	rename("Cell_name" = "V1") %>%
	mutate(Indx = c(1:ncol(mat_raw))) %>%
	left_join(cell_labels, by = "Cell_name") %>%
	filter(!is.na(Group)) %>%
	filter(Group %in% lsc_celltype[[ds_indx]])

gene_indx = read.table(paste0(data_path, "/genes.tsv"), sep = "\t")
colnames(gene_indx) = c("EnsemblID", "GeneName")
lsc3_mat = mat_raw[, barcode_label3$Indx]
rownames(lsc3_mat) = gene_indx$GeneName
colnames(lsc3_mat) = barcode_label3$Cell_name

# sum(rownames(lsc1_mat) != rownames(lsc2_mat)) 
# [1] 0
# sum(rownames(lsc1_mat) != rownames(lsc3_mat))
# [1] 0


tcell_mat = cbind(cbind(lsc1_mat, lsc2_mat), lsc3_mat)
tcell_cinfo = rbind(rbind(barcode_label1, barcode_label2), barcode_label3)

# ------------------------------------------------------------------------------
#                   						Filtering
# ------------------------------------------------------------------------------

PreprocessSelectCells = function(scdata, min_nzg = 300, min_nzg_prop = 0.03){
  nzg_of_each_cell = colSums(scdata > 0)
  ng = nrow(scdata)
  # return submatrix of cells
  return(c(1:ncol(scdata))[(nzg_of_each_cell > min_nzg) & (nzg_of_each_cell/ng > min_nzg_prop)])
}


PreprocessSelectGenes = function(scdata, min_nzc = 100, min_nzc_prop = 0.02){
  nzc_of_each_gene = rowSums(scdata>0)
  nc = ncol(scdata)
  select_gene = c(1:nrow(scdata))[(nzc_of_each_gene > min_nzc) & (nzc_of_each_gene/nc > min_nzc_prop)]
  # return gene index
  return(select_gene)
}



selected_genes = PreprocessSelectGenes(tcell_mat,  min_nzc = 0, min_nzc_prop = 0.02)
tcell_filtergenes = tcell_mat[selected_genes,]
# 7893 genes retained
selected_cells = PreprocessSelectCells(tcell_filtergenes,  min_nzg = 0, min_nzg_prop = 0.02)
tcell_filtercells = tcell_filtergenes[,selected_cells]
# all 2989 cells retained
saveRDS(tcell_filtercells, "Tcell_5type_filtered.rds")
saveRDS(tcell_cinfo, "Tcell_5type_filtered_labels.rds")












