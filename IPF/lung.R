# ------------------------------------------------------------------------------
#                            Libraries 
# ------------------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(dplyr)
library(matrixStats)
library(truncnorm)
library(data.table)
library(mclust)
library(ggsci)
library(cluster)
library(BiocParallel)



# ------------------------------------------------------------------------------
#                           Data process
# ------------------------------------------------------------------------------


# GSE135893
# paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7439444/
# BEC: Notably, we did not observe over batch effects driven by processing site or sequencing batch in our dimensionality reduction and visualization (fig. S2).
# annotation: Seurat(all & sub) & celltype markers
# Data Processing: "Individual sample output files from CellRanger Count v3.0.2 (GRCh38) were read into Seurat v3 and combined into one Seurat object. 
# Cells containing less than 1,000 nFeature_RNA and more than 25% percentage of mitochondrial genes were filtered out. 
# SCTransform with default parameters was used to normalize and scale the data, and dimensionality reduction was performed using PCA and UMAP. 
# To group cells into clusters, FindClusters function was used with a resolution of 0.01, resulting in six major clusters."

info = read.csv('GSE135893_IPF_metadata.csv')
expression_matrix <- ReadMtx(
  mtx = "GSE135893_matrix.mtx", 
  features = "GSE135893_genes.tsv",
  feature.column = 1,
  cells = "GSE135893_barcodes.tsv"
)
seurat_object <- CreateSeuratObject(counts = expression_matrix, min.cells = 100)
# > dim(seurat_object[["RNA"]]@counts)
# [1]  20354 220213
seurat_object@meta.data['barcode'] <- rownames(seurat_object@meta.data)
annotated_object <- subset(seurat_object, barcode %in% info$X)
# > dim(annotated_object[["RNA"]]@counts)
# [1]  20354 114396
rm(seurat_object)


barcode_df = data.frame(barcode = annotated_object@meta.data['barcode'])  
cell_info = left_join(barcode_df, info, by = c('barcode' = 'X'))
annotated_object@meta.data['maintype'] = cell_info$population
annotated_object@meta.data['celltype'] = cell_info$celltype
annotated_object@meta.data['Sample_Name'] = cell_info$Sample_Name
annotated_object@meta.data['Sample_Source'] = cell_info$Sample_Source
annotated_object@meta.data['Status'] = cell_info$Status
annotated_object@meta.data['Diagnosis'] = cell_info$Diagnosis
object.size(annotated_object)
saveRDS(annotated_object, "lung_processed_annotated_robj.rds")



# ------------------------------------------------------------------------------
#                           Generate UMAPs with benchmark labels
# ------------------------------------------------------------------------------

# select a subset, run integration and generate UMAP
seu_obj <- readRDS("lung_processed_annotated_robj.rds")
head(seu_obj@meta.data)
seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
seu_obj <- ScaleData(seu_obj)
seu_obj <- RunPCA(seu_obj, npcs = 30, verbose = FALSE)
seu_obj <- RunUMAP(seu_obj, reduction = "pca", dims = 1:20)
Idents(seu_obj) <- seu_obj@meta.data['maintype']
pdf("lung_maintype_umap.pdf", width = 7, height = 6)
DimPlot(seu_obj, reduction = "umap")
dev.off()

Idents(seu_obj) <- seu_obj@meta.data['celltype']
pdf("lung_subtype_umap.pdf", width = 10, height = 6)
DimPlot(seu_obj, reduction = "umap")
dev.off()



Idents(seu_obj) <- ifelse(seu_obj@meta.data['celltype'] == "Macrophages", "Macrophages", "Others")
pdf("lung_subtype_umap_mac.pdf", width = 8, height = 6)
DimPlot(seu_obj, reduction = "umap")
dev.off()

Idents(seu_obj) <- ifelse(seu_obj@meta.data['celltype'] == "T Cells", "T Cells", "Others")
pdf("lung_subtype_umap_tcell.pdf", width = 8, height = 6)
DimPlot(seu_obj, reduction = "umap")
dev.off()







Start = Sys.time()
seu_obj <- FindNeighbors(seu_obj, dims = 1:20)
# seu_obj <- FindClusters(seu_obj, resolution = 0.5)
End = Sys.time()
difftime(End, Start, units = "secs")
## only takes 57 secs


## Generate labels for Seurat 





BS_resolution2 = function(ideal_K, Seurat.obj, res_seq = seq(0.0001, 0.0099, by = 0.0001)){
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

# ncl =  100 res= 5.97
# ncl =  95 res= 5.93
# ncl =  90 res= 5.12
# ncl =  85 res= 4.71
# ncl =  80 res= 4.231
# ncl =  77 res= 4.2
# ncl =  76 res= 3.86
# ncl =  75 res= 4.22
# ncl =  74 res= 3.87 
# ncl =  70 res= 3.52
# ncl =  68 res= 3.64
# ncl =  65 res= 3.18
# ncl =  58 res= 2.66 
# ncl =  56 res= 2.57 
# ncl =  54 res= 2.53 
# ncl =  51 res= 2.12
# ncl =  50 res= 2.06
# ncl =  49 res= 2.0015 
# ncl =  48 res= 1.9304 
# ncl =  45 res= 1.8
# ncl =  44 res= 1.7664
# ncl =  42 res= 1.6457 
# ncl =  41 res= 1.6
# ncl =  40 res= 1.361
# ncl =  36 res= 1.25 
# ncl =  35 res= 1.091
# ncl =  34 res= 1.011
# ncl =  33 res= 0.933 
# ncl =  32 res= 0.816 
# ncl =  31 res= 0.807
# ncl =  30 res= 0.777
# ncl =  29 res= 0.649
# ncl =  28 res= 0.62 
# ncl =  25 res= 0.464 
# ncl =  21 res= 0.31 
# ncl =  20 res= 0.234 
# ncl =  15 res= 0.16 
# ncl =  12 res= 0.08 
# ncl = 10 res= 0.06055
# ncl =  9 res= 0.04 
# ncl =  8 res= 0.02 
# ncl =  8 res= 0.01   
# ncl =  7 res= 0.005 
# ncl =  6 res= 0.00361
# ncl =  5 res= 0.0025
# ncl = 4 res = 0.0012


ncl_vec = c( 5,     10,       15,   20,    25,    30,     35,   40,          45,  50,    55,    60,     65,    70,   75,    80,    85,    90,    95,    100)
seu_res = c(0.0025,  0.06055, 0.16, 0.234, 0.464, 0.777,  1.091,   1.361,    1.8, 2.06,  2.63, 2.939,   3.18,  3.52, 4.22,  4.231, 4.71,  5.12,  5.93,  5.97)

seurat_lab2 = function(set.indx, seu_obj, res_vec){
  label_result = matrix(NA, nrow = nrow(seu_obj@meta.data), ncol = length(res_vec))
  num_cl = numeric(length(res_vec))
  for(i in 1:length(res_vec)){
    seu_obj = FindClusters(seu_obj, resolution = res_vec[i], verbose = F)
    df = data.frame(Seurat = seu_obj@meta.data$seurat_clusters)
    num_cl[i] = length(unique(df$Seurat))
    #saveRDS(df, paste0(set.indx, "_Seurat_k", num_cl[i], "_labs.rds"))
    saveRDS(df, paste0(set.indx, "_Seurat_k", num_cl[i], "_labs.rds"))
    seu_obj@meta.data$seurat_clusters = NULL
  }
  return(data.frame(ncl = num_cl, res = res_vec))
}


seurat_lab2("lung", seu_obj, res_vec = seu_res)


ncl_vec = c(70,   75,    80,    85)
seu_res = c(3.52, 4.22,  4.231, 4.71)

seurat_lab2("lung", seu_obj, res_vec = seu_res)

# ncl =  100 res= 5.97
# ncl =  95 res= 5.93
# ncl =  90 res= 5.12
ncl_vec = c(90,    95,    100)
seu_res = c(5.12,  5.93,  5.97)

seurat_lab2("lung", seu_obj, res_vec = seu_res)



set.indx = "lung"
method_name = c("Seurat")
lab_df = data.frame(tmp = readRDS(paste0(set.indx, "_", method_name[1], "_k", ncl_vec[1], "_labs.rds"))[,1])
for(k in seq_len(length(method_name))){
  for(i in seq_len(length(ncl_vec))){
    cur_df = readRDS(paste0(set.indx, "_", method_name[k], "_k", ncl_vec[i], "_labs.rds"))
    colnames(cur_df) = paste0(method_name[k], "_k", ncl_vec[i])
    lab_df = cbind(lab_df, cur_df)
  }
}
lab_df = lab_df[,-1]
# saveRDS(lab_df, paste0(set.indx, "_clustering_labels.rds"))



# ----------------------------------------------
#            Write MM files
# ----------------------------------------------

ng = 20354
nc = 114396
ng_file = 102; ng_size = 200
nc_file = 58; nc_size = 2000

seu = readRDS("lung_processed_annotated_robj.rds")
for(i in seq_len(ng_file - 1)){
  writeMM(seu@assays$RNA@counts[((i-1)*ng_size+1):(i*ng_size), ], paste0("./lung_split_gene/lung_gene_p", i, ".txt"))
}
writeMM(seu@assays$RNA@counts[((ng_file-1)*ng_size):ng, ], paste0("./lung_split_gene/lung_gene_p", ng_file, ".txt"))

for(j in seq_len(nc_file - 1)){
  writeMM(seu@assays$RNA@counts[,((j-1)*nc_size+1):(j*nc_size)], paste0("./lung_split_cell/lung_cell_p", j, ".txt"))
}
writeMM(seu@assays$RNA@counts[,((nc_file-1)*nc_size):nc], paste0("./lung_split_cell/lung_cell_p", nc_file, ".txt"))



# ----------------------------------------------
#            Size factor
# ----------------------------------------------

# tried readLines(), it will ignore the sparse matrix form

gene_ref = c()
for(i in seq_len(ng_file)){
  X = readMM(paste0("./lung_split_gene/lung_gene_p", i, ".txt"))
  gene_ref = c(gene_ref, exp(rowMeans(log(pmax(as.matrix(X), 0.5)))))
}
names(gene_ref) = seu@assays$RNA@counts@Dimnames[[1]]
saveRDS(gene_ref, "lung_gene_ref_size.rds")

sc_vec = c()
for(j in seq_len(nc_file)){
  if(j %% 5 == 0){cat("j = ", j, "\n")}
  X = readMM(paste0("./lung_split_cell/lung_cell_p", j, ".txt"))
  sc_vec = c(sc_vec, colMedians(sweep(pmax(as.matrix(X), 0.5), 1, gene_ref, "/")))
}
saveRDS(sc_vec, "lung_size_factor.rds")


# ----------------------------------------------
#            Feature Gene Selection
# ----------------------------------------------

feature_df = data.frame()
for(i in seq_len(ng_file)){
  if(i %% 5 == 0){cat("i = ", i, "\n")}
  X = as.matrix(readMM(paste0("./lung_split_gene/lung_gene_p", i, ".txt")))
  mu = rowMeans(X)
  var = rowVars(X)
  gene_zp = rowMeans(X == 0)
  gene_phi = ifelse(gene_zp > 0.95, NA, (var - mu)/mu^2)
  tmp_feature = data.frame(zp = gene_zp, phi = gene_phi)
  feature_df = rbind(feature_df, tmp_feature)
}
rownames(feature_df)= seu@assays$RNA@counts@Dimnames[[1]]


tmp_phi = ifelse(is.na(feature_df$phi), -10000, feature_df$phi)
phi_rank = rank(-tmp_phi)
feature_df["rank"] = phi_rank
saveRDS(feature_df, "lung_feature_df.rds")

feature_gene_indx = which(feature_df$rank <= 500)

X_sub = seu@assays$RNA@counts[feature_gene_indx, ]
writeMM(X_sub, "lung_X_sub.txt")

# ----------------------------------------------
#            Calculate CDI
# ----------------------------------------------
X_sc = readRDS("lung_size_factor.rds")
cinfo = readRDS("lung_cell_info.rds")
X_sub = as.matrix(readMM("lung_X_sub.txt"))


Start = Sys.time()
celltype_return = calculate_CDI(
  sub_gcmat = X_sub,
  Seurat_obj = NULL,
  cand_lab_df = cinfo$celltype, 
  cell_size_factor = X_sc, 
  batch_label = NULL,
  BPPARAM = MulticoreParam(ncore))
End = Sys.time()
End - Start 

# Time difference of 8.164048 mins
saveRDS(celltype_return, "lung_benchmark_return.rds")



## more efficiently, put this part to cluster with bash script
library(CDI)
ncore = 10
aicbic_df = calculate_CDI(
  sub_gcmat = X_sub,
  Seurat_obj = NULL,
  cand_lab_df = lab_df,
  cell_size_factor = X_sc,
  batch_label = NULL,
  BPPARAM = MulticoreParam(ncore))



saveRDS(aicbic_df, "lung_cdi_df.rds")


# ----------------------------------------------
#            Figures
# ----------------------------------------------

cdi_df = readRDS("lung_cdi_df.rds")
benchmark_return = readRDS("lung_benchmark_return.rds")


pdf("lung_aic_lineplot.pdf", width = 3, height = 3.5)
CDI_lineplot(cdi_dataframe = cdi_df %>% mutate(CDI_AIC = CDI_AIC/10^6, CDI_BIC = CDI_BIC/10^6), 
             cdi_type = "CDI_AIC",
             benchmark_celltype_cdi = lapply(benchmark_return, function(x) x/10^6),  
             benchmark_celltype_ncluster  = 31,
             show_axis_names = FALSE, 
             show_method_legend = FALSE) +
  scale_color_manual(values = c("#BCBD22B2")) 
dev.off()


pdf("lung_bic_lineplot.pdf", width = 3, height = 3.5)
CDI_lineplot(cdi_dataframe = cdi_df %>% mutate(CDI_AIC = CDI_AIC/10^6, CDI_BIC = CDI_BIC/10^6), 
             cdi_type = "CDI_BIC",
             benchmark_celltype_cdi = lapply(benchmark_return, function(x) x/10^6),  
             benchmark_celltype_ncluster  = 31,
             show_axis_names = FALSE, 
             show_method_legend = FALSE) +
  scale_color_manual(values = c("#BCBD22B2")) 
dev.off()



contingency_heatmap2 <- function(benchmark_label, 
                                 candidate_label, 
                                 proportion_size = 4,
                                 show_axis_names = TRUE,
                                 show_legend = TRUE,
                                 rename_candidate_clusters = FALSE,
                                 candidate_cluster_names = NULL){
  ## create contingency table
  count_mat <- as.matrix(table(data.frame(benchmark_label, candidate_label)))
  ## change from counts to proportions
  col_sum_matrix = matrix(rep(colSums(count_mat), nrow(count_mat)), nrow = nrow(count_mat), byrow = TRUE)
  ## each row is a candidate class, and each row sums to 1.
  prop_mat <- t(round(count_mat / col_sum_matrix, 2))
  longData <- reshape2::melt(prop_mat)
  ## To avoid the problem of undifined variable
  PropValue = NULL
  colnames(longData)[3] <- "PropValue"
  ## remove zero proportions
  longData <- longData[longData$PropValue != 0, ]
  ## order benchmark cell types according to the composition of candidate cell types
  order_level <- NULL
  for(cell_type in levels(as.factor(longData$benchmark_label))){
    large_prop_indx <- names(which(prop_mat[,cell_type] >= 0.2))
    large_prop_sort <- large_prop_indx[order(prop_mat[large_prop_indx, cell_type], decreasing = TRUE)]
    distinct_indx <- setdiff(large_prop_sort, order_level)
    order_level <- c(order_level, distinct_indx)
  }
  p1 <- ggplot(longData, 
               aes(x = factor(benchmark_label), y = factor(candidate_label, levels = order_level), fill = PropValue)) + 
    scale_fill_gradientn(limits = c(0,1), 
                         colours = c(rgb(204,204,204, maxColorValue = 255), rgb(25,150,125, maxColorValue = 255))) + 
    labs(x = "Benchmark label", y = "Candidate label", fill = "Proportion") + 
    theme_classic() + geom_tile() + 
    theme(axis.text.x = element_text(angle = 45, hjust=1, size=10), 
          axis.text.y = element_text(size=10, angle=0, vjust=0.3), 
          axis.line = element_line(size = 1))
  if(rename_candidate_clusters){
    p1 <- p1 + scale_x_discrete(labels = candidate_cluster_names)
  }
  if(proportion_size){
    p1 <- p1 + geom_text(label = round(longData$PropValue, 2), 
                         size = proportion_size) 
  }
  if(show_axis_names == FALSE){
    p1 <- p1 + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  }
  if(show_legend == FALSE){
    p1 <- p1 + theme(legend.position = "none")
  }
  return(p1)
}



pdf("lung_seurat80_contingency.pdf", height =  15, width = 10)
contingency_heatmap2(benchmark_label = seu_obj@meta.data$celltype, candidate_label =  seu_obj@meta.data$seurat_clusters, proportion_size = 2.5)
dev.off()


# ------------------------------------------------------------------------------
#                          UMAP for sub-pop
# ------------------------------------------------------------------------------

seu_obj <- readRDS("lung_processed_annotated_robj.rds")
head(seu_obj@meta.data)
seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
seu_obj <- ScaleData(seu_obj)
seu_obj <- RunPCA(seu_obj, npcs = 30, verbose = FALSE)
seu_obj <- RunUMAP(seu_obj, reduction = "pca", dims = 1:20)
seu_obj <- FindNeighbors(seu_obj, dims = 1:20)
seu_obj <- FindClusters(seu_obj, resolution = 4.231)
#seu_obj@meta.data["s80"] <- readRDS("lung_Seurat_k80_labs.rds")[,1]
#Idents(seu_obj) <- readRDS("lung_Seurat_k80_labs.rds")[,1]




# library(ggsci)
info = readRDS("lung_cell_info.rds")
ct_name = c("AT2", "Ciliated", "Endothelial Cells", "Macrophages", "T Cells", "MUC5B+")
info = info %>% mutate(seurat_lab = seu_obj@meta.data$seurat_clusters)

tmp_df1 = info %>% filter(celltype == ct_name[1]) 
sub_names = rownames(tmp_df1)[c(tmp_df1$seurat_lab %in% names(which(table(tmp_df1$seurat_lab) >= 5)))]
sub_seu1 = subset(seu_obj, barcode %in% sub_names)


tmp_df2 = info %>% filter(celltype == ct_name[2]) 
sub_names = rownames(tmp_df2)[c(tmp_df2$seurat_lab %in% names(which(table(tmp_df2$seurat_lab) >= 5)))]
sub_seu2 = subset(seu_obj, barcode %in% sub_names)



tmp_df3 = info %>% filter(celltype == ct_name[3]) 
sub_names = rownames(tmp_df3)[c(tmp_df3$seurat_lab %in% names(which(table(tmp_df3$seurat_lab) >= 5)))]
sub_seu3 = subset(seu_obj, barcode %in% sub_names)

tmp_df4 = info %>% filter(celltype == ct_name[4]) 
sub_names = rownames(tmp_df4)[c(tmp_df4$seurat_lab %in% names(which(table(tmp_df4$seurat_lab) >= 5)))]
sub_seu4 = subset(seu_obj, barcode %in% sub_names)


tmp_df5 = info %>% filter(celltype == ct_name[5]) 
sub_names = rownames(tmp_df5)[c(tmp_df5$seurat_lab %in% names(which(table(tmp_df5$seurat_lab) >= 5)))]
sub_seu5 = subset(seu_obj, barcode %in% sub_names)


tmp_df6 = info %>% filter(celltype == ct_name[6]) 
sub_names = rownames(tmp_df6)[c(tmp_df6$seurat_lab %in% names(which(table(tmp_df6$seurat_lab) >= 5)))]
sub_seu6 = subset(seu_obj, barcode %in% sub_names)


p1 <- DimPlot(sub_seu1, reduction = "umap") + 
  scale_color_d3() + 
  labs(title = ct_name[1]) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "none") +
  guides(shape = guide_legend(override.aes = list(size = 0.2)))

p2 <- DimPlot(sub_seu2, reduction = "umap") + 
  scale_color_d3() + 
  labs(title = ct_name[2]) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "none") +
  guides(shape = guide_legend(override.aes = list(size = 0.2)))

p3 <- DimPlot(sub_seu3, reduction = "umap") + 
  scale_color_d3() + 
  labs(title = ct_name[3]) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "none") +
  guides(shape = guide_legend(override.aes = list(size = 0.2)))

p4 <- DimPlot(sub_seu4, reduction = "umap") + 
  scale_color_manual(values = c(pal_d3(palette = "category20c")(20), pal_d3(palette = "category20b")(20)[c(3,4,5,8,9,10,13, 14, 15, 18,19,20)])) + 
  labs(title = ct_name[4]) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "none") + 
  guides(shape = guide_legend(override.aes = list(size = 0.01)))

p5 <- DimPlot(sub_seu5, reduction = "umap") + 
  scale_color_d3() + 
  labs(title = ct_name[5]) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "none")  +
  guides(shape = guide_legend(override.aes = list(size = 0.2)))


p6 <- DimPlot(sub_seu6, reduction = "umap") + 
  scale_color_d3() + 
  labs(title = ct_name[6]) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        legend.position = "none") +
  guides(shape = guide_legend(override.aes = list(size = 0.2)))


pdf("lung_heter_umap_0429.pdf", width = 12, height = 6)
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)
dev.off()
