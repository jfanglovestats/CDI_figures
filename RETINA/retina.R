library(Matrix)
library(matrixStats)
library(MASS)
library(BiocParallel)
library(truncnorm)
library(data.table)
library(mclust)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(grid)
library(pheatmap)
library(ggsci)
library(dplyr)



set.indx = "retina"
info = readRDS(paste0(set.indx, "_filtered_info.rds"))
X = readRDS(paste0(set.indx, "_filtered.rds"))
para.set = readRDS(paste0(set.indx, "_paraset.rds"))
para.set$ncore = 20

# ------------------------------------------------------------------------------
#                       Batch effect correction
# ------------------------------------------------------------------------------

Seurat.obj = CreateSeuratObject(counts = X, project = set.indx, min.cells = 0, min.features = 0)
Seurat.obj = NormalizeData(Seurat.obj)
Seurat.obj = FindVariableFeatures(Seurat.obj, selection.method = "vst", nfeatures = 2000)
Seurat.obj = ScaleData(Seurat.obj, verbose = FALSE)
Seurat.obj = RunPCA(Seurat.obj, npcs = 30, verbose = FALSE)
Seurat.obj = RunUMAP(Seurat.obj, reduction = "pca", dims = 1:20)

## colnames(X) and rownames(Seurat.obj@meta.data) are the same
Seurat.obj@meta.data["batch"] = para.set$Batch

##
Seurat.obj.list <- SplitObject(Seurat.obj, split.by = "batch")
Seurat.obj.list <- lapply(X = Seurat.obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


Seurat.obj.anchors <- FindIntegrationAnchors(object.list = Seurat.obj.list, 
                                             dims = 1:30, 
                                             anchor.features = 2000 #default setting
)

Seurat.obj.integrated <- IntegrateData(anchorset = Seurat.obj.anchors, 
                                       features.to.integrate = Seurat.obj.anchors@anchor.features, 
                                       dims = 1:30)
DefaultAssay(Seurat.obj.integrated)

saveRDS(Seurat.obj.integrated, "Seurat_integrated_obj.rds")

X_bc_unsorted = as.matrix(Seurat.obj.integrated[["RNA"]]@data)
col_indx = colnames(X)
X_bc = X_bc_unsorted[, col_indx]
saveRDS(X_bc, "retina_bc_countmat.rds")



# Seurat.obj.integrated = readRDS("Seurat_integrated_obj.rds")


# ------------------------------------------------------------------------------
#                       Generate Candidate Labels
# ------------------------------------------------------------------------------

# X_bc = readRDS("retina_bc_countmat.rds")
# Seurat.obj.integrated = readRDS("retina_seurat_integrated_obj.rds")
source("utils/clustering.R")
## Seurat
ncl_vec = c(5,     10,   15,     16,    17,     18,    19,    20,     25,     30,    31,    32,     33,       35)
seu_res = c(0.0035,0.01, 0.1,    0.18,   0.3,    0.5,    0.6,  0.71,   1.2,   2.15,  2.30,  2.60,  2.84,     3.20)
seurat_lab_bc(Seurat.obj.integrated, seu_res, colnames(X))


## SC3
start_time = Sys.time()
sc3_lab_large(X, ncl_vec)
end_time = Sys.time()
end_time - start_time

## CIDR
start_time = Sys.time()
cidr_lab(X, ncl_vec)
end_time = Sys.time()
end_time - start_time


# X_bc = readRDS("retina_bc_countmat.rds")
X_bc_pc = gcmat_pc(X_bc, npc = 200)
# saveRDS(X_bc_pc, "retina_bc_pc200_countmat.rds")
# X_bc_pc = readRDS("retina_bc_pc200_countmat.rds")

## spectral
# start_time = Sys.time()
# spec_lab_pc(X_bc_pc, ncl_vec)
# end_time = Sys.time()
# end_time - start_time
# failed to give output within 24h

## kmeans
start_time = Sys.time()
kmeans_lab_pc(X_bc_pc, ncl_vec)
end_time = Sys.time()
end_time - start_time
#52.65812 secs

## HC
start_time = Sys.time()
hc_lab_pc(X_bc_pc, ncl_vec)
end_time = Sys.time()
end_time - start_time
#4.596798 mins



## combine all labels as a data frame
set.indx = "retina"
method_name = c("CIDR", "KMeans", "HC", "SC3", "Seurat")
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
#                          Calculate CDI
# ------------------------------------------------------------------------------

library(CDI)
## feature selection
selected_feature = feature_gene_selection(
  gcmat = X, 
  Seurat_obj = NULL,
  method = "wds",
  nfeature = 500,
  batch_label = para.set$Batch,
  zp_threshold = 0.95)

X_sub = X[selected_feature, ]

# saveRDS(X_sub, paste0(set.indx, "_X_sub.rds"))


X_sc = size_factor(X)
saveRDS(X_sc, paste0(set.indx, "_size_factor.rds"))

cdi_df = calculate_CDI(
  sub_gcmat = X_sub,
  Seurat_obj = NULL,
  cand_lab_df = lab_df, 
  cell_size_factor = X_sc, 
  batch_label = info$Batch,
  lrt_pval_threshold = 0.01,
  clustering_method = NULL,
  BPPARAM = MulticoreParam(para.set$ncore))
saveRDS(cdi_df, paste0(set.indx, "_cdi_df.rds"))




start_time = Sys.time()	
maintype_return = calculate_CDI(sub_gcmat = X_sub,
                                Seurat_obj = NULL,
                                cand_lab_df = para.set$maintype, 
                                cell_size_factor = X_sc, 
                                batch_label = cell_info$Batch,
                                lrt_pval_threshold = 0.01,
                                clustering_method = NULL,
                                BPPARAM = MulticoreParam(para.set$ncore))
end_time = Sys.time()
difftime(end_time, start_time)

celltype_return =  calculate_CDI(sub_gcmat = X_sub,
                                 Seurat_obj = NULL,
                                 cand_lab_df = para.set$celltype, 
                                 cell_size_factor = X_sc, 
                                 batch_label = cell_info$Batch,
                                 lrt_pval_threshold = 0.01,
                                 clustering_method = NULL,
                                 BPPARAM = MulticoreParam(para.set$ncore))

saveRDS(list(benchmark_main = maintype_return, 
             benchmark_cell = celltype_return), 
        paste0(set.indx, "_benchmark_return.rds"))




# ------------------------------------------------------------------------------
#                         Figures (lineplot & Contingency) 
# ------------------------------------------------------------------------------


## CDI Lineplots

cur_info = info %>% 
  mutate(main_type = recode(cell_type, "RBC" = "Rod BC", 
                            "Muller_glia" = "Muller glia",
                            "Cone_photoreceptors" = "Photoreceptor",
                            "Rod_photoreceptors" = "Photoreceptor", 
                            "BC1A" = "OFF BC", "BC1B" = "OFF BC", 
                            "BC2" = "OFF BC", "BC3A" = "OFF BC",
                            "BC3B" = "OFF BC", "BC4" = "OFF BC",
                            "BC5A" = "ON BC", "BC5B" = "ON BC", 
                            "BC5C" = "ON BC", "BC5D" = "ON BC", 
                            "BC6" = "ON BC", "BC7" = "ON BC", 
                            "BC8_BC9" = "ON BC"),
         cell_type = recode(cell_type, "RBC" = "Rod BC", 
                            "Muller_glia" = "Muller glia", "Rod_PR" = "Rod PR", "Cone_PR" = "Cone PR"))




pdf(paste0(set.indx, "_aic_lineplot.pdf"), width = 3, height = 3.5)
CDI_lineplot(
  cdi_dataframe = cdi_df, 
  cdi_type = "CDI_AIC",
  benchmark_celltype_ncluster = length(unique(cur_info$cell_type)), 
  benchmark_celltype_cdi = lapply(benchmark_return[["benchmark_cell"]], function(x) x/10^6), 
  benchmark_maintype_ncluster = length(unique(cur_info$main_type)), 
  benchmark_maintype_cdi = lapply(benchmark_return[["benchmark_main"]], function(x) x/10^6)) +
  scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
                                rgb(136, 146, 179, max = 255), #light purpleblue -- HC
                                rgb(112, 186, 211, max = 255), #lightblue -- K-means
                                rgb(229, 160, 133, max = 255),#pink -- SC3
                                "#BCBD22B2")) #-- Seurat
dev.off()


pdf(paste0(set.indx, "_bic_lineplot.pdf"), width = 3, height = 3.5)
CDI_lineplot(
  cdi_dataframe = cdi_df, 
  cdi_type = "CDI_BIC",
  benchmark_celltype_ncluster = length(unique(cur_info$cell_type)), 
  benchmark_celltype_cdi = lapply(benchmark_return[["benchmark_cell"]], function(x) x/10^6), 
  benchmark_maintype_ncluster = length(unique(cur_info$main_type)), 
  benchmark_maintype_cdi = lapply(benchmark_return[["benchmark_main"]], function(x) x/10^6)) +
  scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
                                rgb(136, 146, 179, max = 255), #light purpleblue -- HC
                                rgb(112, 186, 211, max = 255), #lightblue -- K-means
                                rgb(229, 160, 133, max = 255),#pink -- SC3
                                "#BCBD22B2")) #-- Seurat
dev.off()



## Contingency table

ordered_maintype = c("OFF BC", "ON BC", "Rod BC", "Amacrine", "Muller glia", "Photoreceptor")
pdf(paste0("retina_bicopt_contingency.pdf"), width = 4, height = 2.5)
p <- contingency_heatmap(factor(cur_info$main_type, levels = ordered_maintype), 
                    as.character(readRDS(paste0("./clustering_labels/", set.indx, "_Seurat_k18_labs.rds"))[,1]))
p1 <- p + scale_x_discrete(labels=c("1", "", "",  "", "5", "", "", "", "", "10", 
                            "", "", "", "",  "15", "", "",  ""))
print(p1)
dev.off()


ordered_celltype = c("BC1A", "BC1B", "BC2", "BC3A", "BC3B", "BC4", "BC5A", "BC5B", "BC5C", 
                     "BC5D", 'BC6', 'BC7', "BC8_BC9",  "Rod BC", "Amacrine", "Muller glia", 
                     "Cone PR", "Rod PR")

pdf(paste0("retina_aicopt_contingency.pdf"), width = 4, height = 8)
p <- contingency_heatmap(factor(cur_info$cell_type, levels = ordered_celltype), 
                       as.character(readRDS(paste0("./clustering_labels/", set.indx, "_KMeans_k33_labs.rds"))[,1]))
p1 <- p + scale_x_discrete(labels=c("1", "", "",  "", "5", "", "", "", "", "10", 
                            "", "", "", "",  "15", "", "",  "", "", "20", "", "", "", "",  "25", 
                            "", "",  "", "", "30", "",  "",  ""))
print(p1)
dev.off()





# ------------------------------------------------------------------------------
#                           UMAP
# ------------------------------------------------------------------------------

Seurat.obj = CreateSeuratObject(counts = X, project = set.indx, min.cells = 0, min.features = 0)
Seurat.obj = NormalizeData(Seurat.obj)
Seurat.obj = FindVariableFeatures(Seurat.obj, selection.method = "vst", nfeatures = 500)
seu_sig_gname = Seurat.obj@assays[["RNA"]]@var.features
seu_sig_indx = as.numeric(sub("^.", "", Seurat.obj@assays[["RNA"]]@var.features))

ordered_ct = c("BC1A", "BC1B", "BC2", "BC3A", "BC3B", "BC4", "BC5A", "BC5B", "BC5C", 
               "BC5D", 'BC6', 'BC7', "BC8_BC9", "Amacrine", "Muller_glia", 
               "RBC", "Cone_PR", "Rod_PR")



library(umap)
X_submat = X[seu_sig_gname,]
X_pseudo = X_submat
X_pseudo[X_pseudo == 0] = 0.1
X_normalize = log(X_pseudo)
X_pc = prcomp(t(X_normalize), scale = TRUE)$x[,1:50]
umap_df = umap(d = X_pc, random_state = 123, min_dist = 0.3,  spread = 1, n_neighbors = 10)
labs = factor(celltype_labels, levels = ordered_ct)
df <- data.frame(x = umap_df$layout[,1],
                 cluster = factor(labs),
                 y = umap_df$layout[,2])



pdf(paste0(set.indx,"_vst_umap1.pdf"), height = 4, width = 4)
p1 = ggplot(df, aes(x, y, colour = cluster)) + geom_point(size = 0.3, alpha = 0.5) + 
  theme_bw()  + scale_color_d3("category20") + 
  theme(axis.title=element_blank(), 
        axis.text = element_blank(), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "none",
        panel.grid.minor = element_line(size = 0.5), 
        panel.border = element_rect(fill=NA, colour = "black", size=2),
        panel.grid.major = element_line(size = 1))

p2 = p1 + guides(color = guide_legend(override.aes = list(size = 5)))
print(p2)
dev.off()


pdf(paste0(set.indx,"_vst_umap2.pdf"), height = 4, width = 4)
p1 = ggplot(df, aes(x, y, colour = cluster)) + geom_point(size = 0.3, alpha = 0.5) + 
  xlim(-15,15) + ylim(-15,8) +  
  theme_bw()  + scale_color_d3("category20") + 
  theme(axis.title=element_blank(), 
        axis.text = element_blank(), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "none",
        panel.grid.minor = element_line(size = 0.5), 
        panel.border = element_rect(fill=NA, colour = "black", size=2),
        panel.grid.major = element_line(size = 1))

p2 = p1 + guides(color = guide_legend(override.aes = list(size = 5)))
print(p2)
dev.off()




X_submat = X[selected_feature,]
X_pseudo = X_submat
X_pseudo[X_pseudo == 0] = 0.1
X_normalize = log(X_pseudo)
X_pc = prcomp(t(X_normalize), scale = TRUE)$x[,1:50]
umap_df = umap(d = X_pc, random_state = 123, min_dist = 0.3,  spread = 1, n_neighbors = 10)
labs = factor(celltype_labels, levels = ordered_ct)
df <- data.frame(x = umap_df$layout[,1],
                 cluster = factor(labs),
                 y = umap_df$layout[,2])



pdf(paste0(set.indx,"_wds_umap1.pdf"), height = 4, width = 4)
p1 = ggplot(df, aes(x, y, colour = cluster)) + geom_point(size = 0.3, alpha = 0.5) + 
  theme_bw()  + scale_color_d3("category20") + 
  theme(axis.title=element_blank(), 
        axis.text = element_blank(), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "none",
        panel.grid.minor = element_line(size = 0.5), 
        panel.border = element_rect(fill=NA, colour = "black", size=2),
        panel.grid.major = element_line(size = 1))

p2 = p1 + guides(color = guide_legend(override.aes = list(size = 5)))
print(p2)
dev.off()


pdf(paste0(set.indx,"_wds_umap2.pdf"), height = 4, width = 4)
p1 = ggplot(df, aes(x, y, colour = cluster)) + geom_point(size = 0.3, alpha = 0.5) + 
  xlim(-15,15) + ylim(-15,8) +  
  theme_bw()  + scale_color_d3("category20") + 
  theme(axis.title=element_blank(), 
        axis.text = element_blank(), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "none",
        panel.grid.minor = element_line(size = 0.5), 
        panel.border = element_rect(fill=NA, colour = "black", size=2),
        panel.grid.major = element_line(size = 1))

p2 = p1 + guides(color = guide_legend(override.aes = list(size = 5)))
print(p2)
dev.off()





