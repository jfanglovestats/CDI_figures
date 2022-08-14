
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
library(dplyr)


set.indx = "cortex"
X = readRDS(paste0(set.indx, "_filtered.rds"))
para.set = readRDS(paste0(set.indx, "_paraset.rds"))


# ------------------------------------------------------------------------------
#                       Generate Candidate Labels
# ------------------------------------------------------------------------------
source("utils/clustering.R")
ncl_vec = c( 7,     8,     9,    10,   15,     19,     20,      21,    25,     30,   31,   32,   33,   34,   35,   36,   37,  40)


## Seurat
seu_res = c(0.005,  0.01,  0.015, 0.02, 0.15, 0.6,   0.71,     0.75,   1.3,    1.9,  2.1,  2.19,  2.3, 2.35, 2.5,  2.58,  2.7,  3.5)
df = data.frame(seu_res, ncluster = seurat_lab(X, seu_res))



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
  Seurat_obj = NULL,
  cand_lab_df = lab_df, 
  cell_size_factor = X_sc, 
  batch_label = NULL,
  BPPARAM = MulticoreParam(para.set$ncore))
saveRDS(cdi_df, paste0(set.indx, "_cdi_df.rds"))

start_time = Sys.time()	
celltype_return = calculate_CDI(
  sub_gcmat = X_sub, 
  cand_lab_df = para.set$celltype, 
  batch_label = NULL,
  BPPARAM = MulticoreParam(para.set$ncore))

subtype_return = calculate_CDI(
  sub_gcmat = X_sub, 
  cand_lab_df = para.set$subtype, 
  batch_label = NULL,
  BPPARAM = MulticoreParam(para.set$ncore))
end_time = Sys.time()
difftime(end_time, start_time)


saveRDS(list(benchmark_cell = celltype_return, 
             benchmark_sub = subtype_return), 
        paste0(set.indx, "_benchmark_return.rds"))



# ------------------------------------------------------------------------------
#                            Figures (lineplot & Contingency) 
# ------------------------------------------------------------------------------


cdi_df <- cdi_df %>% mutate(CDI_AIC = CDI_AIC / 10^6, CDI_BIC = CDI_BIC / 10^6)



pdf(paste0(set.indx, "_aic_lineplot.pdf"),  width = 3, height = 3.5)
CDI_lineplot(
  cdi_dataframe = cdi_df, 
  cdi_type = "CDI_AIC",
  benchmark_maintype_ncluster = length(unique(para.set$celltype)),
  benchmark_maintype_cdi = lapply(benchmark_return[["benchmark_cell"]], function(x) x/10^6),
	benchmark_celltype_ncluster = length(unique(para.set$subtype)),
	benchmark_celltype_cdi = lapply(benchmark_return[["benchmark_sub"]], function(x) x/10^6)) +
	scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
		rgb(136, 146, 179, max = 255), #light purpleblue -- HC
    rgb(112, 186, 211, max = 255), #lightblue -- K-means
    rgb(229, 160, 133, max = 255),#pink -- SC3
    "#BCBD22B2",  #-- Seurat
    rgb(76, 158, 137, max = 255))) #green -- Spectral
dev.off()




pdf(paste0(set.indx, "_aic_lineplot.pdf"), width = 3, height = 3.5)
CDI_lineplot(
  cdi_dataframe = cdi_df, 
  cdi_type = "CDI_BIC",
  benchmark_maintype_ncluster = length(unique(para.set$celltype)),
  benchmark_maintype_cdi = lapply(benchmark_return[["benchmark_cell"]], function(x) x/10^6),
  benchmark_celltype_ncluster = length(unique(para.set$subtype)),
  benchmark_celltype_cdi = lapply(benchmark_return[["benchmark_sub"]], function(x) x/10^6)) +
  scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
                                rgb(136, 146, 179, max = 255), #light purpleblue -- HC
                                rgb(112, 186, 211, max = 255), #lightblue -- K-means
                                rgb(229, 160, 133, max = 255),#pink -- SC3
                                "#BCBD22B2",  #-- Seurat
                                rgb(76, 158, 137, max = 255))) #green -- Spectral
dev.off()







maintype_labels = para.set$celltype
maintype_labels = ifelse(maintype_labels == "Endothelial_SmoothMuscle", "Endothelial", maintype_labels)

ordered_ct2 = c("Endothelial", "Astrocytes", "Excitatory", "Interneurons",
                "Microglia", "Macrophages", "Mural", "Oligodendrocytes")

pdf(paste0(set.indx, "_bicopt_contingency.pdf"), width = 3, height = 6)
contingency_heatmap(factor(maintype_labels, levels = ordered_ct2),
                    as.character(readRDS(paste0("./clustering_labels/", set.indx, "_Seurat_k20_labs.rds"))[,1]))
dev.off()


celltype_labels = para.set$subtype
ordered_ct = c("Endo_1","Endo_2", "SM_1", "SM_2", # Endo 4
               "Astro",
               "ExcL23", "ExcL4", "ExcL5_1", "ExcL5_2", "ExcL5_3", "ExcL6", "Hip", "RSP", "Sub", #Excl 9
               "Int_Cck", "Int_Npy", "Int_Pv", "Int_Sst_1", "Int_Sst_2", "Int_Vip", # Int 6
               "Micro_1", "Micro_2",
               "Macrophage",
               "Pericyte", # Mural
               "Olig_1", "Olig_2", "Olig_3", "Olig_4", "Olig_5", "Olig_6", "Olig_7", "OPC_1", "OPC_2" #Olig 9
)

pdf(paste0("cortex_aicopt_contingency.pdf"), width = 4, height = 8)
p1 <- contingency_heatmap(factor(celltype_labels, levels = ordered_ct), 
                    as.character(readRDS(paste0("./clustering_labels/", set.indx, "_Seurat_k36_labs.rds"))[,1]))
p = p1 + scale_x_discrete(labels=c("1", "", "",  "", "5", "", "", "", "", "10",
		"", "", "", "",  "15", "", "",  "", "", "20", "", "", "", "",  "25",
		"", "",  "", "", "30", "",  "",  "", "", "35",  ""))
print(p)
dev.off()



# ------------------------------------------------------------------------------
#                           UMAP
# ------------------------------------------------------------------------------


ordered_ct = c("Endo_1","Endo_2", "SM_1", "SM_2", # Endo 4
               "Astro",
               "ExcL23", "ExcL4", "ExcL5_1", "ExcL5_2", "ExcL5_3", "ExcL6", "Hip", "RSP", "Sub", #Excl 9
               "Int_Cck", "Int_Npy", "Int_Pv", "Int_Sst_1", "Int_Sst_2", "Int_Vip", # Int 6
               "Micro_1", "Micro_2", 
               "Macrophage",
               "Pericyte", # Mural
               "Olig_1", "Olig_2", "Olig_3", "Olig_4", "Olig_5", "Olig_6", "Olig_7", "OPC_1", "OPC_2" #Olig 9
)


col_panel = c("#CCCCCC", "#969696", "#636363", "#252525", # grey
              "gold2",
              "#FCBBA1", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000", #OrRd
              "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F", "#3F007D", #Purple
              "hotpink4", "violetred3",
              "salmon4",
              "wheat3",
              "#E0F3DB", "#CCEBC5", "#A8DDB5", "#7BCCC4", "#9ECAE1","#4EB3D3", "#2B8CBE", "#0868AC", "#084081" #GnBu
)



## Seurat - VST
Seurat.obj = CreateSeuratObject(counts = X, project = "hh", min.cells = 0, min.features = 0)
Seurat.obj = NormalizeData(Seurat.obj, verbose = FALSE)
Seurat.obj = FindVariableFeatures(Seurat.obj, selection.method = "vst", nfeatures = para.set$nfeature, verbose = FALSE)
seu_sig_gname = Seurat.obj@assays[["RNA"]]@var.features

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
df$y = ifelse(df$y> 20, df$y-9, df$y)


pdf(paste0(set.indx,"_seurat_umap.pdf"), height = 4, width = 4)
p1 = ggplot(df, aes(x, y, colour = cluster)) + geom_point(size = 0.3, alpha = 0.5) + 
  theme_bw()  + scale_color_manual(values = col_panel) + 
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



## CDI - WDS
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
df$y = ifelse(df$y> 20, df$y-9, df$y)



pdf(paste0(set.indx,"_cdi_umap.pdf"), height = 4, width = 4)
p1 = ggplot(df, aes(x, y, colour = cluster)) + geom_point(size = 0.3, alpha = 0.5) + 
  theme_bw()  + scale_color_manual(values = col_panel) + 
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



