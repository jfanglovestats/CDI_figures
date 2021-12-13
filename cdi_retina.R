#### Description:








# ------------------------------------------------------------------------------
#                            Libraries 
# ------------------------------------------------------------------------------
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
library(CDI)
source("utils.R")



# ------------------------------------------------------------------------------
#                        Data & Benchmark Labels & Preprocess
# ------------------------------------------------------------------------------

#### Data source

## Data download
# The RData format dataset was given in 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81904 
# (GSE81904_bipolar_data_Cell2016.Rdata.gz)
#
## Preprocessing & benchmark labels
# Follow the preprocessing steps in 
# https://github.com/broadinstitute/BipolarCell2016/blob/master/BCanalysis.pdf,
# we got a filtered count matrix with 13,166 genes and 27,499 cells as shown in 
# page 2 of BCanalysis.pdf.
#
# These 27,499 cells were clustered by the authors to 26 clusters (page 6 of BCanalysis.pdf), 
# Among these 26 clusters, 18 clusters (26,830 cells) were identified with biomarkers to 18 cell types. 
# The cell type names and corresponding markers are available in the Table S1 of the original paper. 
# Our analysis starts from the count matrix with 26,830 cells from 18 cell types.
# 
## Main type labels
# Based on the paper, we also group the cell types to 5 main types
# OFF cone bipolar cells: BC1A, BC1B, BC2, BC3A, BC3B, BC4
# On cone bipolar cells: BC5A, BC5B, BC5C, BC5D, BC6, BC7, BC8_BC9
# Photo receptors: Rod photo receptors, Cone photo receptors
# Muller glia: Muller glia
# Amacrine: Amacrine
# Rod bipolar cells: Rod bipolar cells


scmat = readRDS("bipolar.RData")


#... from cas_retina2.R get qc2_mat

X = as.matrix(qc2_mat)
info = readRDS("bipolar_cell_info2.rds")



# ------------------------------------------------------------------------------
#                       Generate Candidate Labels
# ------------------------------------------------------------------------------

### Batch effect correction with Seurat
Seurat.obj = CreateSeuratObject(counts = X, project = set.indx, min.cells = 0, min.features = 0)
Seurat.obj = NormalizeData(Seurat.obj)
Seurat.obj = FindVariableFeatures(Seurat.obj, selection.method = "vst", nfeatures = 2000)
Seurat.obj = ScaleData(Seurat.obj, verbose = FALSE)
Seurat.obj = RunPCA(Seurat.obj, npcs = 30, verbose = FALSE)
Seurat.obj = RunUMAP(Seurat.obj, reduction = "pca", dims = 1:20)
Seurat.obj@meta.data["batch"] = para.set$Batch
Seurat.obj.list <- SplitObject(Seurat.obj, split.by = "batch")
Seurat.obj.list <- lapply(X = Seurat.obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
Seurat.obj.anchors <- FindIntegrationAnchors(object.list = Seurat.obj.list, 
                                             dims = 1:30, 
                                             anchor.features = 2000)
Seurat.obj.integrated <- IntegrateData(anchorset = Seurat.obj.anchors, 
                                       features.to.integrate = Seurat.obj.anchors@anchor.features, 
                                       dims = 1:30)
saveRDS(Seurat.obj.integrated, paste0(set.indx, "Seurat_integrated_obj.rds"))
X_bc_unsorted = as.matrix(Seurat.obj.integrated[["RNA"]]@data)
col_indx = colnames(X)
X_bc = X_bc_unsorted[, col_indx]
saveRDS(X_bc, "retina_bc_countmat.rds")


#### Generate candidate label sets with batch effect corrected count


## Seurat
Seurat.obj.integrated = readRDS(paste0(set.indx, "Seurat_integrated_obj.rds"))
ncl_vec = c(5,     10,   15,     16,    17,     18,    19,    20,     25,     30)
seu_res = c(0.0035,0.01, 0.1,    0.18,   0.3,    0.5,    0.6,  0.71,   1.2,    2)
seurat_lab_bc(Seurat.obj.integrated, seu_res, colnames(X))








# ------------------------------------------------------------------------------
#                           Feature Selection 
# ------------------------------------------------------------------------------

## select genes 

feature_gene_indx = FeatureGeneSelection(gcmat = X, 
                                         batch_label = info$Batch, 
                                         method = "wds", 
                                         nfeature = 500)
sub_X = X[feature_gene_indx,]





# ------------------------------------------------------------------------------
#                      UMAPS from Selected Features
# ------------------------------------------------------------------------------








# ------------------------------------------------------------------------------
#                           Calculate CDI
# ------------------------------------------------------------------------------

## calculate size factor
size_factor_vec = SizeFactor(X)

## combine all labels as a data frame
cluster_number = c(5,     10,   15,     16,    17,     18,    19,    20,     25,     30,  31, 32, 33, 35)
ncluster = length(cluster_number)
method_name = c("KMeans", "HC", "CIDR", "SC3", "Seurat")
nmethod = length(method_name)
lab_df = matrix(NA, nrow = para.set$ncell, ncol = ncluster*nmethod)
colnames(lab_df) = paste0("tmp_", c(1:(ncluster*nmethod)))




# ... update with new package functions

# saveRDS(cdi_df, "cdi_df.rds")



	
start_time = Sys.time()	
maintype_return = CalculateCDI(sub_gcmat = sub_X, 
	                             candidate_label_dataframe = cur_info$main_type, 
	                             batch_label = cur_info$Batch, 
	                             cell_size_factor = size_factor_vec, 
	                             ncore = 30)
end_time = Sys.time()
difftime(end_time, start_time)

celltype_return =  CalculateCDI(sub_gcmat = sub_X, 
	                             candidate_label_dataframe = cur_info$cell_type, 
	                             batch_label = cur_info$Batch, 
	                             cell_size_factor = size_factor_vec, 
	                             ncore = 30)

saveRDS(list(benchmark_main = maintype_return, benchmark_cell = celltype_return), "retina_benchmark_return.rds")

benchmark_return = readRDS("~/bipolar/retina_benchmark_return.rds")


# ------------------------------------------------------------------------------
#                        Results Demonstration
# ------------------------------------------------------------------------------


## CDI Lineplots
cdi_df = readRDS("d0822_cdi_df.rds")
cdi_df = cdi_df %>% 
	rename("Cluster_method" =  "Method", "N_cluster" = "Cluster", 
		"CDI_AIC" = "HighResolutionCDI" , "CDI_BIC" = "LowResolutionCDI") %>%
	mutate(CDI_AIC = CDI_AIC/(1e+06), CDI_BIC = CDI_BIC/(1e+06))
saveRDS(cdi_df, "~/bipolar/retina_cdi_df.rds")

cdi_df =  readRDS( "~/bipolar/retina_cdi_df.rds")


		

cur_info = cur_info %>% 
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



                
                

pdf("retina_aic_lineplot.pdf", width = 3, height = 3.5)
CDILineplot(cdi_df, "CDI_AIC",
	benchmark_celltype_ncluster = length(unique(cur_info$cell_type)), 
	benchmark_celltype_cdi = lapply(benchmark_return[["benchmark_cell"]], function(x) x/10^6), 
	benchmark_maintype_ncluster = length(unique(cur_info$main_type)), 
	benchmark_maintype_cdi = lapply(benchmark_return[["benchmark_main"]], function(x) x/10^6), 
	show_axis_names = FALSE, 
	show_method_legend = FALSE) +
	scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
		rgb(136, 146, 179, max = 255), #light purpleblue -- HC
    rgb(112, 186, 211, max = 255), #lightblue -- K-means
    rgb(229, 160, 133, max = 255),#pink -- SC3
    "#BCBD22B2")) #-- Seurat
dev.off()


pdf("retina_bic_lineplot.pdf", width = 3, height = 3.5)
CDILineplot(cdi_df, "CDI_BIC",
	benchmark_celltype_ncluster = length(unique(cur_info$cell_type)), 
	benchmark_celltype_cdi = lapply(benchmark_return[["benchmark_cell"]], function(x) x/10^6), 
	benchmark_maintype_ncluster = length(unique(cur_info$main_type)), 
	benchmark_maintype_cdi = lapply(benchmark_return[["benchmark_main"]], function(x) x/10^6), 
	show_axis_names = FALSE, 
	show_method_legend = FALSE) +
	scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
		rgb(136, 146, 179, max = 255), #light purpleblue -- HC
    rgb(112, 186, 211, max = 255), #lightblue -- K-means
    rgb(229, 160, 133, max = 255),#pink -- SC3
    "#BCBD22B2")) #-- Seurat
dev.off()



## Contingency table



ctname = ctname %>% mutate(ct = recode(ct, "Endothelial_SmoothMuscle" = "Endothelial")) 



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
  	# scale_x_discrete(labels=c(1:length(unique(longData$givlab)))) +
    labs(fill = "Proportion") +  
    theme_classic() + theme(axis.text = element_text(size=15, angle=0, vjust=0.3), 
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank(),
                            axis.line = element_line(size = 1), 
                            #axis.title=element_text(size=17), 
                            legend.position = "none")
  return(list(p, order_level))
}




ordered_maintype = c("OFF BC", "ON BC", "Rod BC", "Amacrine", "Muller glia", "Photoreceptor")
pdf(paste0("retina_bicopt_contingency.pdf"), width = 4, height = 2.5)
labcomp_list = LabComp(factor(cur_info$main_type, levels = ordered_maintype), 
                       as.character(readRDS(paste0(set.indx, "_Seurat_k18_labs.rds"))[,1]))
p = labcomp_list[[1]] +
	scale_x_discrete(labels=c("1", "", "",  "", "5", "", "", "", "", "10", 
		"", "", "", "",  "15", "", "",  ""))
print(p)
dev.off()


ordered_celltype = c("BC1A", "BC1B", "BC2", "BC3A", "BC3B", "BC4", "BC5A", "BC5B", "BC5C", 
                     "BC5D", 'BC6', 'BC7', "BC8_BC9",  "Rod BC", "Amacrine", "Muller glia", 
                     "Cone PR", "Rod PR")

pdf(paste0("retina_aicopt_contingency.pdf"), width = 4, height = 8)
labcomp_list = LabComp(factor(cur_info$cell_type, levels = ordered_celltype), 
                       as.character(readRDS(paste0(set.indx, "_KMeans_k33_labs.rds"))[,1]))
p = labcomp_list[[1]] +
	scale_x_discrete(labels=c("1", "", "",  "", "5", "", "", "", "", "10", 
		"", "", "", "",  "15", "", "",  "", "", "20", "", "", "", "",  "25", 
		"", "",  "", "", "30", "",  "",  ""))
print(p)
dev.off()







pdf(paste0(set.indx, "_maintype_lowres_contingency.pdf"), width = 9, height = 6)
CandidateBenchmarkHeatmap(benchmark_label = factor(cur_info$main_type, levels = ordered_maintype),  
                          candidate_label =  as.character(readRDS(paste0(set.indx, "_Seurat_k18_labs.rds"))[,1]),  
                          proportion_size = 4,
                          show_axis_names = FALSE,
                          show_legend = FALSE)
dev.off()


# ------------------------------------------------------------------------------
#                      Metric comparison
# ------------------------------------------------------------------------------


## CDI vs ARI
for(i in 1:ncluster){
  for(m in 1:nmethod){
    label_indx = (i-1)*nmethod + m
    cluster_label = lab_df[,label_indx]
    cdi_df[label_indx, "Maintype_ARI"] = adjustedRandIndex(cur_info$main_type, cluster_label)
    cdi_df[label_indx, "Celltype_ARI"] = adjustedRandIndex(cur_info$cell_type, cluster_label)
  }
  cat("progress = ", round(i/ncluster,2), "\n")
}




pdf(paste0(set.indx, "_highres_vs_celltype_ari.pdf"), width = 4, height = 4)
ggplot(data = cdi_df, aes(x = HighResolutionCDI/1000000, y = Celltype_ARI)) + 
  theme_bw() + 
  geom_point(alpha = 0.5) + ylim(-0.5,1.1) + 
  stat_smooth(method = "gam", size = 2, color = "red") + 
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=15), 
        axis.text.x = element_text(size=15),
        panel.grid.minor = element_line(size = 0.5), 
        panel.border = element_rect(fill=NA, colour = "black", size=2),
        panel.grid.major = element_line(size = 1))
dev.off()




pdf(paste0(set.indx, "_highres_vs_maintype_ari.pdf"), width = 4, height = 4)
ggplot(data = cdi_df, aes(x = HighResolutionCDI/1000000, y = Maintype_ARI)) + 
  theme_bw() + 
  geom_point(alpha = 0.5) + ylim(-0.5,1.1) + 
  stat_smooth(method = "gam", size = 2, color = "red") + 
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=15), 
        axis.text.x = element_text(size=15),
        panel.grid.minor = element_line(size = 0.5), 
        panel.border = element_rect(fill=NA, colour = "black", size=2),
        panel.grid.major = element_line(size = 1))
dev.off()



pdf(paste0(set.indx, "_lowres_vs_celltype_ari.pdf"), width = 4, height = 4)
ggplot(data = cdi_df, aes(x = LowResolutionCDI/1000000, y = Celltype_ARI)) + 
  theme_bw() + 
  geom_point(alpha = 0.5) + ylim(-0.5,1.1) + 
  stat_smooth(method = "gam", size = 2, color = "blue") + 
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=15), 
        axis.text.x = element_text(size=15),
        panel.grid.minor = element_line(size = 0.5), 
        panel.border = element_rect(fill=NA, colour = "black", size=2),
        panel.grid.major = element_line(size = 1))
dev.off()

pdf(paste0(set.indx, "_lowres_vs_maintype_ari.pdf"), width = 4, height = 4)
ggplot(data = cdi_df, aes(x = LowResolutionCDI/1000000, y = Maintype_ARI)) + 
  theme_bw() + 
  geom_point(alpha = 0.5) + ylim(-0.5,1.1) + 
  stat_smooth(method = "gam", size = 2, color = "blue") + 
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=15), 
        axis.text.x = element_text(size=15),
        panel.grid.minor = element_line(size = 0.5), 
        panel.border = element_rect(fill=NA, colour = "black", size=2),
        panel.grid.major = element_line(size = 1))
dev.off()

