# Description:


## tmp: revised from cas_hrvatin2

# ------------------------------------------------------------------------------
#                            Libraries 
# ------------------------------------------------------------------------------


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
library(RColorBrewer)



### set.seed(10050805) # for random mislabel cells

# ------------------------------------------------------------------------------
#                        Data & Benchmark Labels 
# ------------------------------------------------------------------------------

#### Data source:

## https://github.com/mohuangx/SAVER-paper
## under the SAVER-data folder, it provides the google drive link to datasets:
## https://drive.google.com/drive/folders/1uos3DDymyrh3ZyxveDhqQIcOyD1brxeL?usp=sharing

## ngene = ; ncell = ;

size_factor_vec = SizeFactor(X)
cur_info = data.frame(cellname = colnames(X), main_type = maintype_labels, cell_type = celltype_labels)
saveRDS(cur_info, "cortex_cinfo.RDS")

start_time = Sys.time()	
maintype_return = CalculateCDI(sub_gcmat = X.sub, 
	                             candidate_label_dataframe = cur_info$main_type, 
	                             batch_label = NULL, 
	                             cell_size_factor = size_factor_vec, 
	                             ncore = 30)
end_time = Sys.time()
difftime(end_time, start_time)

celltype_return =  CalculateCDI(sub_gcmat = X.sub, 
	                             candidate_label_dataframe = cur_info$cell_type, 
	                             batch_label = NULL, 
	                             cell_size_factor = size_factor_vec, 
	                             ncore = 30)

saveRDS(list(benchmark_main = maintype_return, 
						 benchmark_cell = celltype_return), 
	"cortex_benchmark_return.rds")

benchmark_return = readRDS("~/hrvatin/cortex_benchmark_return.rds")





# ------------------------------------------------------------------------------
#                        Lineplots 
# ------------------------------------------------------------------------------


aicbic_df = readRDS("~/hrvatin/hrX_aicbic_df.RData")
cdi_df = aicbic_df %>% 
	mutate(AIC = AIC/(1e+06), BIC = BIC/(1e+06)) %>%
	rename("Cluster_method" =  "method", "N_cluster" = "k", 
		"CDI_AIC" = "AIC" , "CDI_BIC" = "BIC")
saveRDS(cdi_df, "~/hrvatin/cortex_cdi_df.rds")




pdf("cortex_aic_lineplot.pdf",  width = 3, height = 3.5)
CDILineplot(cdi_df, "CDI_AIC", 
	benchmark_celltype_ncluster = length(unique(cur_info$cell_type)), 
	benchmark_celltype_cdi = lapply(benchmark_return[["benchmark_cell"]], function(x) x/10^6), 
	benchmark_maintype_ncluster = length(unique(cur_info$main_type)), 
	benchmark_maintype_cdi = lapply(benchmark_return[["benchmark_main"]], function(x) x/10^6), 
	show_axis_names = FALSE, show_method_legend = FALSE) +
	scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
		rgb(136, 146, 179, max = 255), #light purpleblue -- HC
    rgb(112, 186, 211, max = 255), #lightblue -- K-means
    rgb(229, 160, 133, max = 255),#pink -- SC3
    "#BCBD22B2",  #-- Seurat
    rgb(76, 158, 137, max = 255))) #green -- Spectral
dev.off()






pdf("cortex_bic_lineplot.pdf",  width = 3, height = 3.5)
CDILineplot(cdi_df, "CDI_BIC", 
	benchmark_celltype_ncluster = length(unique(cur_info$cell_type)), 
	benchmark_celltype_cdi = lapply(benchmark_return[["benchmark_cell"]], function(x) x/10^6), 
	benchmark_maintype_ncluster = length(unique(cur_info$main_type)), 
	benchmark_maintype_cdi = lapply(benchmark_return[["benchmark_main"]], function(x) x/10^6), 
	show_axis_names = FALSE, show_method_legend = FALSE) +
	scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
		rgb(136, 146, 179, max = 255), #light purpleblue -- HC
    rgb(112, 186, 211, max = 255), #lightblue -- K-means
    rgb(229, 160, 133, max = 255),#pink -- SC3
    "#BCBD22B2",  #-- Seurat
    rgb(76, 158, 137, max = 255))) #green -- Spectral
dev.off()



set.indx = "hrX"
ordered_ct = c("Endo_1","Endo_2", "SM_1", "SM_2", # Endo 4
               "Astro",
               "ExcL23", "ExcL4", "ExcL5_1", "ExcL5_2", "ExcL5_3", "ExcL6", "Hip", "RSP", "Sub", #Excl 9
               "Int_Cck", "Int_Npy", "Int_Pv", "Int_Sst_1", "Int_Sst_2", "Int_Vip", # Int 6
               "Micro_1", "Micro_2", 
               "Macrophage",
               "Pericyte", # Mural
               "Olig_1", "Olig_2", "Olig_3", "Olig_4", "Olig_5", "Olig_6", "Olig_7", "OPC_1", "OPC_2" #Olig 9
)
maintype_labels = ifelse(maintype_labels == "Endothelial_SmoothMuscle", "SmoothMuscle", maintype_labels)

ordered_ct2 = c("Endothelial", "Astrocytes", "Excitatory", "Interneurons",
                "Microglia", "Macrophages", "Mural", "Oligodendrocytes")
pdf(paste0(set.indx, "_bicopt_contingency.pdf"), width = 3, height = 6)
labcomp_list = LabComp(factor(maintype_labels, levels = ordered_ct2), 
                       as.character(readRDS(paste0(set.indx, "_Seurat_k20_labs.rds"))[,1]))
print(labcomp_list[[1]])
dev.off()



pdf(paste0("cortex_aicopt_contingency.pdf"), width = 4, height = 8)
labcomp_list = LabComp(factor(celltype_labels, levels = ordered_ct), 
	as.character(readRDS(paste0(set.indx, "_Seurat_k36_labs.rds"))[,1])) 
p = labcomp_list[[1]] +
	scale_x_discrete(labels=c("1", "", "",  "", "5", "", "", "", "", "10", 
		"", "", "", "",  "15", "", "",  "", "", "20", "", "", "", "",  "25", 
		"", "",  "", "", "30", "",  "",  "", "", "35",  ""))
print(p)
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


ctname = data.frame(ct = maintype_labels)
ctname = ctname %>% mutate(ct = recode(ct, "Endothelial_SmoothMuscle" = "Endothelial")) 

pdf(paste0("cortex_bicopt_contingency.pdf"), width = 4, height = 2.5)
labcomp_list = LabComp(factor(ctname$ct, levels = ordered_ct2), 
                       as.character(readRDS(paste0(set.indx, "_Seurat_k20_labs.rds"))[,1]))
p = labcomp_list[[1]] +
	scale_x_discrete(labels=c("1", "", "",  "", "5", "", "", "", "", "10", 
		"", "", "", "",  "15", "", "",  "", "", "20"))
print(p)
dev.off()



## for legend

library(dplyr)
dat <- data.frame(x = 1:9,
                  y = 1:9,
	Method = c("CIDR", "HC", "KMeans", "SC3", "Seurat", 
		"Spectral", "Benchmark main type", 
		"Benchmark cell type", "CDI optimal"))


pdf("/Users/jfang/Documents/tmp_combine_figs/lineplot_legend.pdf", width = 5, height = 5)
ggplot(dat, aes(x = x,y = y, color = Method, shape = Method)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("CIDR" = rgb(64, 83, 133, max = 255), 
  	"HC" = rgb(136, 146, 179, max = 255), 
  	"KMeans" = rgb(112, 186, 211, max = 255), 
  	"SC3" = rgb(229, 160, 133, max = 255), 
  	"Seurat" = "#BCBD22B2", 
  	"Spectral" = rgb(76, 158, 137, max = 255), 
    "Benchmark main type" = "purple3",
  	"Benchmark cell type" = "#7E6148FF", 
  	"CDI optimal" = rgb(210,85,62, max = 255))) +
  scale_shape_manual(values = c("CIDR" = 19, 
  	"HC" = 19, "KMeans" = 19, "SC3" = 19, 
  	"Seurat" = 19, "Spectral" = 19, 
  	"Benchmark main type" = 18, 
  	"Benchmark cell type" = 18,
  	"CDI optimal" = 17)) +
  theme_classic() +  theme(legend.text = element_text(size=15), 
  	legend.title = element_text(size=20))
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
  	# scale_x_discrete(labels=c(1:length(unique(longData$givlab)))) +
    labs(fill = "Proportion") +  
    theme_classic() + theme(axis.text = element_text(size=15, angle=0, vjust=0.3), 
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank(),
                            axis.line = element_line(size = 1), 
                            #axis.title=element_text(size=17), 
                            #legend.position = "none",
    												legend.text = element_text(size = 15),
    												legend.title = element_text(size = 20))
  return(list(p, order_level))
}

a = c(rep("a",20), rep("b", 30))
b = c(rep("1",20), rep("2", 30))

pdf("/Users/jfang/Documents/tmp_combine_figs/contingency_legend.pdf", width = 5, height = 5)
LabComp(a,b)
dev.off()






