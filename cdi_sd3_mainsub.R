#from cas_sim_mainsub5.R
library(BiocParallel)


aicbic_df = readRDS("~/simulation/mainsub5/mainsub5_aicbic_df.rds")
cdi_df = aicbic_df %>% 
	mutate(AIC = AIC/(1e+06), BIC = BIC/(1e+06)) %>%
	rename("Cluster_method" =  "method", "N_cluster" = "k", 
		"CDI_AIC" = "AIC" , "CDI_BIC" = "BIC")
saveRDS(cdi_df, "~/simulation/mainsub5/sd3_cdi_df.rds")


X.sub = X[FeatureGeneSelection(X),]


size_factor_vec = SizeFactor(X)
	

start_time = Sys.time()	
maintype_return = CalculateCDI(sub_gcmat = X.sub, 
	                             candidate_label_dataframe = para.set$maintype, 
	                             batch_label = NULL, 
	                             cell_size_factor = size_factor_vec, 
	                             ncore = 30)
end_time = Sys.time()
difftime(end_time, start_time)

celltype_return =  CalculateCDI(sub_gcmat = X.sub, 
	                             candidate_label_dataframe = para.set$clust.label.real, 
	                             batch_label = NULL, 
	                             cell_size_factor = size_factor_vec, 
	                             ncore = 30)

saveRDS(list(benchmark_main = maintype_return, 
						 benchmark_cell = celltype_return), 
	"sd3_benchmark_return.rds")

benchmark_return = readRDS("~/simulation/mainsub5/sd3_benchmark_return.rds")






#------- contingency




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
#                        Lineplots 
# ------------------------------------------------------------------------------


aicbic_df = readRDS("~/hrvatin/hrX_aicbic_df.RData")
cdi_df = aicbic_df %>% 
	mutate(AIC = AIC/(1e+06), BIC = BIC/(1e+06)) %>%
	rename("Cluster_method" =  "method", "N_cluster" = "k", 
		"CDI_AIC" = "AIC" , "CDI_BIC" = "BIC")
saveRDS(cdi_df, "~/hrvatin/cortex_cdi_df.rds")




pdf("sd3_aic_lineplot.pdf",  width = 3, height = 3.5)
CDILineplot(cdi_df, "CDI_AIC", 
	benchmark_celltype_ncluster = length(unique(para.set$clust.label.real)), 
	benchmark_celltype_cdi = lapply(benchmark_return[["benchmark_cell"]], function(x) x/10^6), 
	benchmark_maintype_ncluster = length(unique(para.set$maintype)), 
	benchmark_maintype_cdi = lapply(benchmark_return[["benchmark_main"]], function(x) x/10^6), 
	show_axis_names = FALSE, show_method_legend = FALSE) +
	scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
		rgb(136, 146, 179, max = 255), #light purpleblue -- HC
    rgb(112, 186, 211, max = 255), #lightblue -- K-means
    rgb(229, 160, 133, max = 255),#pink -- SC3
    "#BCBD22B2",  #-- Seurat
    rgb(76, 158, 137, max = 255))) #green -- Spectral
dev.off()


pdf("sd3_bic_lineplot.pdf",  width = 3, height = 3.5)
CDILineplot(cdi_df, "CDI_BIC", 
	benchmark_celltype_ncluster = length(unique(para.set$clust.label.real)), 
	benchmark_celltype_cdi = lapply(benchmark_return[["benchmark_cell"]], function(x) x/10^6), 
	benchmark_maintype_ncluster = length(unique(para.set$maintype)), 
	benchmark_maintype_cdi = lapply(benchmark_return[["benchmark_main"]], function(x) x/10^6), 
	show_axis_names = FALSE, show_method_legend = FALSE) +
	scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
		rgb(136, 146, 179, max = 255), #light purpleblue -- HC
    rgb(112, 186, 211, max = 255), #lightblue -- K-means
    rgb(229, 160, 133, max = 255),#pink -- SC3
    "#BCBD22B2",  #-- Seurat
    rgb(76, 158, 137, max = 255))) #green -- Spectral
dev.off()




pdf("sd3_bicopt_contingency.pdf", width = 2.5, height = 2.5)
lc_return = LabComp(as.character(c(rep(1, 1000), rep(2, 1800))), 
                    as.character(readRDS(paste0(set.indx, "_CIDR_k2_labs.rds"))[,1]))
print(lc_return[[1]])
dev.off()

pdf("sd3_aicopt_contingency.pdf", width = 2.5, height = 2.5)
lc_return = LabComp(as.character(para.set$clust.label.real), 
                    as.character(readRDS(paste0(set.indx, "_Seurat_k4_labs.rds"))[,1]))
print(lc_return[[1]])
dev.off()




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





