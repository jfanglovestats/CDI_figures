# cas_lsc123.R



# ------------------------------------------------------------------------------
#                            Figures (lineplot & Contingency) 
# ------------------------------------------------------------------------------


aicbic_df = readRDS("~/lsc12356/lsc_sizeaicbic_df.RData")
# b = readRDS("~/lsc12356/lsc_sizeaicbic_size_df2.RData")
cdi_df = aicbic_df %>% 
	mutate(AIC = AIC/(1e+06), BIC = BIC/(1e+06)) %>%
	rename("Cluster_method" =  "method", "N_cluster" = "k", 
		"CDI_AIC" = "AIC" , "CDI_BIC" = "BIC")
saveRDS(cdi_df, "~/lsc12356/tcell_cdi_df.rds")

cdi_df[which.min(cdi_df$CDI_BIC),]


size_factor_vec = SizeFactor(X)
	
start_time = Sys.time()	
maintype_return = CalculateCDI(sub_gcmat = X.sub, 
	                             candidate_label_dataframe = cinfo$Group, 
	                             batch_label = NULL, 
	                             cell_size_factor = size_factor_vec, 
	                             ncore = 30)
end_time = Sys.time()
difftime(end_time, start_time)


saveRDS(maintype_return, "tcell_benchmark_return.rds")
benchmark_return = readRDS("tcell_benchmark_return.rds")






cur_info = cinfo %>% rename("main_type" = "Group")

	



pdf("tcell_aic_lineplot.pdf", width = 3, height = 3.5)
CDILineplot(cdi_df, "CDI_AIC",
	benchmark_maintype_ncluster = length(unique(cur_info$main_type)), 
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


pdf("tcell_bic_lineplot.pdf", width = 3, height = 3.5)
CDILineplot(cdi_df, "CDI_BIC",
	benchmark_maintype_ncluster = length(unique(cur_info$main_type)), 
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




ctname = data.frame(ct = para.set$Celltype)
ctname = ctname %>% mutate(ct = recode(ct, 
	'Active EM-like Treg' = 'Active EM-like Treg', 
	'CD8 Tcm, IL17RA+ & CD28+' = 'CD8 Tcm',
	'CD8 Trm cells' = 'CD8 Trm',
	'Classical CD4 Tem' = 'Classical CD4 Tem',
	'regulatory Trm cells' = 'Regulatory Trm'))

pdf(paste0(set.indx, "_bicopt_contingency.pdf"), width = 4, height = 2.5)
LabComp(ctname$ct, readRDS(paste0(set.indx, "_SC3_k5_labs.rds"))[,1])
dev.off()





