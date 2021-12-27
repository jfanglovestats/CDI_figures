
library(ggplot2)
library(dplyr)

CDILineplot = function(cdi_dataframe, 
                       cdi_type, 
                       benchmark_celltype_cdi = NULL, 
                       benchmark_celltype_ncluster = NULL, 
                       benchmark_maintype_cdi = NULL, 
                       benchmark_maintype_ncluster = NULL, 
                       clustering_method = NULL,
                       show_axis_names = TRUE, 
                       show_method_legend = TRUE){
  min_indx = which.min(cdi_dataframe[, cdi_type])
  if(!is.null(clustering_method)){
    cdi_dataframe["Cluster_method"] = clustering_method
  } else if(sum(is.na(cdi_dataframe$Cluster_method))){
    cdi_dataframe[is.na(cdi_dataframe$Cluster_method), "Cluster_method"] = "Unknown_method"
  } 
  p1 = ggplot(cdi_dataframe, aes_string(x = "N_cluster", y = cdi_type, group = "Cluster_method")) + 
    geom_point(aes_string(color = "Cluster_method"), size = 4, alpha = 0.7) + 
    geom_line(aes_string(color = "Cluster_method"), size = 2, alpha = 0.7) + 
    geom_point(data = cdi_dataframe[min_indx, ], 
               colour = grDevices::rgb(210,85,62, max = 255), 
               shape = 17, size = 4) + labs(x = "Number of clusters") + 
    # scale_color_d3() + 
  	theme_bw() + 
  	theme(axis.text = element_text(size = 15),
					panel.grid.minor = element_line(size = 0.5), 
					panel.border = element_rect(fill=NA, colour = "black", size=2),
					panel.grid.major = element_line(size = 1)) 
  if(show_axis_names == FALSE){
    p1 = p1 + theme(axis.title.x = element_blank(), 
                    axis.title.y = element_blank())
  }
  if(show_method_legend == FALSE){
    p1 = p1 + theme(legend.position = "none")
  }
  if((missing(benchmark_celltype_cdi) | missing(benchmark_celltype_ncluster)) & (missing(benchmark_maintype_cdi) | missing(benchmark_maintype_ncluster))){
    return(p1)
  } 
  if(cdi_type == "CDI_AIC"){
    if(!is.null(benchmark_celltype_cdi)){
      p1 = p1 + geom_point(data = data.frame(N_cluster = benchmark_celltype_ncluster, 
                                             CDI_AIC = benchmark_celltype_cdi$CDI_AIC, 
                                             Cluster_method = "Benchmark"), 
                           colour="#7E6148FF", 
                           shape = "*", 
                           size = 15) 
    }
    if(!is.null(benchmark_maintype_cdi)){
      p1 = p1 + geom_point(data = data.frame(N_cluster = benchmark_maintype_ncluster, 
                                             CDI_AIC =  benchmark_maintype_cdi$CDI_AIC, 
                                             Cluster_method = "Benchmark"), 
                           colour="purple3", 
                           shape = "*", 
                           size = 15)
    }
    
  } 
  if(cdi_type == "CDI_BIC"){
    if(!is.null(benchmark_celltype_cdi)){
      p1 = p1 + geom_point(data = data.frame(N_cluster = benchmark_celltype_ncluster, 
                                             CDI_BIC = benchmark_celltype_cdi$CDI_BIC, 
                                             Cluster_method = "Benchmark"), 
                           colour="#7E6148FF",
                           shape = "*", 
                           size = 15)
    }
    if(!is.null(benchmark_maintype_cdi)){
      p1 = p1 + geom_point(data = data.frame(N_cluster = benchmark_maintype_ncluster, 
                                             CDI_BIC = benchmark_maintype_cdi$CDI_BIC, 
                                             Cluster_method = "Benchmark"), 
                           colour= "purple3", 
                           shape = "*", 
                           size = 15)
    }
    
  }
  
   return(p1)
}


# ------------------------- hrvatin -------------------------------------------------------
aicbic_df = readRDS("~/hrvatin/hrX_aicbic_df.RData")
cdi_df = aicbic_df %>% 
	mutate(AIC = AIC/(1e+06), BIC = BIC/(1e+06)) %>%
	rename("Cluster_method" =  "method", "N_cluster" = "k", 
		"CDI_AIC" = "AIC" , "CDI_BIC" = "BIC")
saveRDS(cdi_df, "~/hrvatin/cortex_cdi_df.rds")

saveRDS(, "~/hrvatin/benchmark_main_cdi.rds")
saveRDS(, "~/hrvatin/benchmark_cell_cdi.rds")
pdf("cortex_aic_lineplot.pdf",  width = 3, height = 3.5)
CDILineplot(cdi_df, "CDI_AIC", show_axis_names = FALSE, show_method_legend = FALSE) +
	scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
		rgb(136, 146, 179, max = 255), #light purpleblue -- HC
    rgb(112, 186, 211, max = 255), #lightblue -- K-means
    rgb(229, 160, 133, max = 255),#pink -- SC3
    "#BCBD22B2",  #-- Seurat
    rgb(76, 158, 137, max = 255))) #green -- Spectral
dev.off()


pdf("cortex_bic_lineplot.pdf",  width = 3, height = 3.5)
CDILineplot(cdi_df, "CDI_BIC", show_axis_names = FALSE, show_method_legend = FALSE) +
	scale_color_manual(values = c(rgb(64, 83, 133, max = 255),#purpleblue -- CIDR
		rgb(136, 146, 179, max = 255), #light purpleblue -- HC
    rgb(112, 186, 211, max = 255), #lightblue -- K-means
    rgb(229, 160, 133, max = 255),#pink -- SC3
    "#BCBD22B2",  #-- Seurat
    rgb(76, 158, 137, max = 255))) #green -- Spectral
dev.off()



benchmark_main = AICBIC_given_lab_allcell_size(labels = maintype_labels, scmat = X.sub, 
                                         size_vec = allgene_size, ncore = para.set$ncore)

benchmark_cell = AICBIC_given_lab_allcell_size(labels = celltype_labels, scmat = X.sub, 
                                         size_vec = allgene_size, ncore = para.set$ncore)
