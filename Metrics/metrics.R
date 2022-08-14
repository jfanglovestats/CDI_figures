library(NMI) # NMI
library(mclust) # ARI
library(dendextend) # Fowlkes-Mallows Index
library(clusterCrit)
library(cluster) # silhouette
library(ggplot2)
library(ggsci)
library(clValid) # Connectivity
library(doParallel)


## NMI
metric_nmi = function(true_lab, lab_cluster){
  nc = length(true_lab)
  df1 = data.frame(indx = seq_len(nc), true = true_lab)
  df2 = data.frame(indx = seq_len(nc), true = lab_cluster)
  return(as.numeric(NMI(df1, df2)))
}

# FM
metric_fm = function(true_lab, lab_cluster){
  fm_res = FM_index(true_lab, lab_cluster, assume_sorted_vectors = TRUE)
  return(as.numeric(fm_res[1]))
}


get_sup_metrics = function(set_prefix, true_lab, ncl_vec, mname, rds_title){
  nk = length(ncl_vec)
  nm = length(mname)
  sup_metrics = data.frame(
    k = rep(ncl_vec, nm), 
    method = rep(mname, rep(nk, nm)))
  for(m in 1:nm){
    for(k in 1:nk){
      lab_df = readRDS(paste0(set_prefix, "_", mname[m], "_k", ncl_vec[k],"_labs.rds"))
      int_lab = as.integer(as.factor(lab_df[,1]))
      sup_metrics[((m-1)*nk + k), "ARI"] = adjustedRandIndex(true_lab, int_lab)
      sup_metrics[((m-1)*nk + k), "FM"] = metric_fm(true_lab, int_lab)
      sup_metrics[((m-1)*nk + k), "NMI"] = metric_nmi(true_lab, int_lab)
    }
  }
  saveRDS(sup_metrics, rds_title)
  return(sup_metrics)
}



get_uns_metrics = function(set_prefix, tX, ncl_vec, mname, rds_title, cdi_df, cdi_type){
  nk = length(ncl_vec)
  nm = length(mname)
  # DataDist = dist(tX)
  uns_metrics = data.frame(
    k = rep(ncl_vec, nm), 
    method = rep(mname, rep(nk, nm)))
  all_m_res = data.frame()
  for(m in 1:nm){
    cur_m_res <- foreach(k = seq_len(nk), .combine = 'rbind') %dopar% {
      tmp_vec = c()
      lab_df = readRDS(paste0(set_prefix, "_", mname[m], "_k", ncl_vec[k],"_labs.rds"))
      int_lab = as.integer(as.factor(lab_df[,1]))
      row_indx = which((cdi_df$Cluster_method == mname[m]) & (cdi_df$N_cluster == ncl_vec[k]))
      tmp_vec["CDI"] = cdi_df[row_indx, cdi_type]
      tmp_vec["CH"] = as.numeric(intCriteria(tX, int_lab, "Calinski_Harabasz"))
      tmp_vec["Connectivity"] = as.numeric(connectivity(Data = tX, clusters = int_lab))
      tmp_vec["DB"] = as.numeric(intCriteria(tX, int_lab, "Davies_Bouldin"))
      tmp_vec["Dunn"] = as.numeric(intCriteria(tX, int_lab, "Dunn"))
      tmp_vec["Gamma"] = as.numeric(intCriteria(tX, int_lab, "Gamma"))
      tmp_vec["SD_Scat"] = as.numeric(intCriteria(tX, int_lab, "SD_Scat"))
      tmp_vec["Silhouette"] = as.numeric(intCriteria(tX, int_lab, "Silhouette"))
      tmp_vec["XB"] = as.numeric(intCriteria(tX, int_lab, "Xie_Beni"))
      return(tmp_vec)
    }
    all_m_res <- rbind(all_m_res, cur_m_res)
  }
  uns_metrics <- cbind(uns_metrics, all_m_res)
  saveRDS(uns_metrics, rds_title)
  return(uns_metrics)
}

corr_figure = function(sup_df, uns_df, rds_title){
  ex_m = colnames(sup_df)[-c(1,2)]
  in_m = colnames(uns_df)[-c(1,2)]
  nex = length(ex_m)
  nin = length(in_m)
  sp_cor = data.frame(
    ex_metrics = rep(ex_m, rep(nin, nex)), 
    in_metrics = rep(in_m, nex),
    corr = numeric(nex * nin))
  for(e in seq_len(nex)){
    for(i in seq_len(nin)){
      sp_cor_test = suppressWarnings(cor.test(sup_df[,e+2], uns_df[,i+2], method = "spearman"))
      sp_cor[((e-1)*nin + i), "corr"] = as.numeric(sp_cor_test$estimate)
    }
  }
  saveRDS(sp_cor, rds_title)
  p1 <- ggplot(sp_cor, aes(fill = as.factor(in_metrics), y = abs(corr), x= as.factor(ex_metrics))) +
    geom_bar(position="dodge", stat="identity") + 
    theme_classic() + scale_fill_d3() + 
    labs(x = "External metrics", y = "|Spearman correlation|", fill = "Internal metrics")
  return(p1)
}

# ------------------------------------------------------------------------------
#                           SD1
# ------------------------------------------------------------------------------


set.indx = "sd1"
para.set = readRDS(paste0(set.indx, "_paraset.rds"))
ncl_vec = c(4, 6, 7, 8,  9, 10, 11, 12, 14, 16)


get_sup_metrics(
  set_prefix = set.indx,
  true_lab = para.set$clust.label.real, 
  ncl_vec = ncl_vec, 
  mname = c("KMeans", "HC", "Spectral", "CIDR", "SC3", "Seurat"),
  rds_title = paste0(set.indx, "_sup_metrics.rds"))

cl <- makeCluster(length(ncl_vec), type="FORK")
registerDoParallel(cl)
get_uns_metrics(
  set_prefix = set.indx,
  tX = t(readRDS(paste0(set.indx, "_X_sub.rds"))), 
  ncl_vec = ncl_vec, 
  mname = c("KMeans", "HC", "Spectral", "CIDR", "SC3", "Seurat"),
  cdi_df = readRDS(paste0(set.indx, "_cdi_df.rds")), 
  cdi_type = "CDI_BIC",
  rds_title = paste0(set.indx, "_uns_metrics.rds"))
stopCluster(cl)



### Spearman

pdf("sd1_metrics_corr.pdf", height = 2, width = 6)
print(corr_figure(sup_df = readRDS("sd1_sup_metrics.rds"), 
                  uns_df = readRDS("sd1_uns_metrics.rds"), 
                  rds_title = "sd1_metrics_corr.rds"))
dev.off()




# ------------------------------------------------------------------------------
#                                  SD2 (rare)
# ------------------------------------------------------------------------------
set.indx = "sd2"
setwd(paste0("/hpc/home/jf243/", set.indx))
para.set = readRDS(paste0(set.indx, "_paraset.rds"))
lab_prefix = "rare9"
ncl_vec = c(2:8)


sup_metrics <- get_sup_metrics(
  set_prefix = lab_prefix,
  true_lab = para.set$clust.label.real, 
  ncl_vec = ncl_vec, 
  mname = c("KMeans", "HC", "Spectral", "CIDR", "SC3", "Seurat"),
  rds_title = paste0(set.indx, "_sup_metrics.rds"))

cl <- makeCluster(length(ncl_vec), type="FORK")
registerDoParallel(cl)
uns_metrics <- get_uns_metrics(
  set_prefix = lab_prefix,
  tX = t(readRDS(paste0(set.indx, "_X_sub.rds"))), 
  ncl_vec = ncl_vec, 
  mname = c("KMeans", "HC", "Spectral", "CIDR", "SC3", "Seurat"),
  cdi_df = readRDS(paste0(set.indx, "_cdi_df.rds")), 
  cdi_type = "CDI_BIC",
  rds_title = paste0(set.indx, "_uns_metrics.rds"))
stopCluster(cl)

uns_df = readRDS("sd2_uns_metrics.rds")
uns_df_sub <- uns_df %>% select(-c("S_Dbw"))

pdf("sd2_metrics_corr.pdf", height = 2, width = 6)
print(corr_figure(sup_df = readRDS("sd2_sup_metrics.rds"), 
                  uns_df = uns_df_sub, 
                  rds_title = "sd2_metrics_corr.rds"))
dev.off()


# ------------------------------------------------------------------------------
#                                  SD3 (mainsub)
# ------------------------------------------------------------------------------

set.indx = "sd3"
setwd(paste0("/hpc/home/jf243/", set.indx))
para.set = readRDS(paste0(set.indx, "_paraset.rds"))
lab_prefix = "mainsub5"
ncl_vec = c(2:8)

## celltype
sup_metrics <- get_sup_metrics(
  set_prefix = lab_prefix,
  true_lab = para.set$celltype, 
  ncl_vec = ncl_vec, 
  mname = c("KMeans", "HC", "Spectral", "CIDR", "SC3", "Seurat"),
  rds_title = paste0(set.indx, "_celltype_sup_metrics.rds"))

cl <- makeCluster(length(ncl_vec), type="FORK")
registerDoParallel(cl)
uns_metrics <- get_uns_metrics(
  set_prefix = lab_prefix,
  tX = t(readRDS(paste0(set.indx, "_X_sub.rds"))), 
  ncl_vec = ncl_vec, 
  mname = c("KMeans", "HC", "Spectral", "CIDR", "SC3", "Seurat"),
  cdi_df = readRDS(paste0(set.indx, "_cdi_df.rds")), 
  cdi_type = "CDI_BIC",
  rds_title = paste0(set.indx, "_celltype_uns_metrics.rds"))
stopCluster(cl)



uns_df = readRDS("sd3_celltype_uns_metrics.rds")
uns_df_sub <- uns_df %>% select(-c("S_Dbw"))

pdf("sd3_celltype_metrics_corr.pdf", height = 2, width = 6)
print(corr_figure(sup_df = readRDS("sd3_celltype_sup_metrics.rds"), 
                  uns_df = uns_df_sub, 
                  rds_title = "sd3_celltype_metrics_corr.rds"))
dev.off()


## subtype 
sup_metrics <- get_sup_metrics(
  set_prefix = lab_prefix,
  true_lab = para.set$subtype, 
  ncl_vec = ncl_vec, 
  mname = c("KMeans", "HC", "Spectral", "CIDR", "SC3", "Seurat"),
  rds_title = paste0(set.indx, "_subtype_sup_metrics.rds"))

cl <- makeCluster(length(ncl_vec), type="FORK")
registerDoParallel(cl)
uns_metrics <- get_uns_metrics(
  set_prefix = lab_prefix,
  tX = t(readRDS(paste0(set.indx, "_X_sub.rds"))), 
  ncl_vec = ncl_vec, 
  mname = c("KMeans", "HC", "Spectral", "CIDR", "SC3", "Seurat"),
  cdi_df = readRDS(paste0(set.indx, "_cdi_df.rds")), 
  cdi_type = "CDI_AIC",
  rds_title = paste0(set.indx, "_subtype_uns_metrics.rds"))
stopCluster(cl)


uns_df = readRDS("sd3_subtype_uns_metrics.rds")
uns_df_sub <- uns_df %>% select(-c("S_Dbw"))

pdf("sd3_subtype_metrics_corr.pdf", height = 2, width = 6)
print(corr_figure(sup_df = readRDS("sd3_subtype_sup_metrics.rds"), 
                  uns_df = uns_df_sub, 
                  rds_title = "sd3_subtype_metrics_corr.rds"))
dev.off()




# ------------------------------------------------------------------------------
#                                  SD4 (splatter)
# ------------------------------------------------------------------------------
set.indx = "sd4"
setwd(paste0("/hpc/home/jf243/", set.indx))
para.set = readRDS(paste0(set.indx, "_paraset.rds"))
lab_prefix = "splat5cl3"
ncl_vec = c(2:8)


sup_metrics <- get_sup_metrics(
  set_prefix = lab_prefix,
  true_lab = para.set$clust.label.real, 
  ncl_vec = ncl_vec, 
  mname = c("KMeans", "HC", "Spectral", "CIDR", "SC3", "Seurat"),
  rds_title = paste0(set.indx, "_sup_metrics.rds"))

cl <- makeCluster(length(ncl_vec), type="FORK")
registerDoParallel(cl)
uns_metrics <- get_uns_metrics(
  set_prefix = lab_prefix,
  tX = t(readRDS(paste0(set.indx, "_X_sub.rds"))), 
  ncl_vec = ncl_vec, 
  mname = c("KMeans", "HC", "Spectral", "CIDR", "SC3", "Seurat"),
  cdi_df = readRDS(paste0(set.indx, "_cdi_df.rds")), 
  cdi_type = "CDI_BIC",
  rds_title = paste0(set.indx, "_uns_metrics.rds"))
stopCluster(cl)

uns_df = readRDS("sd4_uns_metrics.rds")
uns_df_sub <- uns_df %>% select(-c("S_Dbw"))

pdf("sd4_metrics_corr.pdf", height = 2, width = 6)
print(corr_figure(sup_df = readRDS("sd4_sup_metrics.rds"), 
                  uns_df = uns_df_sub, 
                  rds_title = "sd4_metrics_corr.rds"))
dev.off()

# ------------------------------------------------------------------------------
#                                  TCELL 
# ------------------------------------------------------------------------------


set.indx = "tcell"
setwd(paste0("/hpc/home/jf243/", set.indx))
para.set = readRDS(paste0(set.indx, "_paraset.rds"))
lab_prefix = "lsc_size"
ncl_vec = c(2,   3,    4,    5,    6,    7,    8,    10,   12,   14,   16)


sup_metrics <- get_sup_metrics(
  set_prefix = lab_prefix,
  true_lab = para.set$celltype, 
  ncl_vec = ncl_vec, 
  mname = c("KMeans", "HC", "Spectral", "CIDR", "SC3", "Seurat"),
  rds_title = paste0(set.indx, "_sup_metrics.rds"))

cl <- makeCluster(length(ncl_vec), type="FORK")
registerDoParallel(cl)
uns_metrics <- get_uns_metrics(
  set_prefix = lab_prefix,
  tX = t(readRDS(paste0(set.indx, "_X_sub.rds")))*1.0, 
  ncl_vec = ncl_vec, 
  mname = c("KMeans", "HC", "Spectral", "CIDR", "SC3", "Seurat"),
  cdi_df = readRDS(paste0(set.indx, "_cdi_df.rds")), 
  cdi_type = "CDI_BIC",
  rds_title = paste0(set.indx, "_uns_metrics.rds"))
stopCluster(cl)


uns_df = readRDS("tcell_uns_metrics.rds")
uns_df_sub <- uns_df %>% select(-c("S_Dbw"))

pdf("tcell_metrics_corr.pdf", height = 2, width = 6)
print(corr_figure(sup_df = readRDS("tcell_sup_metrics.rds"), 
                  uns_df = uns_df_sub, 
                  rds_title = "tcell_metrics_corr.rds"))
dev.off()


# ------------------------------------------------------------------------------
#                                  Cortex
# ------------------------------------------------------------------------------
# setwd("/Users/jiyuanfang/Documents/CDI_intermediate_results/tmp_metrics/")
corr_figure(sup_df = readRDS("cortex_subtype_sup_metrics.rds"), 
            uns_df = readRDS("cortex_subtype_uns_metrics.rds"), 
            rds_title = "cortex_subtype_metrics_corr.rds")

corr_figure(sup_df = readRDS("cortex_celltype_sup_metrics.rds"), 
            uns_df = readRDS("cortex_celltype_uns_metrics.rds"), 
            rds_title = "cortex_maintype_metrics_corr.rds")


# ------------------------------------------------------------------------------
#                                 Retina 
# ------------------------------------------------------------------------------
# results obtained from retina_cluster.R
k5 = readRDS("/Users/jiyuanfang/Desktop/metrics_d0405/retina_maintype_uns_metrics_k5.rds") 
k10 = readRDS("/Users/jiyuanfang/Desktop/metrics_d0405/retina_maintype_uns_metrics_k10.rds")
k15 = readRDS("/Users/jiyuanfang/Desktop/metrics_d0405/retina_maintype_uns_metrics_k15.rds")
k16 = readRDS("/Users/jiyuanfang/Desktop/metrics_d0405/retina_maintype_uns_metrics_k16.rds") 
k17 = readRDS("/Users/jiyuanfang/Desktop/metrics_d0405/retina_maintype_uns_metrics_k17.rds") 
k18 = readRDS("/Users/jiyuanfang/Desktop/metrics_d0405/retina_maintype_uns_metrics_k18.rds") 
k19 = readRDS("/Users/jiyuanfang/Desktop/metrics_d0405/retina_maintype_uns_metrics_k19.rds") 
k20 = readRDS("/Users/jiyuanfang/Desktop/metrics_d0405/retina_maintype_uns_metrics_k20.rds") 
k25 = readRDS("/Users/jiyuanfang/Desktop/metrics_d0405/retina_maintype_uns_metrics_k25.rds") 
k30 = readRDS("/Users/jiyuanfang/Desktop/metrics_d0405/retina_maintype_uns_metrics_k30.rds") 
k31 = readRDS("/Users/jiyuanfang/Desktop/metrics_d0405/retina_maintype_uns_metrics_k31.rds") 
k32 = readRDS("/Users/jiyuanfang/Desktop/metrics_d0405/retina_maintype_uns_metrics_k32.rds") 
k33 = readRDS("/Users/jiyuanfang/Desktop/metrics_d0405/retina_maintype_uns_metrics_k33.rds") 
k35 = readRDS("/Users/jiyuanfang/Desktop/metrics_d0405/retina_maintype_uns_metrics_k35.rds") 


retina_uns = do.call(rbind, list(k5, k10, k15, k16, k17, k18, k19, 
                                 k20, k25, k30, k31, k32, k33, k35))
retina_uns <- retina_uns %>%
  arrange(factor(method, levels = c("KMeans", "HC", "CIDR", "SC3", "Seurat"))) %>%
  select(c(k, method, CDI_AIC, CDI_BIC, CH, Connectivity, DB, Dunn, Gamma, SD_Scat, Silhouette, XB))
saveRDS(retina_uns, "/Users/jiyuanfang/Documents/CDI_intermediate_results/tmp_metrics/retina_uns_metrics.rds")
retina_sup = readRDS("/Users/jiyuanfang/Documents/CDI_intermediate_results/tmp_metrics/retina_maintype_sup_metrics.rds")
retina_sup <- arrange(retina_sup, factor(method, levels = c("KMeans", "HC", "CIDR", "SC3", "Seurat"))) 


retina_uns = readRDS("/Users/jiyuanfang/Documents/CDI_intermediate_results/tmp_metrics/retina_uns_metrics.rds")
retina_uns_aic <- retina_uns %>% select(-c(CDI_BIC)) %>% rename(CDI = CDI_AIC) 
retina_uns_bic <- retina_uns %>% select(-c(CDI_AIC)) %>% rename(CDI = CDI_BIC) 

corr_figure(sup_df = retina_sup, 
            uns_df = retina_uns_aic, 
            rds_title = "retina_subtype_metrics_corr.rds")

corr_figure(sup_df = retina_sup, 
            uns_df = retina_uns_bic, 
            rds_title = "retina_maintype_metrics_corr.rds")

# ------------------------------------------------------------------------------
#                                  Figure together 
# ------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggsci)
setwd("/Users/jiyuanfang/Documents/CDI_intermediate_results/tmp_metrics")
sd1 = readRDS("sd1_metrics_corr.rds")
sd1 <- sd1 %>% mutate(dataset = "sd1", lab_type = "main")

sd2 = readRDS("sd2_metrics_corr.rds")
sd2 <- sd2 %>% mutate(dataset = "sd2", lab_type = "main")

sd3m = readRDS("sd3_celltype_metrics_corr.rds")
sd3m <- sd3m %>% mutate(dataset = "sd3m", lab_type = "main")

sd3s = readRDS("sd3_subtype_metrics_corr.rds")
sd3s <- sd3s %>% mutate(dataset = "sd3s", lab_type = "sub")

sd4 = readRDS("sd4_metrics_corr.rds")
sd4 <- sd4 %>% mutate(dataset = "sd4", lab_type = "main")

tcell = readRDS("tcell_metrics_corr.rds")
tcell <- tcell %>% mutate(dataset = "tcell", lab_type = "main")

corm = readRDS("cortex_maintype_metrics_corr.rds")
corm <- corm %>% mutate(dataset = "cortexm", lab_type = "main") %>% filter(in_metrics != "S_Dbw")

cors = readRDS("cortex_subtype_metrics_corr.rds")
cors <- cors %>% mutate(dataset = "cortexs", lab_type = "sub") %>% filter(in_metrics != "S_Dbw")

retm = readRDS("retina_maintype_metrics_corr.rds")
retm <- retm %>% mutate(dataset = "retinam", lab_type = "main") 

rets = readRDS("retina_subtype_metrics_corr.rds")
rets <- rets %>% mutate(dataset = "retinas", lab_type = "sub")

df_combined1 = do.call(rbind, list(sd1, sd2, sd3m, sd3s, sd4, tcell, corm, cors, retm, rets))

sd1[sd1$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"] = - sd1[sd1$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"]
sd2[sd2$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"] = - sd2[sd2$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"]
sd3m[sd3m$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"] = - sd3m[sd3m$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"]
sd3s[sd3s$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"] = - sd3s[sd3s$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"]
sd4[sd4$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"] = - sd4[sd4$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"]
tcell[tcell$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"] = - tcell[tcell$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"]
corm[corm$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"] = - corm[corm$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"]
cors[cors$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"] = - cors[cors$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"]
retm[retm$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"] = - retm[retm$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"]
rets[rets$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"] = - rets[rets$in_metrics %in% c("CDI", "Connectivity", "DB", "SD_Scat", "XB"), "corr"]


df_combined = do.call(rbind, list(sd1, sd2, sd3m, sd3s, sd4, tcell, corm, cors, retm, rets))
write.csv(df_combined, "/Users/jiyuanfang/Documents/CDI_review/new_analysis/spearman_correaltion.csv")
# pdf("/Users/jiyuanfang/Documents/CDI_review/new_analysis/metrics_boxplot_original2.pdf", width = 5, height =15)
# ggplot(df_combined, aes(x = factor(in_metrics, levels = c("CDI", "Connectivity", "DB", "SD_Scat", "XB", "CH", "Dunn", "Gamma", "Silhouette")), y = corr, fill = in_metrics)) + 
#   geom_boxplot() + facet_wrap(~ ex_metrics) + 
#   # ylim(0,1.05) + 
#   geom_point(aes(shape = lab_type), alpha = 0.7) + scale_shape_manual(values=c(16, 17)) +
#   facet_grid(rows = vars(ex_metrics)) + 
#   theme_classic() + scale_fill_d3(alpha = 0.6) + 
#   theme(legend.position="none") + 
#   labs(x = "Internal metrics", 
#        y = "Spearman correlation", 
#        fill = "Internal metrics")
# dev.off()

aa = df_combined %>% filter((in_metrics == "Dunn") & (ex_metrics == "FM"))
quantile(aa$corr, 0.75) - quantile(aa$corr, 0.25)

aa = df_combined %>% filter((in_metrics == "CDI") & (ex_metrics == "FM"))
quantile(aa$corr, 0.75) - quantile(aa$corr, 0.25)

aa = df_combined %>% filter((in_metrics == "CH") & (ex_metrics == "FM"))
quantile(aa$corr, 0.75) - quantile(aa$corr, 0.25)

aa = df_combined %>% filter((in_metrics == "DB") & (ex_metrics == "FM"))
quantile(aa$corr, 0.75) - quantile(aa$corr, 0.25)

aa = df_combined %>% filter((in_metrics == "Dunn") & (ex_metrics == "ARI"))
quantile(aa$corr, 0.75) - quantile(aa$corr, 0.25)

aa = df_combined %>% filter((in_metrics == "CDI") & (ex_metrics == "ARI"))
quantile(aa$corr, 0.75) - quantile(aa$corr, 0.25)

aa = df_combined %>% filter((in_metrics == "CH") & (ex_metrics == "ARI"))
quantile(aa$corr, 0.75) - quantile(aa$corr, 0.25)

aa = df_combined %>% filter((in_metrics == "DB") & (ex_metrics == "ARI"))
quantile(aa$corr, 0.75) - quantile(aa$corr, 0.25)


ggplot(df_combined1, aes(x = in_metrics, y = corr, fill = in_metrics)) + 
  geom_boxplot() + facet_wrap(~ ex_metrics) + 
  # ylim(0,1.05) + 
  geom_point(aes(shape = lab_type), alpha = 0.7) + scale_shape_manual(values=c(16, 17)) +
  facet_grid(rows = vars(ex_metrics)) + 
  theme_classic() + scale_fill_d3(alpha = 0.6) + 
  theme(legend.position="none") + 
  labs(x = "Internal metrics", 
       y = "Absolute value of Spearman correlation", 
       fill = "Internal metrics")



pdf("/Users/jiyuanfang/Documents/CDI_review/results_for_review/metrics_boxplot.pdf", width = 6, height = 7)
ggplot(df_combined, aes(x = in_metrics, y = corr, fill = in_metrics)) + 
  geom_boxplot() + facet_wrap(~ ex_metrics) + 
  # ylim(0,1.05) + 
  geom_point(aes(shape = lab_type), alpha = 0.7) + scale_shape_manual(values=c(16, 17)) +
  facet_grid(rows = vars(ex_metrics)) + 
  theme_classic() + scale_fill_d3(alpha = 0.6) + 
  theme(legend.position="none") + 
  labs(x = "Internal metrics", 
       y = "Spearman correlation", 
       fill = "Internal metrics")
dev.off()



# ------------------------------------------------------------------------------
#                                 Optimal label analysis
# ------------------------------------------------------------------------------

opt_info = function(uns_df, sup_df, data_name, benchmark_type){
  int = c("CDI", "CH", "Connectivity", "DB", "Dunn", "Gamma", "SD_Scat", "Silhouette", "XB")
  if(!("CDI" %in% colnames(uns_df))){
    if(benchmark_type == "main"){
      colnames(uns_df) = ifelse(colnames(uns_df) == "CDI_BIC", "CDI", colnames(uns_df))
    } else{
      colnames(uns_df) = ifelse(colnames(uns_df) == "CDI_AIC", "CDI", colnames(uns_df))
    }
  }
  uns_df[, c("CDI", "Connectivity", "DB", "SD_Scat", "XB")] = - uns_df[, c("CDI", "Connectivity", "DB", "SD_Scat", "XB")]
  n_int = length(int)
  ext = c("ARI", "FM", "NMI")
  n_ext = length(ext)
  df = data.frame(dataset = rep(data_name, n_int * n_ext), 
                  type = rep(benchmark_type, n_int * n_ext),
                  int = rep(int, rep(n_ext, n_int)), 
                  k = numeric(n_int * n_ext), 
                  method = rep(NA, n_int * n_ext), 
                  ext = rep(ext, n_int), 
                  ext_val = numeric(n_int * n_ext))
  for(int_m in int){
    uns_target = uns_df[which.max(uns_df[,int_m]), ]
    df[(df$int == int_m), "k"] = uns_target$k
    df[(df$int == int_m), "method"] = uns_target$method
    sup_target = sup_df[(sup_df$k == uns_target$k) & (sup_df$method == uns_target$method), ]
    df[(df$int == int_m) & (df$ext == "ARI"), "ext_val"] = sup_target$ARI
    df[(df$int == int_m) & (df$ext == "FM"), "ext_val"] = sup_target$FM
    df[(df$int == int_m) & (df$ext == "NMI"), "ext_val"] = sup_target$NMI
  }
  return(df)
}






sd1_opt = opt_info(uns_df = readRDS("sd1_uns_metrics.rds"), 
                   sup_df = readRDS("sd1_sup_metrics.rds"), 
                   data_name = "SD1", 
                   benchmark_type = "main")
sd2_opt = opt_info(uns_df = readRDS("sd2_uns_metrics.rds"), 
                   sup_df = readRDS("sd2_sup_metrics.rds"), 
                   data_name = "SD2", 
                   benchmark_type = "main")
sd3_opt_m = opt_info(uns_df = readRDS("sd3_celltype_uns_metrics.rds"), 
                     sup_df = readRDS("sd3_celltype_sup_metrics.rds"), 
                     data_name = "SD3", 
                     benchmark_type = "main")
sd3_opt_s = opt_info(uns_df = readRDS("sd3_subtype_uns_metrics.rds"), 
                     sup_df = readRDS("sd3_subtype_sup_metrics.rds"), 
                     data_name = "SD3", 
                     benchmark_type = "subtype")
sd4_opt = opt_info(uns_df = readRDS("sd4_uns_metrics.rds"), 
                   sup_df = readRDS("sd4_sup_metrics.rds"), 
                   data_name = "SD4", 
                   benchmark_type = "main")
tcell_opt = opt_info(uns_df = readRDS("tcell_uns_metrics.rds"), 
                     sup_df = readRDS("tcell_sup_metrics.rds"), 
                     data_name = "T-CELL", 
                     benchmark_type = "main")
cortex_opt_m = opt_info(uns_df = readRDS("cortex_celltype_uns_metrics.rds"), 
                        sup_df = readRDS("cortex_celltype_sup_metrics.rds"), 
                        data_name = "CORTEX", 
                        benchmark_type = "main")
cortex_opt_s = opt_info(uns_df = readRDS("cortex_subtype_uns_metrics.rds"), 
                        sup_df = readRDS("cortex_subtype_sup_metrics.rds"), 
                        data_name = "CORTEX", 
                        benchmark_type = "subtype")
retina_opt_m = opt_info(uns_df = readRDS("retina_uns_metrics.rds"), 
                        sup_df = readRDS("retina_maintype_sup_metrics.rds"), 
                        data_name = "RETINA", 
                        benchmark_type = "main")
retina_opt_s = opt_info(uns_df = readRDS("retina_uns_metrics.rds"), 
                        sup_df = readRDS("retina_subtype_sup_metrics.rds"), 
                        data_name = "RETINA", 
                        benchmark_type = "subtype")

opt_df =do.call(rbind, list(sd1_opt, sd2_opt, sd3_opt_m, sd3_opt_s, sd4_opt, 
                            tcell_opt, cortex_opt_m, cortex_opt_s,
                            retina_opt_m, retina_opt_s))



opt_df2 <- opt_df %>% 
  mutate(title = factor(ifelse(opt_df$type == "main",  as.character(opt_df$dataset), paste0(opt_df$dataset, "(", opt_df$type, ")")), 
                        levels = c("SD1", "SD2", "SD3", "SD3(subtype)", "SD4", "T-CELL", "CORTEX", "CORTEX(subtype)", "RETINA", "RETINA(subtype)"))) %>%
  filter(ext == "ARI") %>% mutate(int_title = factor(int, levels = c("Benchmark", as.character(unique(opt_df$int))))) %>% select(title, int_title, k) 


opt_df3 = rbind(opt_df2, 
                data.frame(title = c("SD1", "SD2", "SD3", "SD3(subtype)", "SD4", "T-CELL", "CORTEX", "CORTEX(subtype)", "RETINA", "RETINA(subtype)"), 
                           int_title = rep("Benchmark", length(unique(opt_df2$title))), 
                           k = c(10, 4, 2, 4, 5, 5, 8, 33, 6, 18)))


pdf("/Users/jiyuanfang/Documents/CDI_review/results_for_review/opt_ncluster_barplot.pdf", width = 5, height = 14)
ggplot(opt_df3, aes(x = int_title, y = k, fill = int_title)) + 
  geom_bar(stat="identity") + facet_wrap(~ title, ncol = 3) +  theme_classic() + 
  labs(y = "Number of clusters", fill = "Benchmark/Internal metrics") + 
  theme(
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    # legend.position = "none",
    legend.justification = c("right", "bottom"),
    legend.box.just = "left",
    legend.margin = margin(t = 0, 0, b = 0, l = -200),
    legend.spacing.x = unit(0.1, 'cm'),
    legend.spacing.y = unit(0.1, 'cm')) +
  scale_fill_manual(values = c("black", pal_d3()(9))) + 
  guides(fill=guide_legend(ncol=2,nrow=5,byrow=TRUE)) 

dev.off()




pdf("/Users/jiyuanfang/Documents/CDI_review/results_for_review/opt_metrics_boxplot.pdf", width = 6, height = 7)
ggplot(opt_df, aes(x = int, y = ext_val, fill = int)) +
  geom_boxplot() + facet_wrap(~ ext) +
  # ylim(0,1.05) +
  geom_point(aes(shape = type), alpha = 0.7) + scale_shape_manual(values=c(16, 17)) +
  facet_grid(rows = vars(ext)) +
  theme_classic() + scale_fill_d3(alpha = 0.6) +
  theme(legend.position="none") +
  labs(x = "Internal metrics",
       y = "Optimal candidate label external metric value",
       fill = "Internal metrics")
dev.off()

## Cortex

cortex_s_ari = subset(cortex_opt_s, ext == "ARI")
cortex_s_ari[order(-cortex_s_ari$ext_val),]

cortex_s_fm = subset(cortex_opt_s, ext == "FM")
cortex_s_fm[order(-cortex_s_fm$ext_val),]

cortex_s_nmi = subset(cortex_opt_s, ext == "NMI")
cortex_s_nmi[order(-cortex_s_nmi$ext_val),]



cortex_m_ari = subset(cortex_opt_m, ext == "ARI")
cortex_m_ari[order(-cortex_m_ari$ext_val),]

cortex_m_fm = subset(cortex_opt_m, ext == "FM")
cortex_m_fm[order(-cortex_m_fm$ext_val),]

cortex_m_nmi = subset(cortex_opt_m, ext == "NMI")
cortex_m_nmi[order(-cortex_m_nmi$ext_val),]



## Retina

retina_s_ari = subset(retina_opt_s, ext == "ARI")
retina_s_ari[order(-retina_s_ari$ext_val),]

retina_s_fm = subset(retina_opt_s, ext == "FM")
retina_s_fm[order(-retina_s_fm$ext_val),]

retina_s_nmi = subset(retina_opt_s, ext == "NMI")
retina_s_nmi[order(-retina_s_nmi$ext_val),]



retina_m_ari = subset(retina_opt_m, ext == "ARI")
retina_m_ari[order(-retina_m_ari$ext_val),]

retina_m_fm = subset(retina_opt_m, ext == "FM")
retina_m_fm[order(-retina_m_fm$ext_val),]

retina_m_nmi = subset(retina_opt_m, ext == "NMI")
retina_m_nmi[order(-retina_m_nmi$ext_val),]


## t-cell

tcell_ari = subset(tcell_opt, ext == "ARI")
tcell_ari[order(-tcell_ari$ext_val),]

tcell_fm = subset(tcell_opt, ext == "FM")
tcell_fm[order(-tcell_fm$ext_val),]

tcell_nmi = subset(tcell_opt, ext == "NMI")
tcell_nmi[order(-tcell_nmi$ext_val),]











