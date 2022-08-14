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
      tmp_vec["S_Dbw"] = as.numeric(intCriteria(tX, int_lab, "S_Dbw"))
      tmp_vec["SD_Scat"] = as.numeric(intCriteria(tX, int_lab, "SD_Scat"))
      # tmp_vec["Silhouette"] = summary(silhouette(int_lab, DataDist))$avg.width
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






set.indx = "cortex"
setwd("CORTEX")
para.set = readRDS(paste0(set.indx, "_paraset.rds"))
lab_prefix = "cortex"
ncl_vec = c( 7:10,   15,     19:21,    25,    30:37,   40)


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



# sup_metrics <- get_sup_metrics(
#   set_prefix = lab_prefix,
#   true_lab = para.set$subtype, 
#   ncl_vec = ncl_vec, 
#   mname = c("KMeans", "HC", "Spectral", "CIDR", "SC3", "Seurat"),
#   rds_title = paste0(set.indx, "_subtype_sup_metrics.rds"))
# 
# cl <- makeCluster(length(ncl_vec), type="FORK")
# registerDoParallel(cl)
# uns_metrics <- get_uns_metrics(
#   set_prefix = lab_prefix,
#   tX = t(readRDS(paste0(set.indx, "_X_sub.rds"))), 
#   ncl_vec = ncl_vec, 
#   mname = c("KMeans", "HC", "Spectral", "CIDR", "SC3", "Seurat"),
#   cdi_df = readRDS(paste0(set.indx, "_cdi_df.rds")), 
#   cdi_type = "CDI_AIC",
#   rds_title = paste0(set.indx, "_subtype_uns_metrics.rds"))
# stopCluster(cl)


