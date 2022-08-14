library(Matrix)
library(matrixStats)
library(MASS)
library(cluster)
library(BiocParallel)
library(stringr)
library(clusterCrit)
library(clValid) # Connectivity
library(dplyr)





neg_nb_logsum <- function(input_parameter, x_vec, sc_vec){
  mu <- input_parameter[1]
  r <- input_parameter[2]
  return(-sum(dnbinom(x_vec, size = r, mu = sc_vec * mu, log = TRUE)))
}


## Calculate MLE of NB distribution with size factor 
## per cell-type per gene
nb_size_mle <- function(x_vec, sc_vec){
  if(sum(x_vec) == 0){
    return(list(param = c(0.001, 0.2), 
                value = neg_nb_logsum(c(0.001, 0.2), x_vec, sc_vec)))
  }
  avg <- mean(x_vec); s2 <- var(x_vec)
  init_mu <- avg
  tmp_init_r <- avg^2/(s2 - avg)
  init_r_0 <- ifelse(is.na(tmp_init_r)|(tmp_init_r < 0.1), 0.1, tmp_init_r)
  init_r <- ifelse(init_r_0 > 50, 50, init_r_0)
  nb_est <- nlminb(c(init_mu, init_r), 
                   objective = neg_nb_logsum, 
                   gradient = NULL, 
                   lower = c(1e-6, 1e-6), 
                   upper = c(1e3, 1e6),
                   x_vec = x_vec,
                   sc_vec = sc_vec)
  return(list(param = nb_est$par, value = nb_est$objective))
}



## Calculate the sum of negative likelihood
## all cell-type per gene per batch
single_batch_one_gene_likelihood <- function(gvec, 
                                             candidate_label, 
                                             ncluster, 
                                             cell_size_factor){
  gvec_list <- split(gvec, f = factor(candidate_label))
  sc_list <- split(cell_size_factor, f = factor(candidate_label))
  neg_llk_g <- 0
  for(k in seq_len(ncluster)){
    gvec_ct <- gvec_list[[k]]
    sc_ct <- sc_list[[k]]
    neg_llk_g  <- neg_llk_g  + nb_size_mle(gvec_ct, sc_ct)$value
  }
  return(neg_llk_g)
  
}



## Calculate the sum of negative likelihood
## all cell-types per gene all batches
multi_batch_one_gene_likelihood <- function(gvec,
                                            candidate_label, 
                                            ncluster,
                                            size_ct_list,
                                            batch_ct_list,
                                            lrt_pval_threshold = 0.01){
  gvec_ct_list <- split(gvec, f = factor(candidate_label))
  lrt_test_pval <- numeric(ncluster)
  total_nllk <- 0
  for(k in seq_len(ncluster)){
    cur_ct_gvec <- gvec_ct_list[[k]]
    cur_ct_size <- size_ct_list[[k]]
    cur_ct_batch <- batch_ct_list[[k]]
    cur_ct_nbatch <- length(unique(cur_ct_batch))
    
    ## fit one NB for all batches within this cluster
    common_nb_nllk <- nb_size_mle(cur_ct_gvec, cur_ct_size)$value
    if(cur_ct_nbatch == 1){
      total_nllk <- total_nllk + common_nb_nllk 
      lrt_test_pval[k] <- 1
    } else{
      sep_nllk <- 0
      gvec_ct_batch_list <- split(cur_ct_gvec, f = cur_ct_batch)
      size_ct_batch_list <- split(cur_ct_size, f = cur_ct_batch)
      
      ## fit separate NB for each batch within this cluster
      for(b in seq_len(cur_ct_nbatch)){
        sep_nllk <- sep_nllk + nb_size_mle(gvec_ct_batch_list[[b]], 
                                           size_ct_batch_list[[b]])$value
      }
      ## likelihood ratio test (LRT) to decide which one to choose
      lrt_test_pval_cur_ct <- round(pchisq(2*(common_nb_nllk - sep_nllk), 
                                           2*(cur_ct_nbatch - 1), 
                                           lower.tail = FALSE), 4)
      lrt_test_pval[k] <- lrt_test_pval_cur_ct
      total_nllk <- total_nllk + ifelse(lrt_test_pval_cur_ct < lrt_pval_threshold,  
                                        sep_nllk, common_nb_nllk)
    }
  }
  return(list(NegLLK = total_nllk, nReject = sum(lrt_test_pval < lrt_pval_threshold)))
  
}


## Calculate the CDI values for one label set
calculate_CDI_oneset <- function(sub_gcmat, 
                                 candidate_label, 
                                 batch_label = NULL, 
                                 cell_size_factor, 
                                 BPPARAM,
                                 lrt_pval_threshold = 0.01){
  original_ncluster <- length(unique(candidate_label))
  
  #### One-batch scenario
  if(is.null(batch_label) | (length(unique(batch_label)) == 1)){
    ## filter clusters with small number of cells
    min_ncell <- min(table(candidate_label))
    if(min_ncell < 3){
      sub_indx <- c(seq_len(length(candidate_label)))[(candidate_label %in% names(which(table(candidate_label) > 2)))]
      candidate_label <- candidate_label[sub_indx]
      sub_gcmat <- sub_gcmat[,sub_indx]
      cell_size_factor <- cell_size_factor[sub_indx]
    }
    ng <- nrow(sub_gcmat); nc <- ncol(sub_gcmat)
    ## after filtering, it is possible ncluster < original_cluster
    ## for fair comparison, use original_cluster in penalty
    ncluster <- length(unique(candidate_label))
    
    ## calculating log-likelihood
    if(is.null(rownames(sub_gcmat))){
      rownames(sub_gcmat) <- paste0("g", seq_len(ng))
    }
    sub_gclist <- split(sub_gcmat, f = rownames(sub_gcmat))
    neg_llk_list <- bplapply(sub_gclist, 
                             single_batch_one_gene_likelihood,  
                             candidate_label = candidate_label, 
                             ncluster = ncluster, 
                             cell_size_factor = cell_size_factor, 
                             BPPARAM = BPPARAM)
    neg_llk <- sum(unlist(neg_llk_list))
    npara <- ng * original_ncluster * 2
    
    #### Multi-batch scenario
  } else{
    combine_label <- paste0("ct_", candidate_label, "_b_", batch_label)
    min_combine_ncell <- min(table(combine_label))
    
    ## filter clusters with small number of cells
    if(min_combine_ncell < 3){
      sub_indx <- c(seq_len(length(combine_label)))[(combine_label %in% names(which(table(combine_label) > 2)))]
      candidate_label <- candidate_label[sub_indx]
      batch_label <- batch_label[sub_indx]
      sub_gcmat <- sub_gcmat[,sub_indx]
      cell_size_factor <- cell_size_factor[sub_indx]
    }
    ng <- nrow(sub_gcmat); nc <- ncol(sub_gcmat)
    ## after filtering, it is possible ncluster < original_cluster
    ## for fair comparison, use original_cluster in penalty
    ncluster <- length(unique(candidate_label))
    if(is.null(rownames(sub_gcmat))){
      rownames(sub_gcmat) <- paste0("g", seq_len(ng))
    }
    batch_ct_list <- split(batch_label, f = candidate_label)
    size_ct_list <- split(cell_size_factor, f = candidate_label)
    # bp = BiocParallel::MulticoreParam(ncore)
    sub_gclist <- split(sub_gcmat, f = rownames(sub_gcmat))
    neg_llk_list <- bplapply(sub_gclist, 
                             multi_batch_one_gene_likelihood,  
                             candidate_label = candidate_label, 
                             ncluster = ncluster, 
                             batch_ct_list = batch_ct_list, 
                             size_ct_list  = size_ct_list,
                             lrt_pval_threshold = lrt_pval_threshold,
                             BPPARAM = BPPARAM)
    neg_llk <- sum(unlist(lapply(neg_llk_list, '[[', 'NegLLK')))
    total_rej <- sum(unlist(lapply(neg_llk_list, '[[', 'nReject')))
    npara <- (ng * original_ncluster + total_rej) * 2
  }
  return(list(CDI_AIC = 2*neg_llk + 2*npara,
              CDI_BIC = 2*neg_llk + npara*log(nc),
              neg_llk_val = neg_llk))
  
}



calculate_CDI <- function(
  sub_gcmat = NULL,
  Seurat_obj = NULL,
  cand_lab_df, 
  cell_size_factor, 
  batch_label = NULL,
  lrt_pval_threshold = 0.01,
  clustering_method = NULL,
  ncore = 1, ...){
  # Initialize a BiocParallel object
  BPPARAM = MulticoreParam(ncore, ...)
  if(is.null(sub_gcmat)){
    sub_gcmat = as.matrix(Seurat_obj@assays$RNA@counts)
  }
  if((is.null(sub_gcmat)) & (is.null(Seurat_obj))){
    stop("One of sub_gcmat and Seurat_obj needs to be non-empty!")
  } 
  ## check format of arguments
  if((!is.vector(cand_lab_df)) & (!is.data.frame(cand_lab_df))){
    stop("Please input a vector or a data frame for cand_lab_df.")
  }
  if(!is.matrix(sub_gcmat)){
    stop("Please input a matrix for sub_gcmat.")
  }
  if(!is.vector(cell_size_factor)){
    stop("Please input a numerical vector for cell_size_factor.")
  }
  if((!is.null(batch_label)) & (!is.vector(batch_label))){
    stop("Please input a vector for batch_label.")
  }
  if((!is.null(clustering_method)) & (!is.vector(clustering_method))){
    stop("Please input a vector for batch_label.")
  }
  ## if cand_lab_df is a vector or a data frame with one column
  vec_or_1col_df = ifelse(is.vector(cand_lab_df), TRUE, dim(cand_lab_df)[2] == 1)
  if(vec_or_1col_df){
    return(calculate_CDI_oneset(sub_gcmat = sub_gcmat, 
                                candidate_label = unlist(cand_lab_df), 
                                batch_label = batch_label, 
                                cell_size_factor = cell_size_factor, 
                                BPPARAM = BPPARAM, 
                                lrt_pval_threshold = lrt_pval_threshold))
    
    ## if cand_lab_df is a a data frame with more than one column
  } else {
    lab_name <- colnames(cand_lab_df)
    cdi_return_df <- data.frame(Label_name = paste0("Label", seq_len(ncol(cand_lab_df))))
    if(!is.null(lab_name)){
      cdi_return_df["Label_name"] <- lab_name
      cdi_return_df["Cluster_method"] <- ifelse(
        grepl(pattern = "^(\\w+)(_k)(\\d+)$", x = lab_name, ignore.case = TRUE), 
        unlist(lapply(strsplit(lab_name, "_"), "[", 1)), 
        NA)
      
    }
    if(!is.null(clustering_method)){
      cdi_return_df["Cluster_method"] <- clustering_method
    }
    cdi_return <- apply(cand_lab_df, 
                        2, 
                        calculate_CDI_oneset,
                        sub_gcmat = sub_gcmat, 
                        batch_label = batch_label, 
                        cell_size_factor = cell_size_factor,
                        BPPARAM = BPPARAM, 
                        lrt_pval_threshold = lrt_pval_threshold) 
    cdi_return_df["CDI_AIC"] <- unlist(lapply(cdi_return, "[[", 1))
    cdi_return_df["CDI_BIC"] <- unlist(lapply(cdi_return, "[[", 2))
    cdi_return_df["Negllk"] <- unlist(lapply(cdi_return, "[[", 3))
    cdi_return_df["N_cluster"] <- apply(cand_lab_df, 2, function(x) length(unique(x)))
    return(cdi_return_df)
  }
}



# ------------------------------------------------------------------------------
#                          Lung
# ------------------------------------------------------------------------------
setwd("/hpc/group/xielab/jf243/lung/")
X_sc = readRDS("lung_size_factor.rds")
info = readRDS("lung_cell_info.rds")
X_sub = as.matrix(readMM("lung_X_sub.txt"))
data_name = "lung"

# recode since intCriteria only allows integer cluster
bc_lab = recode(info$celltype, 
                "AT1" = "1", "AT2" = "2", "Basal" = "3", "B Cells" = "4", 
                "cDCs" = "5", "Ciliated" = "6", "Differentiating Ciliated" = "7", 
                "Endothelial Cells" = "8", "Fibroblasts" = "9", "HAS1 High Fibroblasts" = "10", 
                "KRT5-/KRT17+" = "11", "Lymphatic Endothelial Cells" = "12", "Macrophages" = "13", 
                "Mast Cells" = "14", "Mesothelial Cells" = "15", "Monocytes" = "16", 
                "MUC5AC+ High" = "17", "MUC5B+" = "18", "Myofibroblasts" = "19", 
                "NK Cells" = "20", "pDCs" = "21", "Plasma Cells" = "22", 
                "PLIN2+ Fibroblasts" = "23", "Proliferating Epithelial Cells" = "24", 
                "Proliferating Macrophages" = "25", "Proliferating T Cells" = "26", "SCGB3A2+" = "27", 
                "SCGB3A2+ SCGB1A1+" = "28", "Smooth Muscle Cells" = "29", "T Cells" = "30", "Transitional AT2" = "31")

# paste0(paste0("\"", names(table(info$celltype)), "\" = \"", seq_len(length(names(table(info$celltype)))), "\", "), collapse = '')



mnames = c("CDI", "CH", "Connectivity", "DB", "Dunn", "Gamma", "SD_Scat","SD_Dis", "Silhouette", "XB")
nm = length(mnames)
# time_comp = data.frame(dataset = rep(data_name, nm), method = mnames, time = numeric(nm), ncell = rep(ncol(X_sub), nm))
time_comp = readRDS(paste0("/hpc/home/jf243/time_complexity_", data_name, ".rds"))


## CDI
# Start = Sys.time()
# aicbic_df = calculate_CDI(
#   sub_gcmat = X_sub,
#   Seurat_obj = NULL,
#   cand_lab_df = bc_lab, 
#   cell_size_factor = X_sc, 
#   batch_label = NULL,
#   lrt_pval_threshold = 0.01,
#   clustering_method = NULL,
#   ncore = 1)
# End = Sys.time()
# time_comp[time_comp$method == "CDI", "time"] = as.double(difftime(End, Start, units = "secs"), units = "secs")



tX = t(X_sub) *1.0 
rm(X_sub)

## CH
# Start = Sys.time()
# print(Start)
# intCriteria(tX, as.integer(bc_lab), "Calinski_Harabasz")
# End = Sys.time()
# time_comp[time_comp$method == "CH", "time"] = as.double(difftime(End, Start, units = "secs"), units = "secs")
# saveRDS(time_comp, paste0("/hpc/home/jf243/time_complexity_", data_name, ".rds"))


## DB
# Start = Sys.time()
# print(Start)
# intCriteria(tX, as.integer(bc_lab), "Davies_Bouldin")
# End = Sys.time()
# time_comp[time_comp$method == "DB", "time"] = as.double(difftime(End, Start, units = "secs"), units = "secs")
# saveRDS(time_comp, paste0("/hpc/home/jf243/time_complexity_", data_name, ".rds"))

## Dunn
# Start = Sys.time()
# print(Start)
# intCriteria(tX, as.integer(bc_lab), "Dunn")
# End = Sys.time()
# time_comp[time_comp$method == "Dunn", "time"] = as.double(difftime(End, Start, units = "secs"), units = "secs")
# saveRDS(time_comp, paste0("/hpc/home/jf243/time_complexity_", data_name, ".rds"))

## Gamma
# Start = Sys.time()
# print(Start)
# intCriteria(tX, as.integer(bc_lab), "Gamma")
# End = Sys.time()
# time_comp[time_comp$method == "Gamma", "time"] = as.double(difftime(End, Start, units = "secs"), units = "secs")
# saveRDS(time_comp, paste0("/hpc/home/jf243/time_complexity_", data_name, ".rds"))

## SD_Scat
Start = Sys.time()
print(Start)
intCriteria(tX, as.integer(bc_lab), "SD_Scat")
End = Sys.time()
time_comp[time_comp$method == "SD_Scat", "time"] = as.double(difftime(End, Start, units = "secs"), units = "secs")
saveRDS(time_comp, paste0("/hpc/home/jf243/time_complexity_", data_name, ".rds"))


## SD_Dis
Start = Sys.time()
print(Start)
intCriteria(tX, as.integer(bc_lab), "SD_Dis")
End = Sys.time()
time_comp[time_comp$method == "SD_Dis", "time"] = as.double(difftime(End, Start, units = "secs"), units = "secs")
saveRDS(time_comp, paste0("/hpc/home/jf243/time_complexity_", data_name, ".rds"))


## Silhouette
Start = Sys.time()
print(Start)
intCriteria(tX, as.integer(bc_lab), "Silhouette")
End = Sys.time()
time_comp[time_comp$method == "Silhouette", "time"] = as.double(difftime(End, Start, units = "secs"), units = "secs")
saveRDS(time_comp, paste0("/hpc/home/jf243/time_complexity_", data_name, ".rds"))



## XB
Start = Sys.time()
print(Start)
intCriteria(tX, as.integer(bc_lab), "Xie_Beni")
End = Sys.time()
time_comp[time_comp$method == "XB", "time"] = as.double(difftime(End, Start, units = "secs"), units = "secs")
saveRDS(time_comp, paste0("/hpc/home/jf243/time_complexity_", data_name, ".rds"))


## Connectivity
Start = Sys.time()
print(Start)
connectivity(Data = tX, clusters = bc_lab)
End = Sys.time()
time_comp[time_comp$method == "Connectivity", "time"] = as.double(difftime(End, Start, units = "secs"), units = "secs")

saveRDS(time_comp, paste0("/hpc/home/jf243/time_complexity_", data_name, ".rds"))
