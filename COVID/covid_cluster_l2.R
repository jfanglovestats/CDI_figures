library(Seurat)
library(Matrix)
library(dplyr)
library(matrixStats)
library(truncnorm)
library(data.table)
library(mclust)
library(ggsci)
library(cluster)
library(BiocParallel)



# ------------------------------------------------------------------------------
#                            functions
# ------------------------------------------------------------------------------

one_batch_feature_gene_rank <- function(
  gcmat = NULL,
  Seurat_obj = NULL,
  method = "wds", 
  nfeature = 500,
  zp_threshold = 0.95){
  ## Gene selection via wds
  if(method == "wds"){
    if(is.null(gcmat)){
      gcmat <- as.matrix(Seurat_obj@assays$RNA@counts)
    }
    ng <- nrow(gcmat)
    nc <- ncol(gcmat)
    zp <- rowMeans(gcmat == 0)
    zp_indx <- c(seq_len(ng))[zp < zp_threshold]
    mu_g <- rowMeans(gcmat); var_g <- rowVars(gcmat)
    phi_est <- (var_g-mu_g)/mu_g^2
    phi_df <- data.frame(indx = seq_len(ng), phi_est = phi_est)
    phi_df[zp_indx, "phi_rank"] <- rank(-phi_est[zp_indx])
    ## Assign genes with NA rank to a large value to remove NA
    phi_df$phi_rank[is.na(phi_df$phi_rank)] <- ng + 1000
    return(phi_df$phi_rank)
  }
  ## Gene selection via vst
  if(method == "vst"){
    ## if no Seurat object input: create a Seurat object
    if(is.null(Seurat_obj)){
      ng <- nrow(gcmat)
      nc <- ncol(gcmat)
      if(is.null(colnames(gcmat))){
        colnames(gcmat) <- paste0("c", seq_len(nc))
      }
      if(is.null(rownames(gcmat))){
        rownames(gcmat) <- paste0("g", seq_len(ng))
      }
      gene_indx_name <- data.frame(indx = c(seq_len(ng)), grank = rep(ng, ng))
      rownames(gene_indx_name) <- rownames(gcmat)
      Seurat_obj <- CreateSeuratObject(counts = as.data.frame(gcmat))
    }
    Seurat_obj <- NormalizeData(Seurat_obj, verbose = FALSE)
    Seurat_obj <- FindVariableFeatures(Seurat_obj, selection.method = "vst", nfeatures = nfeature, verbose = FALSE)
    tmp_df <- data.frame(gname = Seurat_obj@assays[["RNA"]]@var.features, grank = seq_len(nfeature))
    gene_indx_name[as.character(tmp_df$gname), "grank"] <- tmp_df$grank
    return(gene_indx_name$grank)
  }
}
feature_gene_selection <- function(
  gcmat = NULL, 
  Seurat_obj = NULL,
  method = "wds",
  nfeature = 500,
  batch_label = NULL,
  zp_threshold = 0.95){
  if((method != "wds") & (method != "vst")){
    stop("Method not available!")
  }
  if((is.null(gcmat)) & (is.null(Seurat_obj))){
    stop("One of gcmat and Seurat_obj needs to be non-empty!")
  } else{
    ng <- ifelse(is.null(Seurat_obj), nrow(gcmat), nrow(Seurat_obj@assays$RNA@counts))
    if(ng < nfeature){
      stop("The number of feature genes to select exceeds the number of genes.")
    }
  }
  if(is.null(batch_label) | (length(unique(batch_label)) == 1)){
    feature_bool = (one_batch_feature_gene_rank(gcmat = gcmat, 
                                                Seurat_obj = Seurat_obj, 
                                                method = method, 
                                                nfeature = nfeature, 
                                                zp_threshold = zp_threshold) <= nfeature)
    return(c(seq_len(ng))[feature_bool])
  } else{
    if(is.null(gcmat)){
      gcmat <- as.matrix(Seurat_obj@assays$RNA@counts)
    }
    ## split gcmat columns according to the batch labels
    mat_list <- lapply(split(seq_along(batch_label), batch_label), 
                       function(indx, mat) mat[,indx], 
                       mat = gcmat)
    ## rank genes in each matrix
    gene_rank_list <- lapply(mat_list, 
                             one_batch_feature_gene_rank, 
                             method = method, 
                             nfeature = nfeature, 
                             zp_threshold = zp_threshold)
    ## return top nfeature genes
    gene_rank_mat <- do.call(cbind, gene_rank_list)
    gene_min_rank <- rowMins(gene_rank_mat)
    names(gene_min_rank) <- seq_len(nrow(gcmat))
    as.numeric(names(gene_min_rank)[order(gene_min_rank)[seq_len(nfeature)]])
    return(sort(as.numeric(names(gene_min_rank)[order(gene_min_rank)[seq_len(nfeature)]])))
  }
  
}


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
  nReject = sum(lrt_test_pval < lrt_pval_threshold)
  # cat("llk = ", total_nllk, "nrej = ", nReject, "\n")
  return(list(NegLLK = total_nllk, nReject = nReject))
  
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


size_factor <- function(gcmat){
  gcmat[gcmat == 0] <- 0.5
  nc <- ncol(gcmat)
  log_gcmat <- log(gcmat)
  ref_size <- exp(rowMeans(log_gcmat))
  ratio_to_ref <- sweep(gcmat, 1, ref_size, "/")
  cell_size_factor <- colMedians(ratio_to_ref)
  return(cell_size_factor)
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
#                            Data & label 
# ------------------------------------------------------------------------------

# setwd("/hpc/group/xielab/jf243/covid_immune/")
set.indx = "covid"
info = readRDS("covid_cellinfo.rds")
X_sub = as.matrix(readMM("covid_X_sub.txt"))
X_sc = rep(1, ncol(X_sub))


Start = Sys.time()
cdi_return = calculate_CDI(sub_gcmat = X_sub,
                          Seurat_obj = NULL,
                          cand_lab_df = info$predicted.celltype.l2,
                          cell_size_factor = X_sc,
                          batch_label = info$dataset,
                          lrt_pval_threshold = 0.01,
                          clustering_method = NULL,
                          ncore = 10)

End = Sys.time()
difftime(End, Start, units = "secs")
print(cdi_return)
saveRDS(cdi_return, "cdi_covid_benchmark_l2.rds")


