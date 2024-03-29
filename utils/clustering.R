## SC3 
# BiocManager::install("SC3")
#devtools::install_github("davismcc/scater", build_vignettes = FALSE)
library(SingleCellExperiment)
library(SC3)
library(scater)

## Seurat
library(Seurat)

## Spectral
library(kernlab)

## CIDR
# devtools::install_github("VCCRI/CIDR")
library(cidr)


## PCA
gcmat_pc = function(gcmat, npc){
  dat_pr = prcomp(t(gcmat), scale = TRUE)[[5]][,1:npc]
  return(dat_pr)
}



## Seurat with default setting
seurat_lab = function(gcmat, res_vec){
  Seurat.obj = CreateSeuratObject(counts = gcmat, project = "sc20a", min.cells = 0, min.features = 0)
  Seurat.obj = NormalizeData(Seurat.obj, verbose = F)
  Seurat.obj = FindVariableFeatures(Seurat.obj, selection.method = "vst", nfeatures = 0.2*nrow(gcmat), verbose = F)
  Seurat.obj = ScaleData(Seurat.obj, verbose = F)
  Seurat.obj = RunPCA(Seurat.obj, npcs = 30, verbose = F)
  Seurat.obj = FindNeighbors(Seurat.obj, reduction = "pca", dims = 1:20)
  label_result = matrix(NA, nrow = ncol(gcmat), ncol = length(res_vec))
  num_cl = numeric(length(res_vec))
  for(i in 1:length(res_vec)){
    Seurat.obj = FindClusters(Seurat.obj, resolution = res_vec[i], verbose = F)
    df = data.frame(Seurat = Seurat.obj@meta.data$seurat_clusters)
    num_cl[i] = length(unique(df$Seurat))
    saveRDS(df, paste0(set.indx, "_Seurat_k", num_cl[i], "_labs.rds"))
    Seurat.obj@meta.data$seurat_clusters = NULL
  }
  return(num_cl)
}


# seurat_lab2 = function(seu_obj, res_vec){
#   label_result = matrix(NA, nrow = nrow(seu_obj@meta.data), ncol = length(res_vec))
#   num_cl = numeric(length(res_vec))
#   for(i in 1:length(res_vec)){
#     seu_obj = FindClusters(seu_obj, resolution = res_vec[i], verbose = F)
#     df = data.frame(Seurat = seu_obj@meta.data$seurat_clusters)
#     num_cl[i] = length(unique(df$Seurat))
#     #saveRDS(df, paste0(set.indx, "_Seurat_k", num_cl[i], "_labs.rds"))
#     saveRDS(df, paste0(set.indx, "_Seurat_k", num_cl[i], "_labs.rds"))
#     seu_obj@meta.data$seurat_clusters = NULL
#   }
#   return(num_cl)
# }




## Seurat from batch effect corrected object
seurat_lab_bc = function(seu_obj, res_vec, cell_name){
  DefaultAssay(seu_obj) <- "integrated"
  label_result = matrix(NA, nrow = ncol(seu_obj[["RNA"]]), ncol = length(res_vec))
  num_cl = numeric(length(res_vec))
  for(i in 1:length(res_vec)){
    seu_obj = FindClusters(seu_obj, resolution = res_vec[i], verbose = FALSE)
    df = data.frame(Seurat = seu_obj@meta.data$seurat_clusters)
    num_cl[i] = length(unique(df$Seurat))
    rownames(df) = rownames(Seurat.obj.integrated@meta.data)
    df_sorted = data.frame(Seurat = df[cell_name, ])
    rownames(df_sorted) = cell_name
    saveRDS(df_sorted, paste0(set.indx, "_Seurat_k", num_cl[i], "_labs.rds"))
    seu_obj@meta.data$seurat_clusters = NULL
    cat("finish res = ", res_vec[i],"\n")
  }
  return(num_cl)
}


## Spectral clustering on PCs
spec_lab_pc = function(gcmat_pc, K_vec){
  for(i in 1:length(K_vec)){
    sc = specc(gcmat_pc, centers = K_vec[i])
    df = data.frame(Spectral = sc@.Data)
    saveRDS(df, paste0(set.indx, "_Spectral_k", K_vec[i], "_labs.rds"))
    cat("K=", K_vec[i],"\n")
  }
}

## Hierarchical clustering on PCs
hc_lab_pc = function(gcmat_pc, K_vec){
  dist.obj = dist(gcmat_pc)
  hc.obj = hclust(dist.obj, method = "complete")
  for(i in 1:length(K_vec)){
    ## linkage (how to define new distance when merging)
    # complete: maximum
    # single: minimum
    df = data.frame(HC = cutree(hc.obj, k = K_vec[i]))
    saveRDS(df, paste0(set.indx, "_HC_k", K_vec[i], "_labs.rds"))
    cat("K=", K_vec[i],"\n")
  }
}

## K-Means on PCs
kmeans_lab_pc = function(gcmat_pc, K_vec, n_start = 3){
  for(i in 1:length(K_vec)){
    km.obj = kmeans(gcmat_pc, centers = K_vec[i], nstart = n_start)
    df = data.frame(KMeans = km.obj$cluster)
    saveRDS(df, paste0(set.indx, "_KMeans_k", K_vec[i], "_labs.rds"))
    cat("K=", K_vec[i],"\n")
  }
}

## CIDR on top 200 PCs
cidr_lab = function(gcmat, K_vec, numPCs = 200){
  sData = scDataConstructor(gcmat, tagType = "raw")
  cat("finish constructor \n")
  sData = determineDropoutCandidates(sData)
  cat("finish dropout \n")
  sData = wThreshold(sData)
  cat("finish wthreshold \n")
  sData = scDissim(sData)
  cat("finish scdissim \n")
  sData = scPCA(sData)
  cat("finish pca \n")
  sData = nPC(sData)
  for(i in 1:length(K_vec)){
    sData = scCluster(sData, nCluster = K_vec[i], nPC = numPCs)
    df = data.frame(CIDR = sData@clusters)
    saveRDS(df, paste0(set.indx, "_CIDR_k", K_vec[i], "_labs.rds"))
    cat("K=", K_vec[i],"\n")
  }
}

## SC3 with default settings
sc3_lab = function(gcmat, K_vec, n_cores = 40){
  ## build a sc-object for SC3
  sc.obj = SingleCellExperiment(assays = list(counts = as.matrix(gcmat),
                                              logcounts = log2(as.matrix(gcmat) + 1)))
  rowData(sc.obj)$feature_symbol = rownames(gcmat)
  sc.obj = runPCA(sc.obj)
  for(i in 1:length(K_vec)){
    sc3.obj = sc3(sc.obj, ks = K_vec[i], n_cores = n_cores, biology = FALSE, kmeans_nstart = 3)
    df = data.frame(SC3 = as.vector(colData(sc3.obj)[[1]]))
    saveRDS(df, paste0(set.indx, "_SC3_k", K_vec[i], "_labs.rds"))
    cat("K=", K_vec[i],"\n")
  }
}




sc3_lab_large = function(gcmat, K_vec){
  ## build a sc-object for SC3
  sc.obj = SingleCellExperiment(assays = list(counts = as.matrix(gcmat),logcounts = log2(as.matrix(gcmat) + 1)))
  rowData(sc.obj)$feature_symbol = rownames(gcmat)
  sc.obj = runPCA(sc.obj)
  sc3.obj = sc3(sc.obj, ks = K_vec, n_cores = ncl_vec, biology = FALSE, kmeans_nstart = 3, svm_num_cells = 5000)
  sc3.obj = sc3_run_svm(sc3.obj, ks = K_vec)
  for(i in 1:length(K_vec)){
    df = data.frame(SC3 = as.vector(colData(sc3.obj)[[i]]))
    #cat(table(df[,1]), "\n")
    saveRDS(df, paste0(set.indx, "_SC3_k", K_vec[i], "_labs.rds"))
  }
}





BS_resolution = function(ideal_K, gcmat, res_seq = seq(0.01, 5, by = 0.01)){
  Seurat.obj = CreateSeuratObject(counts = gcmat, project = "sc20a", 
                                  min.cells = 0, min.features = 0)
  Seurat.obj = NormalizeData(Seurat.obj, verbose = FALSE)
  Seurat.obj = FindVariableFeatures(Seurat.obj, selection.method = "vst", 
                                    nfeatures = 0.2*nrow(gcmat), verbose = FALSE)
  Seurat.obj = ScaleData(Seurat.obj, verbose = FALSE)
  Seurat.obj = RunPCA(Seurat.obj, npcs = 30, verbose = FALSE)
  Seurat.obj = FindNeighbors(Seurat.obj, reduction = "pca", 
                             dims = 1:20, verbose = FALSE)
  nres = length(res_seq)
  res_indx = round(nres / 2)
  max_iter = 20
  iter = 1
  while((iter <= max_iter) & (res_indx >= 1) & (res_indx <= nres)){
    cur_res = res_seq[res_indx]
    Seurat.obj = FindClusters(Seurat.obj, resolution = cur_res, verbose = F)
    cur_ncl = length(unique(Seurat.obj@meta.data$seurat_clusters))
    if(cur_ncl == ideal_K){
      return(cur_res)
    } else if(cur_ncl > ideal_K){
      cat("ncl = ", cur_ncl, "res=", cur_res,"\n")
      res_indx = round(res_indx / 2)
    } else if(cur_ncl < ideal_K){
      cat("ncl = ", cur_ncl, "res=", cur_res,"\n")
      res_indx = round(res_indx * 1.5)
    }
    iter = iter + 1
  }
  if(iter == max_iter){
    cat("last resolution = ", cur_res, ". Exceed maximum iteration!")
  }
  return("Not found!")
}