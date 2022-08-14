library(Matrix)
library(matrixStats)
library(MASS)
library(doParallel) 
library(Rcpp)
library(truncnorm)
library(data.table)
library(mclust)
library(Seurat)
library(gridExtra)
library(ggsci)
library(RColorBrewer)
library(grid) 


ncore = 40

# ------------------------------------------------
#            part 1.  functions                      
# ------------------------------------------------

## size factor for each cell
size_factor = function(scmat){
  scmat[scmat == 0] = 0.5
  nc = ncol(scmat)
  log_scmat = log(scmat)
  ref_size = exp(rowMeans(log_scmat))
  ratio_to_ref = sweep(scmat, 1, ref_size, "/")
  cell_size_factor = colMedians(ratio_to_ref)
  return(cell_size_factor)
}

#------------------- gene-specific NB ----------------------------------#
## negative sum of log-likelihood for NB distribution with size factor 
neg_nb_size_logsum = function(input_parameter, x_vec, sc_vec){
  mu = input_parameter[1]
  r = input_parameter[2]
  s_mu = sc_vec * mu
  return(sum(-lgamma(x_vec + r) + lgamma(r) + lgamma(x_vec + 1)- x_vec * log(s_mu) 
             - r*log(r) + (x_vec+r) *log(r + s_mu) ))
}
## Note: the prob in rnbinom() is not the same as the commonly used def. 
## prob = 1- p

## MLE of NB distribution parameter with the existence of size factor
nb_size_mle = function(x_vec, sc_vec){
  if(sum(x_vec) == 0){
    return(c(0.001, 0.2))
  }
  avg = mean(x_vec); s2 = var(x_vec)
  init_mu = avg
  tmp_init_r = avg^2/(s2 - avg)
  init_r_0 = ifelse(is.na(tmp_init_r)|(tmp_init_r < 0.1), 0.1, tmp_init_r)
  init_r = ifelse(init_r_0 > 50, 50, init_r_0)
  nb_est = nlminb(c(init_mu, init_r), 
                  objective = neg_nb_size_logsum, 
                  gradient = NULL, 
                  lower = c(1e-6, 1e-6), 
                  upper = c(1e3, 1e6),
                  x_vec = x_vec,
                  sc_vec = sc_vec)
  return(nb_est$par)
}

#------------------- gene-specific ZINB ----------------------------------#
neg_zinb_size_logsum = function(input_parameter, x_vec, sc_vec){
  mu = input_parameter[1]
  r = input_parameter[2]
  pi = input_parameter[3]
  s_mu = sc_vec * mu
  zero_prob = pi + (1-pi)*(r/(s_mu + r))^r
  neg_log_non_zero_prob = -lgamma(x_vec+r)+lgamma(r)+lgamma(x_vec+1)-x_vec*log(s_mu)-r*log(r)+(x_vec+r)*log(r+s_mu)
  return(sum(ifelse(x_vec == 0, -log(zero_prob), -log(1-pi) + neg_log_non_zero_prob)))
}


## MLE of NB distribution parameter with the existence of size factor
zinb_size_mle = function(x_vec, sc_vec){
  if(sum(x_vec) == 0){
    return(c(0.001, 0.2, 0.9))
  }
  zp = sum(x_vec ==0)/length(x_vec)
  avg = mean(x_vec); s2 = var(x_vec)
  init_pi = zp/2
  init_mu = avg
  tmp_init_r = avg^2/(s2 - avg)
  init_r_0 = ifelse(is.na(tmp_init_r)|(tmp_init_r < 0.1), 0.1, tmp_init_r)
  init_r = ifelse(init_r_0 > 50, 50, init_r_0)
  nb_est = nlminb(c(init_mu, init_r, init_pi), 
                  objective = neg_zinb_size_logsum, 
                  gradient = NULL, 
                  lower = c(1e-6, 1e-6, 0), 
                  upper = c(1e3, 1e6, zp),
                  x_vec = x_vec,
                  sc_vec = sc_vec)
  return(nb_est$par)
}


## gene specific NB
gs_nb_size_zp = function(data){
  sc_vec = size_factor(data)
  ng = nrow(data)
  gsp_zcp <- foreach(g = seq_len(ng), .combine = 'c') %dopar% {
    gvec = data[g,]
    mle2 = nb_size_mle(gvec, sc_vec)
    median((mle2[2]/(sc_vec * mle2[1] + mle2[2]))^mle2[2])
    #mean((mle2[2]/(sc_vec * mle2[1] + mle2[2]))^mle2[2])
  }
  return(gsp_zcp)
}



## gene specific ZINB
gs_zinb_size_zp = function(data){
  sc_vec = size_factor(data)
  #sc_vec = rep(1, ncol(data))
  ng = nrow(data)
  gsp_zcp <- foreach(g = seq_len(ng), .combine = 'c') %dopar% {
    gvec = data[g,]
    mle3 = zinb_size_mle(gvec, sc_vec)
    median((1-mle3[3])*(mle3[2]/(sc_vec * mle3[1] + mle3[2]))^mle3[2] + mle3[3])
  }
  return(gsp_zcp)
}


# -----------------------------------------
#          part 3.  data                   
# -----------------------------------------

set.seed(100) # for random selection of genes

hrvatin_info = read.csv("hrvatin_celltypes.csv")
hrvatin = readRDS("hrvatin.rds")
# dim(hrvatin) = 19155, 10000
hrvatin_cell_name = data.frame(X = colnames(hrvatin))
merged_info = merge(hrvatin_cell_name, hrvatin_info, sort = F)
labeled_maintype_indx = c(1:nrow(merged_info))[!is.na(merged_info$maintype)]
labeled_maintype_info = merged_info[labeled_maintype_indx,]
rownames(labeled_maintype_info) = as.character(labeled_maintype_info$X)
scmat = hrvatin[, rownames(labeled_maintype_info)]
info = labeled_maintype_info



celllabel = info[, "maintype"]
celltype = unique(celllabel)
ns = length(celltype)
cell_indx = c(1:ns)
cellsample = list()
for(i in 1:ns){
  cellsample[[i]] <- assign(celltype[i], scmat[, celllabel == celltype[i]])
}



qc1ht = function(scdata){
  nzg_of_each_cell = colSums(scdata>0)
  ng = nrow(scdata)
  return(scdata[, (nzg_of_each_cell > 300) | (nzg_of_each_cell/ng > 0.03)])
}


qc1_cellsample = lapply(cellsample, qc1ht)
qc1_mat = do.call(cbind, qc1_cellsample)


qc2ht = function(scdata){
  nzc_of_each_gene = rowSums(scdata>0)
  nc = ncol(scdata)
  select_gene = c(1:nrow(scdata))[(nzc_of_each_gene > 50) | (nzc_of_each_gene/nc > 0.01)]
  return(select_gene)
}

filter_gene = qc2ht(qc1_mat)
sample_gene = sample(filter_gene, size = 2000)


qc2_list = lapply(qc1_cellsample, function(data){ data[sample_gene,]})
qc2_mat = do.call(cbind, qc2_list)

## qc2_mat is the same as cortex_filtered.rds if random selection of genes is not applied

train_indx = function(mat) {
  nc = ncol(mat)
  return(sample(1:nc, round(nc/2)))
}



train_indx_list = lapply(qc2_list, train_indx)
train_list = list(); test_list = list()

for(i in 1:ns){
  train_list[[i]] = qc2_list[[i]][,train_indx_list[[i]]]
  test_list[[i]] = qc2_list[[i]][,-train_indx_list[[i]]]
}



ng = nrow(test_list[[1]])
cell_type = NULL
for(i in 1:ns){
  cell_type = c(cell_type, rep(celltype[i], ng))
}

train_mat = do.call(cbind, train_list)


obs_zp = unlist(lapply(test_list, function(mat) rowSums(mat == 0) / ncol(mat)))



cl <- makeCluster(25, type="FORK")
registerDoParallel(cl)
est_zp_ctc_nb = rep(gs_nb_size_zp(train_mat), ns)
est_zp_ctc_zinb = rep(gs_zinb_size_zp(train_mat), ns)
est_zp_ctsp_nb = unlist(lapply(train_list, function(mat) gs_nb_size_zp(mat)))
est_zp_ctsp_zinb = unlist(lapply(train_list, function(mat) gs_zinb_size_zp(mat)))
stopCluster(cl)

# saveRDS(est_zp_ctc_nb, "est_zp_ctc_nb.rds")
# saveRDS(est_zp_ctc_zinb, "est_zp_ctc_zinb.rds")
# saveRDS(est_zp_ctsp_nb, "est_zp_ctsp_nb.rds")
# saveRDS(est_zp_ctsp_zinb, "est_zp_ctsp_zinb.rds")


zp_ctc_nb_df = data.frame(obs_zp = obs_zp, cell_type = cell_type, est_zp = est_zp_ctc_nb)
zp_ctc_zinb_df = data.frame(obs_zp = obs_zp, cell_type = cell_type, est_zp = est_zp_ctc_zinb)
zp_ctsp_nb_df = data.frame(obs_zp = obs_zp, cell_type = cell_type, est_zp = est_zp_ctsp_nb)
zp_ctsp_zinb_df = data.frame(obs_zp = obs_zp, cell_type = cell_type, est_zp = est_zp_ctsp_zinb)




zp_ct_scatter = function(org_df){
  
  ## compress figures
  df = data.frame(obs_zp = round(org_df$obs_zp, 2),
                  cell_type = org_df$cell_type, 
                  est_zp = round(org_df$est_zp, 2))
  df = df[!duplicated(df),]
  
  ## scatter plot
  p1 = ggplot(data = df, mapping = aes(x = est_zp, y = obs_zp - est_zp, shape = cell_type)) + 
    geom_point(aes(colour = cell_type), alpha = 0.5, size = 1) + ylim(-0.7, 0.7) + xlim(0, 1) +
    geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "black", size = 1.5) +
    theme_bw() + scale_shape_manual(values=1:8, name = "Cell-types") +
    scale_color_npg() +
    theme(axis.text=element_text(size=15), 
          axis.title = element_text(size=15),
          axis.text = element_blank(),
          # legend.position = "none",
          axis.title = element_blank())
  p2 = p1 +  guides(colour = guide_legend(override.aes = list(shape = 1:length(unique(org_df$cell_type)))))
  return(p1)
  
}

pdf(paste0(set.indx, "_zp_fitting.pdf"),  width = 8, height = 5)
p1 = zp_ct_scatter(zp_ctc_nb_df)
p2 = zp_ct_scatter(zp_ctc_zinb_df)
p3 = zp_ct_scatter(zp_ctsp_nb_df)
p4 = zp_ct_scatter(zp_ctsp_zinb_df)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()


