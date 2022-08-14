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
library(scales)
library(ggsci)

# -----------------------------------------
#          part 1   functions                    
# -----------------------------------------


qc1 = function(scdata){
  nzg_of_each_cell = colSums(scdata>0)
  ng = nrow(scdata)
  mt_pnt = colSums(scdata[18773:ng,])/colSums(scdata)
  return(scdata[, (mt_pnt<0.1) & ((nzg_of_each_cell > 300) | (nzg_of_each_cell/ng > 0.03))])
}

qc2 = function(scdata){
  nzc_of_each_gene = rowSums(scdata>0)
  nc = ncol(scdata)
  select_gene = c(1:nrow(scdata))[(nzc_of_each_gene > 50) | (nzc_of_each_gene/nc > 0.01)]
  return(select_gene)
}




## FPKM
## https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-021-02936-w

FPKM = function(count, glen_vec){
  FPM = sweep(count* (10^6), 2, FUN = "/", colSums(count))
  return(sweep(FPM, 1, FUN = "/", (glen_vec / 10^3)))
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

#------------------- gene-common NB ----------------------------------#
neg_nb_size_logsum_fix_r = function(input_parameter, x_vec, sc_vec, r){
  mu = input_parameter
  s_mu = sc_vec * mu
  return(sum(-lgamma(x_vec + r) + lgamma(r) + lgamma(x_vec + 1)- x_vec * log(s_mu) 
             - r*log(r) + (x_vec+r) *log(r + s_mu) ))
}


nb_size_mle_fix_r = function(x_vec, sc_vec, iter_r){
  if(sum(x_vec) == 0){
    return(0.001)
  }
  avg = mean(x_vec); s2 = var(x_vec)
  init_mu = avg
  nb_est = nlminb(init_mu, 
                  objective = neg_nb_size_logsum_fix_r, 
                  gradient = NULL, 
                  lower = 1e-6, 
                  upper = 1e3,
                  x_vec = x_vec,
                  sc_vec = sc_vec,
                  r = iter_r)
  return(nb_est$par)
}


neg_nb_size_logsum_fix_mu = function(input_parameter, gcmat, sc_vec, mu_vec){
  r = input_parameter
  nllk_sum = 0
  for(g in 1:length(mu_vec)){
    s_mu = sc_vec * mu_vec[g]
    x_vec = gcmat[g,]
    nllk_sum = nllk_sum + sum(-lgamma(x_vec + r) + lgamma(r) + lgamma(x_vec + 1)- x_vec * log(s_mu) - r*log(r) + (x_vec+r) *log(r + s_mu) )
  }
  return(nllk_sum)
}


nb_size_mle_fix_mu = function(gcmat, sc_vec, iter_mu){
  nb_est = nlminb(5, 
                  objective = neg_nb_size_logsum_fix_mu, 
                  gradient = NULL, 
                  lower = 1e-6, 
                  upper = 1e3,
                  gcmat = gcmat,
                  sc_vec = sc_vec,
                  mu_vec = iter_mu)
  return(nb_est$par)
  
}

#------------------- gene-common ZINB ----------------------------------#

neg_zinb_size_logsum_fix_rpi = function(input_parameter, x_vec, sc_vec, r, pi){
  mu = input_parameter[1]
  s_mu = sc_vec * mu
  zero_prob = pi + (1-pi)*(r/(s_mu + r))^r
  neg_log_non_zero_prob = -lgamma(x_vec+r)+lgamma(r)+lgamma(x_vec+1)-x_vec*log(s_mu)-r*log(r)+(x_vec+r)*log(r+s_mu)
  return(sum(ifelse(x_vec == 0, -log(zero_prob), -log(1-pi) + neg_log_non_zero_prob)))
  
}


zinb_size_mle_fix_rpi = function(x_vec, sc_vec, iter_r, iter_pi){
  if(sum(x_vec) == 0){
    return(c(1e-6, 1-1e-6))
  }
  avg = mean(x_vec); s2 = var(x_vec)
  init_mu = avg
  zinb_est = nlminb(init_mu, 
                    objective = neg_zinb_size_logsum_fix_rpi, 
                    gradient = NULL, 
                    lower = 1e-6, 
                    upper = 1e3,
                    x_vec = x_vec,
                    sc_vec = sc_vec,
                    r = iter_r,
                    pi = iter_pi)
  return(zinb_est$par)
}


neg_zinb_size_logsum_fix_mu = function(input_parameter, gcmat, sc_vec, mu_vec){
  r = input_parameter[1]
  pi = input_parameter[2]
  nllk_sum = 0
  for(g in 1:length(mu_vec)){
    s_mu = sc_vec * mu_vec[g]
    x_vec = gcmat[g,]
    zero_prob = pi + (1-pi)*(r/(s_mu + r))^r
    neg_log_non_zero_prob = -lgamma(x_vec+r)+lgamma(r)+lgamma(x_vec+1)-x_vec*log(s_mu)-r*log(r)+(x_vec+r)*log(r+s_mu)
    nllk_sum = nllk_sum + sum(ifelse(x_vec == 0, -log(zero_prob), -log(1-pi) + neg_log_non_zero_prob))
  }
  return(nllk_sum)
}


zinb_size_mle_fix_mu = function(gcmat, sc_vec, iter_mu){
  zinb_est = nlminb(c(5, 0.5), 
                    objective = neg_zinb_size_logsum_fix_mu, 
                    gradient = NULL, 
                    lower = c(1e-6,0),  
                    upper = c(1e6, 1),
                    gcmat = gcmat,
                    sc_vec = sc_vec,
                    mu_vec = iter_mu)
  return(zinb_est$par)
  
}


obs_point_mass = function(input_vec, max_val){
  obs_count = numeric(max_val + 1)
  for(i in 1:max_val){
    obs_count[i] = sum(input_vec == i-1)
  }
  obs_count[(max_val + 1)] = max(0, sum(table(input_vec)) - sum(obs_count))
  return(obs_count)
}


est_point_mass_nb_size = function(sc_vec, input_para, max_val){
  mu = input_para[1]
  r = input_para[2]
  s_mu = sc_vec * mu
  exp_count = numeric(max_val + 1)
  for(i in 1:max_val){
    exp_count[i] = sum(dnbinom(i-1, mu = s_mu, size = r))
  }
  exp_count[(max_val + 1)] = max(0, length(sc_vec) - sum(exp_count))
  return(exp_count)
}


est_point_mass_zinb_size = function(sc_vec, input_para, max_val){
  mu = input_para[1]
  r = input_para[2]
  pi = input_para[3]
  s_mu = sc_vec * mu
  exp_count = numeric(max_val + 1)
  exp_count[1] = sum(c(pi + (1-pi)*dnbinom(0, mu = s_mu, size = r)))
  for(i in 2:max_val){
    exp_count[i] = (1-pi)*sum(dnbinom(i-1, mu = s_mu, size = r))
  }
  exp_count[(max_val + 1)] = max(0, length(sc_vec) - sum(exp_count))
  return(exp_count)
}


gc_nb_size_mle = function(train_data, train_sc, niter = 20){
  data = as.matrix(train_data)
  ng = nrow(data)
  r_est = 1
  mu_est = numeric(ng)
  for(i in 1:niter){
    cat("i=", i, "\n")
    last_mu = mu_est; last_r = r_est
    mu_est = foreach(g = seq_len(ng), .combine = 'c') %do% {
      nb_size_mle_fix_r(data[g,], train_sc, iter_r = r_est)
    }
    r_est = nb_size_mle_fix_mu(gcmat = data, sc_vec = train_sc, iter_mu = mu_est)
    if((mean((last_mu - mu_est)^2) < 1e-6) & (last_r -r_est)^2 <1e-4) break
  }
  return(list(mu = mu_est, r = r_est))
}


gc_zinb_size_mle = function(train_data, train_sc, niter = 20){
  data = as.matrix(train_data)
  ng = nrow(data)
  r_est = 1; pi_est = 0.5
  mu_est = numeric(ng)
  for(i in 1:niter){
    cat("i=", i, "\n")
    last_mu = mu_est; last_r = r_est; last_pi = pi_est
    mu_est = foreach(g = seq_len(ng), .combine = 'c') %dopar% {
      zinb_size_mle_fix_rpi(data[g,], sc_vec = train_sc, iter_r = r_est, iter_pi = pi_est)
    }
    zinb_return = zinb_size_mle_fix_mu(gcmat = data, sc_vec = train_sc, iter_mu = mu_est)
    r_est = zinb_return[1]; pi_est = zinb_return[2]
    if((mean((last_mu - mu_est)^2) < 1e-6) & (mean((last_r - r_est)^2) < 1e-4) &(mean((last_pi - pi_est)^2) < 1e-4)) break
  }
  return(list(mu = mu_est, r = r_est, pi = pi_est))
}



# -----------------------------------------------------------
#          part 2   GOF test                    
# -----------------------------------------------------------



input_data = readRDS("WT1_nzbulk.rds")
pre_cell = qc1(input_data)
process_data = as.matrix(pre_cell[qc2(pre_cell),])


## Obtain gene length from EDASeq (with ensembl ID)
library (EDASeq)
ensembl_list <- rownames(input_data)
gene_length <- getGeneLengthAndGCContent(ensembl_list, "mmu")

intersect_gene = intersect(rownames(gene_length)[!is.na(gene_length[, "length"])], rownames(process_data))

gof_data = round(FPKM(count = process_data[intersect_gene,], 
                      glen_vec = gene_length[intersect_gene, "length"]))


set.seed(1)
train_indx = sample(ncol(gof_data), round(ncol(gof_data)/2))
train_set = as.matrix(gof_data[,train_indx])
test_set = as.matrix(gof_data[, -train_indx])

sc_vec = rep(1, ncol(gof_data))
ng = nrow(gof_data)
train_sc_vec = sc_vec[train_indx]
test_sc_vec = sc_vec[-train_indx]
bin_max_val = 3


cl <- makeCluster(20, type="FORK")
registerDoParallel(cl)





# ------------------ gene-specific NB -------------------------------------- #
gs_nb_pval_vec <- foreach(g = seq_len(ng), .combine = 'c') %dopar% {
  test_gvec = test_set[g,]
  train_gvec =  train_set[g,]
  obs_freq_vec = obs_point_mass(test_gvec, bin_max_val)
  input_parameter = nb_size_mle(train_gvec, train_sc_vec)
  est_freq_vec = est_point_mass_nb_size(test_sc_vec, input_parameter, bin_max_val)
  dof = (bin_max_val + 1) - 1
  # dof = 1
  test_stat = sum((obs_freq_vec- est_freq_vec)^2 / est_freq_vec)
  pval = pchisq(test_stat, df = dof, lower.tail = F)
  return(pval)
}
gs_nb_pval_df = data.frame(pval = gs_nb_pval_vec)
saveRDS(gs_nb_pval_df, "wt1_fpkm_gs_nb_pval_df.rds")

# ------------------ gene-specific ZINB -------------------------------------- #
gs_zinb_pval_vec <- foreach(g = seq_len(ng), .combine = 'c') %dopar% {
  test_gvec = test_set[g,]
  train_gvec =  train_set[g,]
  obs_freq_vec = obs_point_mass(test_gvec, bin_max_val)
  input_parameter = zinb_size_mle(train_gvec, train_sc_vec)
  est_freq_vec = est_point_mass_zinb_size(train_sc_vec, input_parameter, bin_max_val)
  dof = (bin_max_val + 1) - 2
  test_stat = sum((obs_freq_vec- est_freq_vec)^2 / est_freq_vec)
  pval = pchisq(test_stat, df = dof, lower.tail = F)
  return(pval)
}
gs_zinb_pval_df = data.frame(pval = gs_zinb_pval_vec)
saveRDS(gs_zinb_pval_df, "wt1_fpkm_gs_zinb_pval_df.rds")


# ------------------ gene-common NB -------------------------------------- #
gc_nb_para = gc_nb_size_mle(train_set, train_sc_vec)
gc_nb_pval_vec <- foreach(g = seq_len(ng), .combine = 'c') %dopar% {
  test_gvec = test_set[g,]
  train_gvec =  train_set[g,]
  obs_freq_vec = obs_point_mass(test_gvec, bin_max_val)
  input_parameter = c(gc_nb_para$mu[g], gc_nb_para$r)
  est_freq_vec = est_point_mass_nb_size(train_sc_vec, input_parameter, bin_max_val)
  dof = (bin_max_val + 1) - 1
  test_stat = sum((obs_freq_vec- est_freq_vec)^2 / est_freq_vec)
  pval = pchisq(test_stat, df = dof, lower.tail = F)
  return(pval)
}
gc_nb_pval_df = data.frame(pval = gc_nb_pval_vec)
saveRDS(gc_nb_pval_df, "wt1_fpkm_gc_nb_pval_df.rds")


# ------------------ gene-common ZINB -------------------------------------- #
gc_zinb_para = gc_zinb_size_mle(train_set, train_sc_vec)
gc_zinb_pval_vec <- foreach(g = seq_len(ng), .combine = 'c') %dopar% {
  test_gvec = test_set[g,]
  train_gvec =  train_set[g,]
  obs_freq_vec = obs_point_mass(test_gvec, bin_max_val)
  input_parameter = c(gc_zinb_para$mu[g], gc_zinb_para$r, gc_zinb_para$pi)
  est_freq_vec = est_point_mass_zinb_size(train_sc_vec, input_parameter, bin_max_val)
  dof = (bin_max_val + 1) - 1
  test_stat = sum((obs_freq_vec- est_freq_vec)^2 / est_freq_vec)
  pval = pchisq(test_stat, df = dof, lower.tail = F)
  return(pval)
}
gc_zinb_pval_df = data.frame(pval = gc_zinb_pval_vec)
saveRDS(gc_zinb_pval_df, "wt1_fpkm_gc_zinb_pval_df.rds")

stopCluster(cl)


# -----------------------------------------------------------
#          part 3   Figure                    
# -----------------------------------------------------------

gc_nb_pval_df = readRDS("wt1_fpkm_gc_nb_pval_df.rds")
gc_zinb_pval_df = readRDS("wt1_fpkm_gc_zinb_pval_df.rds")
gs_nb_pval_df = readRDS("wt1_fpkm_gs_nb_pval_df.rds")
gs_zinb_pval_df = readRDS("wt1_fpkm_gs_zinb_pval_df.rds")


pval_range = c(0.005, 0.01, 0.05, 0.1)
model_list = list(gc_nb_pval_df[,1], gc_zinb_pval_df[,1], gs_nb_pval_df[,1], gs_zinb_pval_df[,1])
pval_df = data.frame(
  Pval = paste0("<", rep(pval_range, rep(4, 4))),
  Model = rep(c("Gene-common NB", "Gene-common ZINB", "Gene-specific NB", "Gene-specific ZINB"), 4),
  Percentage = NA)
for(i in 1:4){
  for(j in 1:4){
    pval_df$Percentage[(i-1)*4 + j] = round(mean(model_list[[j]] < pval_range[i]),4)
  }
}

pdf("wt1_fpkm_pval_dist.pdf")
ggplot(data = pval_df, aes(group = as.factor(Model))) +
  labs(x = "P-value") + theme_bw() +
  geom_bar(aes(x= as.factor(Pval), y = Percentage,
               fill = as.factor(Model)),  stat = "identity", position=position_dodge()) +
  scale_fill_npg(name = "Model")
dev.off()
discoveries = data.frame(
  model = c("gcnb", "gczinb", "gsnb", "gszinb"),
  pval_rej = c(sum(gc_nb_pval_df[,1] < 0.05), sum(gc_zinb_pval_df[,1] < 0.05),
               sum(gs_nb_pval_df[,1] < 0.05), sum(gs_zinb_pval_df[,1] < 0.05)),
  padj_rej = c(sum(p.adjust(gc_nb_pval_df[,1], method = "BH") <0.05),
               sum(p.adjust(gc_zinb_pval_df[,1], method = "BH") <0.05),
               sum(p.adjust(gs_nb_pval_df[,1], method = "BH") <0.05),
               sum(p.adjust(gs_zinb_pval_df[,1], method = "BH") <0.05))
)
saveRDS(discoveries, "wt1_fpkm_discoveries.rds")





