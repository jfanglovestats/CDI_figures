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


# ----------------------------------------------------
#          part 1.  functions for estimation                   
# ----------------------------------------------------

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


bin_nb_nloglik = function(para, obs_freq, sc_vec, max_val){
  exp_freq = est_point_mass_nb_size(sc_vec, para, max_val)
  exp_prob = exp_freq / length(sc_vec)
  return(-sum(obs_freq * log(exp_prob)))
}


nb_est_para = function(obs_freq_tab, sc_vec, max_val){
  ## for some celltypes, it is possible that the
  ## number of non-zero counts is 0
  ## it will be hard to optimize (also might give unstable result)
  ## approximate with NB(p=0.999, r = 0.2) mu will be less than 0.01
  if(sum(obs_freq_tab>0) == 0){
    return(c(0.999, 0.2))
  }
  else{
    obs_count = c(0:max_val)
    avg = sum(obs_count * obs_freq_tab)/sum(obs_freq_tab)
    nb_bin_est = nlminb(c(avg, 1), 
                        objective = bin_nb_nloglik, 
                        gradient = NULL, 
                        lower = c(1e-3, 1e-5), 
                        upper = c(1e3, 1e5), 
                        obs_freq = obs_freq_tab,
                        sc_vec = sc_vec,
                        max_val = max_val)
    return(nb_bin_est$par)
  }
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


bin_zinb_nloglik = function(para, obs_freq, sc_vec, max_val){
  exp_freq = est_point_mass_zinb_size(sc_vec, para, max_val)
  exp_prob = exp_freq / length(sc_vec)
  return(-sum(obs_freq * log(exp_prob)))
}
zinb_est_para(c(123,0,0,0), rep(1,123), 3)

zinb_est_para = function(obs_freq_tab, sc_vec, max_val){
  ## for some celltypes, it is possible that the
  ## number of non-zero counts is 0
  ## it will be hard to optimize (also might give unstable result)
  ## approximate with NB(p=0.999, r = 0.2) mu will be less than 0.01
  if(sum(obs_freq_tab[-1]>0) == 0){
    return(c(0.0001, 0.2, 0.5))
  }
  else{
    obs_count = c(0:max_val)
    avg = sum(obs_count * obs_freq_tab)/sum(obs_freq_tab)
    nb_bin_est = nlminb(c(avg, 1, obs_freq_tab[1]/2), 
                        objective = bin_zinb_nloglik, 
                        gradient = NULL, 
                        lower = c(1e-3, 1e-5,0), 
                        upper = c(1e3, 1e5, obs_freq_tab[1]/sum(obs_freq_tab)), 
                        obs_freq = obs_freq_tab,
                        sc_vec = sc_vec,
                        max_val = max_val)
    return(nb_bin_est$par)
  }
}

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






# -----------------------------------------------------------
#                   part 2   data                  
# -----------------------------------------------------------

set.indx = "retina"
info = readRDS(paste0(set.indx, "_filtered_info.rds"))
qc2_mat = readRDS(paste0(set.indx, "_filtered.rds"))
scmat = qc2_mat[, info$Batch == "Batch1"]
info_b1 = info[info$Batch == "Batch1",]

mt_name = c("Amacrine", "Muller_glia", "OFF", "ON", "PR", "RBC")
ns = length(mt_name)

bin_max_val = 4
nsplit = 30
npar_core = 30
ngene_per_core = 400


for(i in 1:nsplit){
  saveRDS(scmat[((i-1)*ngene_per_core+1):(i*ngene_per_core),], paste0(set.indx, "qc2_mat_split",i,".rds"))
}

sc_vec = size_factor(qc2_mat)
labs = info_b1$main_type
ng = nrow(qc2_mat)


# ------------------ cell-type-common NB -------------------------------------- #
cl <- makeCluster(npar_core, type="FORK")
registerDoParallel(cl)

start = Sys.time()
foreach(i = seq_len(npar_core)) %dopar% {
  qc2_mat = readRDS(paste0(set.indx, "qc2_mat_split",i,".rds"))
  pval = numeric(ngene_per_core)
  for(g in 1:ngene_per_core){
    gvec = qc2_mat[g,]
    gvec_list = split(gvec, f = factor(labs))
    sc_list =  split(sc_vec, f = factor(labs))
    cur_obs_freq = obs_point_mass(gvec, bin_max_val)
    para_est = nb_est_para(cur_obs_freq, sc_vec, bin_max_val)
    test_stat = 0
    dof = ((bin_max_val + 1) - 1) * ns - 2
    for(k in 1:ns){
      ct_obs_freq = obs_point_mass(gvec_list[[k]], bin_max_val)
      ct_est_freq = est_point_mass_nb_size(sc_list[[k]], para_est, max_val = bin_max_val)
      test_stat = test_stat + sum((ct_obs_freq - ct_est_freq)^2 / ct_est_freq)
    }
    pval[g] = pchisq(test_stat, df = dof, lower.tail = F)
  }
  saveRDS(pval, paste0(set.indx, "_ctc_nb_pval_",i,".rds"))
}

stopCluster(cl)
end = Sys.time()
end-start

## 2 genes (/core) time: 29 secs
## 400 gene (/core) estimate: 1.61 hours -- real 40mins

pval_vec = NULL
for(j in 1:nsplit){
  pval_vec = c(pval_vec, readRDS(paste0(set.indx, "_ctc_nb_pval_",j,".rds")))
}
saveRDS(pval_vec, paste0(set.indx, "_ctc_nb_pval_df.rds"))


# ------------------ cell-type-common ZINB -------------------------------------- #
cl <- makeCluster(npar_core, type="FORK")
registerDoParallel(cl)

start = Sys.time()
foreach(i = seq_len(npar_core)) %dopar% {
  qc2_mat = readRDS(paste0(set.indx, "qc2_mat_split",i,".rds"))
  pval = numeric(ngene_per_core)
  for(g in 1:ngene_per_core){
    gvec = qc2_mat[g,]
    
    gvec_list = split(gvec, f = factor(labs))
    sc_list =  split(sc_vec, f = factor(labs))
    
    cur_obs_freq = obs_point_mass(gvec, bin_max_val)
    para_est = zinb_est_para(cur_obs_freq, sc_vec, bin_max_val)
    
    test_stat = 0
    dof = ((bin_max_val + 1) - 1) * ns - 3
    for(k in 1:ns){
      ct_obs_freq = obs_point_mass(gvec_list[[k]], bin_max_val)
      ct_est_freq = est_point_mass_zinb_size(sc_list[[k]], para_est, max_val = bin_max_val)
      test_stat = test_stat + sum((ct_obs_freq - ct_est_freq)^2 / ct_est_freq)
    }
    pval[g] = pchisq(test_stat, df = dof, lower.tail = F)
  }
  saveRDS(pval, paste0(set.indx, "_ctc_zinb_pval_",i,".rds"))
}

stopCluster(cl)

end = Sys.time()
end-start
#1.48 hour

pval_vec = NULL
for(j in 1:nsplit){
  pval_vec = c(pval_vec, readRDS(paste0(set.indx, "_ctc_zinb_pval_",j,".rds")))
}

saveRDS(pval_vec, paste0(set.indx, "_ctc_zinb_pval_df.rds"))





# ------------------ cell-type-specific NB -------------------------------------- #
cl <- makeCluster(npar_core, type = "FORK")
registerDoParallel(cl)

start = Sys.time()
foreach(i = seq_len(npar_core)) %dopar% {
  qc2_mat = readRDS(paste0(set.indx, "qc2_mat_split",i,".rds"))
  nb_pval = numeric(ngene_per_core)
  zinb_pval = numeric(ngene_per_core)
  
  for(g in 1:ngene_per_core){
    gvec = qc2_mat[g,]
    gvec_list = split(gvec, f = factor(labs))
    sc_list =  split(sc_vec, f = factor(labs))
    
    nb_test_stat = 0
    zinb_test_stat = 0
    nb_dof = ((bin_max_val + 1) - 2 -1) * ns
    zinb_dof = ((bin_max_val + 1) - 3 -1) * ns
    for(k in 1:ns){
      ct_obs_freq = obs_point_mass(gvec_list[[k]], bin_max_val)
      
      nb_para_est = nb_est_para(ct_obs_freq, sc_vec, bin_max_val)
      nb_ct_est_freq = est_point_mass_nb_size(sc_list[[k]], nb_para_est, max_val = bin_max_val)
      nb_test_stat = nb_test_stat + sum((ct_obs_freq - nb_ct_est_freq)^2 / nb_ct_est_freq)
      
      zinb_para_est = zinb_est_para(ct_obs_freq, sc_vec, bin_max_val)
      zinb_ct_est_freq = est_point_mass_zinb_size(sc_list[[k]], zinb_para_est, max_val = bin_max_val)
      zinb_test_stat = zinb_test_stat + sum((ct_obs_freq - zinb_ct_est_freq)^2 / zinb_ct_est_freq)
      
    }
    nb_pval[g] = pchisq(nb_test_stat, df = nb_dof, lower.tail = F)
    zinb_pval[g] = pchisq(zinb_test_stat, df = zinb_dof, lower.tail = F)
  }
  saveRDS(nb_pval, paste0(set.indx, "_cts_nb_pval_",i,".rds"))
  saveRDS(zinb_pval, paste0(set.indx, "_cts_zinb_pval_",i,".rds"))
}


stopCluster(cl)
end = Sys.time()
end-start
## 2 genes (/core) time: 3.47min
## 400 gene (/core) estimate: 11.56hours -- real <9h

pval_vec = NULL
for(j in 1:nsplit){
  pval_vec = c(pval_vec, readRDS(paste0(set.indx, "_cts_nb_pval_",j,".rds")))
}
saveRDS(pval_vec, paste0(set.indx, "_cts_nb_pval_df.rds"))



pval_vec = NULL
for(j in 1:nsplit){
  pval_vec = c(pval_vec, readRDS(paste0(set.indx, "_cts_zinb_pval_",j,".rds")))
}
saveRDS(pval_vec, paste0(set.indx, "_cts_zinb_pval_df.rds"))


ctc_nb_pval_df = data.frame(pval = ctc_nb_pval_vec)
ctc_zinb_pval_df = data.frame(pval = ctc_zinb_pval_vec)
cts_nb_pval_df = data.frame(pval = cts_nb_pval_vec)
cts_zinb_pval_df = data.frame(pval = cts_zinb_pval_vec)


# ctc_nb_pval = readRDS(paste0(set.indx, "_ctc_nb_pval_df.rds"))
# ctc_zinb_pval = readRDS(paste0(set.indx, "_ctc_zinb_pval_df.rds"))
# cts_nb_pval = readRDS(paste0(set.indx, "_cts_nb_pval_df.rds"))
# cts_zinb_pval = readRDS(paste0(set.indx, "_cts_zinb_pval_df.rds"))


pval_range = c(0.005, 0.01, 0.05, 0.1)
model_list = list(ctc_nb_pval, ctc_zinb_pval, cts_nb_pval, cts_zinb_pval)
pval_df = data.frame(
  Pval = paste0("<", rep(pval_range, rep(4, 4))),
  Model = rep(c("Cell-type-common NB", "Cell-type-common ZINB", 
                "Cell-type-specific NB", "Cell-type-specific ZINB"), 4),
  Percentage = NA)
for(i in 1:4){
  for(j in 1:4){
    pval_df$Percentage[(i-1)*4 + j] = round(mean(model_list[[j]] < pval_range[i]),4)
  }
}

pdf(paste0(set.indx, "_pval_dist.pdf"))
ggplot(data = pval_df, aes(group = as.factor(Model))) + 
  labs(x = "P-value") + theme_bw() + 
  geom_bar(aes(x= as.factor(Pval), y = Percentage, 
               fill = as.factor(Model)),  stat = "identity", position=position_dodge()) +
  scale_fill_npg(name = "Model") 
dev.off()




