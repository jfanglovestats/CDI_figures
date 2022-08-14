library(sfsmisc) # for plot axis 
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

set.seed(1100155)
parallel_ncore = 40

# This code requires the output (WT1_nzbulk.rds) from ct26wt_filter.R


# -------------------------------------------------------------
#          part 1   shared functions                      
# -------------------------------------------------------------

## sample zinb data
sample_zinb = function(n, mu, r, pi){
  ber_vec = rbinom(n, size = 1, prob = 1-pi)
  nb_vec = rnbinom(n, size = r, mu = mu)
  return(ifelse(ber_vec == 0, rep(0,n), nb_vec))
}


qc1 = function(scdata){
  nzg_of_each_cell = colSums(scdata>0)
  ng = nrow(scdata)
  # remove genes with mitochondrial proportion greater than 10%
  mt_pnt = colSums(scdata[18773:ng,])/colSums(scdata)
  return(scdata[, (mt_pnt<0.1) & ((nzg_of_each_cell > 300) | (nzg_of_each_cell/ng > 0.03))])
}

qc2 = function(scdata){
  nzc_of_each_gene = rowSums(scdata>0)
  nc = ncol(scdata)
  select_gene = c(1:nrow(scdata))[(nzc_of_each_gene > 50) | (nzc_of_each_gene/nc > 0.01)]
  return(select_gene)
}


zp_scatter = function(org_df){
  ## compress figures
  df = data.frame(obs_zp = round(org_df$obs_zp, 3),
                  est_zp = round(org_df$est_zp, 3))
  df = df[!duplicated(df),]
  
  ## scatter plot
  p1 = ggplot(data = df, mapping = aes(x = est_zp, y = obs_zp - est_zp)) + 
    geom_point(colour = "#00A087FF", alpha = 0.5, size = 2) + ylim(-0.1, 0.25) + xlim(0, 1) +
    geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "black", size = 1.5) +
    theme_bw() + 
    labs(x = "Estimated zero proportion", y = "Zero proportion difference (observed - estimated)") +
    theme(axis.text=element_text(size=15), 
          axis.title = element_blank(),
          legend.position = "none") 
  return(p1)
  
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




# ---------------------------------------------------------------------
#          part 2.   four main functions                      
# ---------------------------------------------------------------------

## gene specific NB
gs_nb_size_zp = function(data){
  sc_vec = size_factor(data)
  ng = nrow(data)
  gsp_zcp <- foreach(g = seq_len(ng), .combine = 'c') %dopar% {
    gvec = data[g,]
    mle2 = nb_size_mle(gvec, sc_vec)
    median((mle2[2]/(sc_vec * mle2[1] + mle2[2]))^mle2[2])
  }
  return(gsp_zcp)
}



## gene specific ZINB
gs_zinb_size_zp = function(data){
  sc_vec = size_factor(data)
  ng = nrow(data)
  gsp_zcp <- foreach(g = seq_len(ng), .combine = 'c') %dopar% {
    gvec = data[g,]
    mle3 = zinb_size_mle(gvec, sc_vec)
    median((1-mle3[3])*(mle3[2]/(sc_vec * mle3[1] + mle3[2]))^mle3[2] + mle3[3])
  }
  return(gsp_zcp)
}



## gene common NB
gc_nb_size_zp = function(data, niter = 20){
  data = as.matrix(data)
  r_est = 1
  ng = nrow(data)
  gc_zcp = numeric(ng)
  sc_vec = size_factor(data)
  for(i in 1:niter){
    cat("i=", i, "\n")
    mu_est = foreach(g = seq_len(ng), .combine = 'c') %dopar% {
      nb_size_mle_fix_r(data[g,], sc_vec, iter_r = r_est)
    }
    r_est = nb_size_mle_fix_mu(gcmat = data, sc_vec = sc_vec, iter_mu = mu_est)
    
    last_zcp = gc_zcp
    gc_zcp <- foreach(g = seq_len(ng), .combine = 'c') %dopar% {
      gvec = data[g,]
      median((r_est/(sc_vec * mu_est[g] + r_est))^r_est)
    }
    if((mean((last_zcp - gc_zcp)^2) < 1e-6)) break
  }
  return(gc_zcp)
}


## gene common ZINB
gc_zinb_size_zp = function(data, niter = 20){
  ng = nrow(data)
  gc_zcp = numeric(ng)
  sc_vec = size_factor(data)
  # sc_vec = rep(1, ncol(data))
  r_est = 1; pi_est = 0.5
  for(i in 1:niter){
    cat("i=", i, "\n")
    mu_est = foreach(g = seq_len(ng), .combine = 'c') %dopar% {
      zinb_size_mle_fix_rpi(data[g,], sc_vec, iter_r = r_est, iter_pi = pi_est)
    }
    zinb_return = zinb_size_mle_fix_mu(gcmat = data, sc_vec = sc_vec, iter_mu = mu_est)
    r_est = zinb_return[1]; pi_est = zinb_return[2]
    
    last_zcp = gc_zcp
    gc_zcp <- foreach(g = seq_len(ng), .combine = 'rbind') %dopar% {
      gvec = data[g,]
      median((1-pi_est)*(r_est/(sc_vec * mu_est[g] + r_est))^r_est + pi_est)
    }
    if((mean((last_zcp - gc_zcp)^2) < 1e-6)) break
  }
  return(gc_zcp)
}









# ------------------------------------------------------------------
#                   part 3.   data result                     
# ------------------------------------------------------------------


input_data = readRDS("WT1_nzbulk.rds")
pre_cell = qc1(input_data)
pre_data = as.matrix(pre_cell[qc2(pre_cell),])
train_indx = sample(ncol(pre_data), round(ncol(pre_data)/2))
train_set = as.matrix(pre_data[,train_indx])
test_set = as.matrix(pre_data[, train_indx])
obs_zp = rowSums(test_set == 0)/ncol(test_set)

cl <- makeCluster(parallel_ncore, type="FORK")
registerDoParallel(cl)
df_gs_nb = data.frame(obs_zp = obs_zp, est_zp = gs_nb_size_zp(train_set))
df_gs_zinb = data.frame(obs_zp = obs_zp, est_zp = gs_zinb_size_zp(train_set))
df_gc_nb = data.frame(obs_zp = obs_zp, est_zp = gc_nb_size_zp(train_set))
df_gc_zinb = data.frame(obs_zp = obs_zp, est_zp = gc_zinb_size_zp(train_set))
stopCluster(cl)

# saveRDS(df_gs_nb, "df_gs_nb.rds")
# saveRDS(df_gs_zinb, "df_gs_zinb.rds")
# saveRDS(df_gc_nb, "df_gc_nb.rds")
# saveRDS(df_gc_zinb, "df_gc_zinb.rds")
# 
# 
# df_gs_nb = readRDS("df_gs_nb.rds")
# df_gs_zinb = readRDS("df_gs_zinb.rds")
# df_gc_nb = readRDS("df_gc_nb.rds")
# df_gc_zinb = readRDS("df_gc_zinb.rds")


pdf("wt1_zp_fitting.pdf",  width = 8, height = 5)
p1 = zp_scatter(df_gc_nb)
p2 = zp_scatter(df_gc_zinb)
p3 = zp_scatter(df_gs_nb)
p4 = zp_scatter(df_gs_zinb)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

