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
library(ggsci)
library(RColorBrewer)
library(grid) 


ncore = 40

# -----------------------------------------------------
#          part 1.  functions                      
# -----------------------------------------------------
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
#          part 2.  data                   
# -----------------------------------------

# for random mislabel cells

set.indx = "retina"
info = readRDS(paste0(set.indx, "_filtered_info.rds"))
qc2_mat = readRDS(paste0(set.indx, "_filtered.rds"))
scmat = qc2_mat[, info$Batch == "Batch1"]
info_b1 = info[info$Batch == "Batch1",]

mt_name = c("Muller_glia", "OFF", "ON", "RBC")

cellsample = list()
ns = length(mt_name)
for(i in 1:ns){
  cellsample[[i]] <- assign(mt_name[i], scmat[, info_b1$main_type == mt_name[i]])
}
names(cellsample) = mt_name

set.seed(1)

sample_gene = sample(nrow(scmat), size = 2000)
qc2_list = lapply(cellsample, function(data){ data[sample_gene,]})
qc2_mat = do.call(cbind, qc2_list)


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
cur_mt = NULL
for(i in 1:ns){
  cur_mt = c(cur_mt, rep(mt_name[i], ng))
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

zp_ctc_nb_df = data.frame(obs_zp = obs_zp, cell_type = cur_mt, est_zp = est_zp_ctc_nb)
zp_ctc_zinb_df = data.frame(obs_zp = obs_zp, cell_type = cur_mt, est_zp = est_zp_ctc_zinb)
zp_ctsp_nb_df = data.frame(obs_zp = obs_zp, cell_type = cur_mt, est_zp = est_zp_ctsp_nb)
zp_ctsp_zinb_df = data.frame(obs_zp = obs_zp, cell_type = cur_mt, est_zp = est_zp_ctsp_zinb)


zp_ct_scatter = function(org_df){
  
  ## compress figures
  df = data.frame(obs_zp = round(org_df$obs_zp, 2),
                  cell_type = org_df$cell_type, 
                  est_zp = round(org_df$est_zp, 2))
  df = df[!duplicated(df),]
  
  ## scatter plot
  p1 = ggplot(data = df, mapping = aes(x = est_zp, y = obs_zp - est_zp, color = cell_type, shape = cell_type)) + 
    scale_color_npg() + geom_point(alpha = 1, size = 1) +
    ylim(-0.6, 0.6) + xlim(0, 1) +
    geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "black", size = 1.5) +
    theme_bw() + 
    theme(axis.text=element_text(size=15), 
          axis.title = element_text(size=15),
          axis.text = element_blank(),
          # legend.position = "none",
          axis.title = element_blank())
  p2 = p1 +  guides(colour = guide_legend(override.aes = list(shape = 1:length(unique(org_df$cell_type)))))
  return(p2)
  
}


pdf(paste0(set.indx, "_zp_fitting.pdf"),  width = 8, height = 5)
p1 = zp_ct_scatter(zp_ctc_nb_df)
p2 = zp_ct_scatter(zp_ctc_zinb_df)
p3 = zp_ct_scatter(zp_ctsp_nb_df)
p4 = zp_ct_scatter(zp_ctsp_zinb_df)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()
