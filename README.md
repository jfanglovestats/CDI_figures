# Clustering Deviation Index (CDI): A robust and accurate internal measure for evaluating scRNA-seq data clustering

Jiyuan Fang, Cliburn Chan, Kouros Owzar, Liuyang Wang, Diyuan Qin, Qi-Jing Li, and Jichun Xie

## Description


This repository contains codes to reproduce analysis. For details see the [biorxiv preprint](https://www.biorxiv.org/content/10.1101/2022.01.03.474840v1). Please create a GitHub issue if you have questions. 

Clustering deviation index (CDI) package is now available as an [R package](https://github.com/jichunxie/CDI). 


## Simulation

**SDX/sdX.R** ($X = 1, 2, 3, 4$) - Codes for simulating setting SD1 (10 equal-sized clusters), SD2 (rare cell populations), SD3 (main type and subtypes), and SD4 (splatter simulated dataset).

- Simulation setting

  Intermediate results: sdX_paraset.rds

- Generate candidate clustering labels

  Intermediate results: ./clustering_labels/

- Calculate CDI

  Intermediate results: sdX_cdi_df.rds, sdX_benchmark_return.rds

- Lineplot & Heatmap -- Fig.4A and Fig.4B

- UMAPS from selected Features -- Fig.3


`**SD2/sd2_rare_sensitivity.R**` Codes for sensitivity analysis of the number of cells in rare cell population -- Fig. S11


## CT26WT

**ct26wt_filter.R** - Filter genes with zero bulk-seq values. 

Output: WT1_nzbulk.rds

**ct26wt_zp_fitting.R** - Scatter plot in Fig.1A

Input: WT1_nzbulk.rds

**ct26wt_gof.R** - Chi-squared goodness-of-fit test in Fig.S1A

Input: WT1_nzbulk.rds

Intermediate results: wt1_gs_nb_pval_df.rds, wt1_gc_zinb_pval_df.rds, wt1_gc_nb_pval_df.rds, wt1_gc_zinb_pval_df.rds, wt1_discoveries.rds

**ct26wt_fpkm_gof.R**  - Chi-squared goodness-of-fit test in Fig.S1B. FPKM normalized counts. 

Input: WT1_nzbulk.rds

Intermediate results: wt1_fpkm_gs_nb_pval_df.rds, wt1_fpkm_gc_zinb_pval_df.rds, wt1_fpkm_gc_nb_pval_df.rds, wt1_fpkm_gc_zinb_pval_df.rds, wt1_fpkm_discoveries.rds

**ct26wt_tpm_gof.R**  - Chi-squared goodness-of-fit test in Fig.S1C. TPM normalized counts. 

Input: WT1_nzbulk.rds

Intermediate results: wt1_tpm_gs_nb_pval_df.rds, wt1_tpm_gc_zinb_pval_df.rds, wt1_tpm_gc_nb_pval_df.rds, wt1_tpm_gc_zinb_pval_df.rds, wt1_tpm_discoveries.rds


## TCELL (3k)

`**tcell_filter.R** - Select five distinct T-CELL groups and filter genes and cells`

Output: Tcell_5type_filtered.rds, Tcell_5type_filtered_labels.rds

**tcell_zp_fitting.R** - Scatter plot in Fig.1B

Input: Tcell_5type_filtered.rds, Tcell_5type_filtered_labels.rds


**tcell_gof.R** - Chi-squared goodness-of-fit test in Fig.S2

Input: Tcell_5type_filtered.rds, Tcell_5type_filtered_labels.rds

Intermediate results: tcell_ctc_nb_pval_vec.rds, tcell_ctc_zinb_pval_vec.rds, tcell_cts_nb_pval_vec.rds, tcell_cts_zinb_pval_vec.rds


**tcell.R** Codes for evaluating T-CELL clustering with CDI. 

Input: Tcell_5type_filtered.rds, Tcell_5type_filtered_labels.rds


- Generate candidate clustering labels

  Intermediate results: ./clustering_labels/

- Calculate CDI

  Intermediate results: tcell_cdi_df.rds, tcell_benchmark_return.rds

- Lineplot & Heatmap -- Fig.4A and Fig.4B

- UMAPS from selected Features -- Fig.3

- WDS- and VST-selected feature genes lead to different CDI-BIC scores -- Fig.S8

- CDI performance on T-CELL using different numbers of WDS-selected feature genes -- Fig.S9

- Heatmaps for six five-cluster candidate label sets -- Fig.S10

## CORTEX (7k)


**cortex_filter.R** - Filter genes and cells

Output: cortex_filtered.rds (large, not included in this folder), cortex_paraset.rds

**cortex_zp_fitting.R** - Scatter plot in Fig.S4


**cortex_gof.R** - Chi-squared goodness-of-fit test in Fig.S3

Input: cortex_filtered.rds, cortex_paraset.rds

Intermediate results: cortex_ctc_nb_pval_vec.rds, cortex_ctc_zinb_pval_vec.rds, cortex_cts_nb_pval_vec.rds, cortex_cts_zinb_pval_vec.rds


**cortex.R** Codes for evaluating CORTEX clustering with CDI. 

Input: cortex_filtered.rds, cortex_paraset.rds

- Generate candidate clustering labels

  Intermediate results: ./clustering_labels/

- Calculate CDI

  Intermediate results: tcell_cdi_df.rds, tcell_benchmark_return.rds

- Lineplot & Heatmap -- Fig.4A and Fig.4B

- UMAPS from selected Features -- Fig.3



## RETINA (27k)


**retina_filter.R** - Filter genes and cells

Output: retina_filtered.rds (large, not included in this folder), retina_paraset.rds, retina_filtered_info.rds

**retina_zp_fitting.R** - Scatter plot in Fig.S6


**retina_gof.R** - Chi-squared goodness-of-fit test in Fig.S5

Input: retina_filtered.rds, retina_paraset.rds, retina_filtered_info.rds

Intermediate results: retina_ctc_nb_pval_vec.rds, retina_ctc_zinb_pval_vec.rds, retina_cts_nb_pval_vec.rds, retina_cts_zinb_pval_vec.rds


**retina.R** Codes for evaluating RETINA clustering with CDI. 

Input: retina_filtered.rds, retina_paraset.rds, retina_filtered_info.rds

- Generate candidate clustering labels

  Intermediate results: ./clustering_labels/

- Calculate CDI

  Intermediate results: tcell_cdi_df.rds, tcell_benchmark_return.rds

- Lineplot & Heatmap -- Fig.4A and Fig.4B

- UMAPS from selected Features -- Fig.3


## LUNG (114k)

**lung.R** Codes for evaluating LUNG clustering with CDI. 

- Lineplot & Heatmap -- Fig.S11 and Fig.S12

- UMAPS for subpopulations -- Fig.S13

**lung_cluster.R** Codes to run in cluster for evaluating LUNG clustering with CDI. 

**lung_cluster_bash.q** Bash script to run lung_cluster.R 



## COVID (1,251k)


**covid.R** Codes for evaluating COVID clustering with CDI. 

- Generate candidate clustering labels

  Intermediate results: covid_lab_df.rds (large, not included in this folder)

- Lineplot -- Fig.S16

- Heatmap of reference labels -- Fig.S15

**covid_cluster_l2.R** Codes to run in cluster for evaluating COVID reference layer 2 label with CDI. 

**covid_cluster_l2_bash.q** Bash script to run covid_cluster_l2.R

**covid_cluster_s1.R** Codes to run in cluster for evaluating COVID clustering label 1 (Seurat with 12 clusters) with CDI. 

**covid_cluster_l2_bash.q** Bash script to run covid_cluster_s1.R

Intermediate results: covid_cdi_df.rds, covid_benchmark_return_l2.rds, covid_benchmark_return_l3.rds


## Metrics

Corresponds to Fig.5

**metrics.R** Codes for calculating all internal and external metrics.

**retina_cluster.R** Codes to calculate internal and external metrics for RETINA in cluster. 

Input: RETINA/retina_filtered_info.rds, cortex_X_sub.rds (large, not included in this folder)

**retina_cluster_bash.R** Bash script to run retina_cluster.R

**cortex_cluster.R** Codes to calculate internal and external metrics for CORTEX in cluster. 

Input: CORTEX/cortex_paraset.rds, cortex_X_sub.rds (large, not included in this folder)

**cortex_cluster_bash.R** Bash script to run cortex_cluster.R

Results: Metrics/Results

## Scalability

Corresponds to Fig.6

**time.R** Codes to record running time for internal indices of SD1-SD4, T-CELL, and CORTEX

**time_retina_cluster.R** Codes to record running time for internal indices of RETINA in cluster. 

**time_retina_cluster_bash.q** Bash script to run time_retina_cluster.R

**time_lung_cluster.R** Codes to record running time for internal indices of LUNG in cluster. 

**time_lung_cluster_bash.q** Bash script to run time_lung_cluster.R


**time_lung_cluster_all_labels.R** Codes to record running time for CDI of LUNG across all candidate labels in cluster. 

**time_lung_cluster_all_labels_bash.q** Bash script to run time_lung_cluster_all_labels.R

The running time for COVID data was recorded in COVID/covid_cluster_l2.R
## Utils

**clustering.R** - wrappers of clustering functions (CDI, K-Means, HC, Seurat, SC3, and spectral)


