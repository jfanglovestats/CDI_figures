## 
library(Matrix)
bulk_seq = read.table("CT26_bulk_30k.txt", header=T, row.names=1)
#dim 30805     5
sc_ginfo = read.table("scRNAseq/features.tsv", sep="\t")
#dim 31053     3
intersect_gid = intersect(rownames(bulk_seq), sc_ginfo[,1])
# > length(intersect_gid)
# [1] 30805 
## this intersection contains all genes in bulk_seq
nz_bulk_gid = rownames(bulk_seq)[(bulk_seq[,1]!=0)|(bulk_seq[,2]!=0)|(bulk_seq[,3]!=0)|(bulk_seq[,4]!=0)]
# 18785 genes
WT1 =  readMM("scRNAseq/matrix.mtx")
rownames(WT1) = sc_ginfo[,1]
# saveRDS(WT1[nz_bulk_gid, ], "WT1_nzbulk.rds")

