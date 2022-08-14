#### Data source

## Data download
# The RData format dataset was given in 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81904 
# (GSE81904_bipolar_data_Cell2016.Rdata.gz)
#
## Preprocessing & benchmark labels
# Follow the preprocessing steps in 
# https://github.com/broadinstitute/BipolarCell2016/blob/master/BCanalysis.pdf,
# we got a filtered count matrix with 13,166 genes and 27,499 cells as shown in 
# page 2 of BCanalysis.pdf.
#
# These 27,499 cells were clustered by the authors to 26 clusters (page 6 of BCanalysis.pdf), 
# Among these 26 clusters, 18 clusters (26,830 cells) were identified with biomarkers to 18 cell types. 
# The cell type names and corresponding markers are available in the Table S1 of the original paper. 
# Our analysis starts from the count matrix with 26,830 cells from 18 cell types.
# 
## Main type labels
# Based on the paper, we also group the cell types to 5 main types
# OFF cone bipolar cells: BC1A, BC1B, BC2, BC3A, BC3B, BC4
# On cone bipolar cells: BC5A, BC5B, BC5C, BC5D, BC6, BC7, BC8_BC9
# Photo receptors: Rod photo receptors, Cone photo receptors
# Muller glia: Muller glia
# Amacrine: Amacrine
# Rod bipolar cells: Rod bipolar cells

set.seed(155501)


bipolar = readRDS("retina.rds")
# renamed from filtered_count.rds 
# dim(bipolar)
# 13166 27499

info = readRDS("retina_cell_info.rds")
# renamed from cell_info.rds
cell_info = info %>% 
  mutate(main_type = recode(cell_type, "RBC" = "Rod_BC", 
                            "Cone_photoreceptors" = "Photoreceptor",
                            "Rod_photoreceptors" = "Photoreceptor", 
                            "BC1A" = "OFF_BC", "BC1B" = "OFF_BC", 
                            "BC2" = "OFF_BC", "BC3A" = "OFF_BC",
                            "BC3B" = "OFF_BC", "BC4" = "OFF_BC",
                            "BC5A" = "ON_BC", "BC5B" = "ON_BC", 
                            "BC5C" = "ON_BC", "BC5D" = "ON_BC", 
                            "BC6" = "ON_BC", "BC7" = "ON_BC", 
                            "BC8_BC9" = "ON_BC"),
         sub_type = recode(cell_type, "RBC" = "Rod_BC"))
#renamed as retina_cell_info.rds
# dim(info)
# 27499 4
cell_label = as.character(info[,2])


## select verified cell types
verified_celltype = c("RBC","Muller_glia","BC1A","BC1B","BC2","BC3A","BC3B","BC4",
                      "BC5A","BC5B","BC5C","BC5D","BC6","BC7","BC8_BC9", "Amacrine", 
                      "Rod_photoreceptors","Cone_photoreceptors")

cell_indx = c(1:length(cell_label))[cell_label %in% verified_celltype]
info = info[cell_indx,]
select_bipolar = bipolar[,cell_indx]
select_celllabel = info[,2]




PreprocessSelectGenes = function(scdata, min_nzc = 100, min_nzc_prop = 0.02){
  nzc_of_each_gene = rowSums(scdata>0)
  nc = ncol(scdata)
  select_gene = c(1:nrow(scdata))[(nzc_of_each_gene > min_nzc) | (nzc_of_each_gene/nc > min_nzc_prop)]
  # return gene index
  return(select_gene)
}


filter_gene = PreprocessSelectGenes(select_bipolar, min_nzc = 100, min_nzc_prop = 0.02)

X = as.matrix(select_bipolar[filter_gene,])
saveRDS(para.set, paste0(set.indx, "_paraset.rds"))

para.set = list(K = length(unique(cell_info$cell_type)),
                ncell = ncol(X),
                ngene = nrow(X),
                nbatch = length(unique(cell_info$Batch)),
                maintype = cell_info$main_type,
                subtype = cell_info$sub_type,
                Batch = cell_info$Batch,
                nfeature = 500,
                ncore = 20)
set.indx = "retina"

saveRDS(select_celllabel, paste0(set.indx, "_filtered_info.rds"))
saveRDS(X, paste0(set.indx, "_filtered.rds"))
saveRDS(para.set, paste0(set.indx, "_paraset.rds"))

