#!/usr/bin/env Rscript
##
#.libPaths("/project/uvm_mckay/shamima.akter/stereo_liver/R/RStudio-x86_64-pc-linux-gnu-library/4.4") ## remote server
setwd("/project/uvm_mckay/shamima.akter/stereo_liver/data")
getwd()
library(SingleR)
library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
library(reticulate)
library(dplyr)
library(scRNAseq)
adata<- readH5AD("B03022E1_bin50out.h5ad")
dim(adata)
##getting gene expression matrix
exprs<- assay(adata) 
rownames(exprs) 
colnames(exprs)
rownames(exprs)
dim(exprs)
if (class(exprs)[1] == "dgTMatrix") {
  exprs <- as(exprs, "dgCMatrix")
}
test_se<- SummarizedExperiment(list(logcounts = exprs)) ## preprocessing and lognormalized 
dim(test_se)
#reference liver single cell data
ref <- readRDS("Sc_liverRef.RDS")##sc reference data in .rds format
dim(ref)
head(ref)
ref1<-ref$CellType
ref
ref_expr <- GetAssayData(ref) ## extract expression data
##Here preprocessing was done for dgCMatrix of reference sc gene expression data in a similar coding for expression  data of stereo-seq data earlier.
rownames(ref_expr) 
colnames(ref_expr)

ref_se <- SummarizedExperiment(list(logcounts = ref_expr))## preprocessing and lognormalized 
dim(ref_se) 
##prediction of cell types
pred <- SingleR(test = test_se, ref = ref_se,labels =ref$CellType)
pred

colData(adata)$cell_type <- pred$labels # prepare a table of predictions
# Save predictions
write.csv(data.frame(Cell=rownames(pred), Label=pred$labels), 'cell_annotations_onlyLiver4815.csv', row.names = FALSE)
adata$cell_type
df <- as.data.frame(colData(adata))
#data as df can be savesd as .RDS   or .csv file
# Write to .h5ad
writeH5AD(adata, "4815_Refannotation.h5ad")## can be implemented for other analysis
head(df[, c("x", "y", "cell_type")])
library(ggplot2)
library(ggimage)
ggplot(df, aes(x = x, y = y, color = cell_type)) +
  geom_point(size = 4.0, alpha = 0.8) +
  coord_fixed() +  # 👈 Coordinate setting goes here
  theme_void(base_size = 14) +  # 👈 Then apply the theme
  labs(color = "Cell Type") +
  ggtitle("Spatial Distribution of Cell Types") 
  #This plot can be generated with custom color palate with cell types
ggsave("spatial_plot14815livercells_SingleR.png", width = 12, height = 10, dpi = 600, bg = "transparent")

