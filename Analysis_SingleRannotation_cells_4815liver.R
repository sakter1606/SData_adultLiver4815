##
#.libPaths("/project/uvm_mckay/shamima.akter/stereo_liver/R/RStudio-x86_64-pc-linux-gnu-library/4.4")
setwd("/project/uvm_mckay/shamima.akter/stereo_liver/data")
getwd()
library(SingleR)
library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
library(reticulate)
library(IRIS)
library(dplyr)
library(scRNAseq)
adata<- readH5AD("B03022E1_bin50out.h5ad")
dim(adata)
#adata<-adata[1:5000,1:5000]
#  getting gene expression matrix
exprs<- assay(adata) 

#coords <- colData(adata)[, c("x", "y")]  
rownames(exprs) 
colnames(exprs)
rownames(exprs)
dim(exprs)
if (class(exprs)[1] == "dgTMatrix") {
  exprs <- as(exprs, "dgCMatrix")
}
#coords <- colData(Spatial)[,c("x","y")]
#dim(coords)
test_se<- SummarizedExperiment(list(logcounts = exprs))
dim(test_se)
#reference liver single cell data
#sc_input_liver_meta <- readRDS("Liver_anno_meta.RDS")##reference scRNA-Seq data for liver cell.metadata
#sc_input_liver_meta 
#sc_input_liver_count <- readRDS("Liver_count.RDS")##reference scRNA-Seq data for liver cells count
#sc_input_liver_count 
ref <- readRDS("Sc_liverRef.RDS")##sc reference data in .rds format
dim(ref)
head(ref)
#ref<-ref[1:21125,1:25000]
#ref1<-ref$CellType
ref
ref_expr <- GetAssayData(ref)
if (class(ref_expr)[1] == "dgTMatrix") {
  ref_expr <- as(ref_expr, "dgCMatrix")
}

rownames(ref_expr) 
colnames(ref_expr)

ref_se <- SummarizedExperiment(list(logcounts = ref_expr))
dim(ref_se) 
##prediction of cell types
pred <- SingleR(test = test_se, ref = ref_se,labels =ref$CellType)
pred
pred$labels  # this should be a character vector
colData(adata)$cell_type <- pred$labels
# Save predictions
write.csv(data.frame(Cell=rownames(pred), Label=pred$labels), 'cell_annotations_onlyLiver4815.csv', row.names = FALSE)
adata$cell_type
df <- as.data.frame(colData(adata))
saveRDS(df,file =  "refannotation_Liver4815.RDS")
df
write.csv(df, file = "liver_annotation_4815.csv")
#df2 <- as.SingleCellExperiment(ref)
head(adata)
# Write to .h5ad
writeH5AD(adata, "4815_Refannotation.h5ad")## can be implemented for other analysis
#df1=read.csv("liver_annotation_4815.csv")
#head(df1[, c("x", "y", "cell_type")])
head(df[, c("x", "y", "cell_type")])
library(ggplot2)
library(ggimage)
# Define custom colors for each cell type
my_colors <- c(
  "B cells" = "#2f3635", "T cells" = "#ff7f0e", "NK cells" = "#287d0f","Plasma cells"="#ed4e88","NK T cells"="#eda84e",
  "Dendritic cells" ="#b09f3e" , "Myofibroblasts" = "#bf454f",
  "Kupffer cells" = "#594a6b","Cholangiocytes" = "#7d12c9",
  "Hepatocytes" = "#79b7c7","Unknown cells" = "#b6e043",
  "Blood vascular endothelial cells"= "#6e2704", "Liver sinusoidal endothelial cells"="#a38577","Neutrophils"="#3043bf","Proliferative dendritic cells"= "#e0ba7b"
  # Add al cell types here with distinct colors
)
ggplot(df, aes(x = x, y = y, color = cell_type)) +
  #geom_point(size = 0.8, alpha = 1.0) +
  geom_point(size = 4.0, alpha = 0.8) +
  coord_fixed() +  # 👈 Coordinate setting goes here
  theme_void(base_size = 14) +  # 👈 Then apply the theme
  scale_color_manual(values = my_colors) +
  labs(color = "Cell Type") +
  ggtitle("Spatial Distribution of Cell Types") +
  
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
#ggsave("spatial_plot14815liver_Hepaqtocytes.png", width = 12, height = 10, dpi = 600, bg = "transparent")

ggsave("spatial_plot14815livercells_SingleR.png", width = 12, height = 10, dpi = 600, bg = "transparent")

