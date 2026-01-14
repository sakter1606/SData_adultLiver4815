#!/usr/bin/env Rscript

setwd("/project/uvm_mckay/shamima.akter/stereo_liver/data")## set working directory where the data are
getwd()
library(zellkonverter)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(harmony)
library(HDF5Array)
library(SeuratDisk)
library(patchwork)
# 1. Read .h5ad
sce <- readH5AD("B03022E1_bin50out.h5ad")
# 2. Ensure counts are accessible
if (!"logcounts" %in% assayNames(sce)) {
  logcounts(sce) <- log1p(assay(sce, "X"))
}

# 3. Convert to Seurat object
seurat_obj <- as.Seurat(sce, counts = "X", data = "logcounts")
seurat_obj
Assays(seurat_obj)
Assays(seurat_obj)
##need to make sure Assays would be "RNA"
##Preprocessing/filtration
gene_filter <- Matrix::rowSums(seurat_obj[["RNA"]]@counts > 0) >= 100
seurat_obj <- seurat_obj[gene_filter, ]
cell_filter <- Matrix::colSums(seurat_obj[["RNA"]]@counts > 0) >= 5
seurat_obj <- seurat_obj[, cell_filter]
seurat_obj
# 4. Optionally: Add spatial coordinates if present
if ("x" %in% colnames(seurat_obj@meta.data) & "y" %in% colnames(seurat_obj@meta.data)) {
  cat("✅ Found spatial coordinates: using 'x' and 'y'\n")
} else {
  cat("⚠️ No spatial coordinates found. Only UMAP will be plotted.\n")
}

# 5. Standard Seurat pipeline
options(future.globals.maxSize = 5 * 1024^3)  # 5GB ##other the code will not run

seurat_obj <- SCTransform(seurat_obj,verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj<- FindNeighbors(seurat_obj,dims = 1:20, k.param = 35)##dims = 1:20/1:30,k.param = >20
seurat_obj <- FindClusters(seurat_obj, resolution = 0.55)#resolution => 0.3
seurat_obj <- RunUMAP(seurat_obj,dims = 1:20)##dims = 1:20/1:30

# 6. Save outputs
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(seurat_obj, file = "seurat_preprocessed_stereo-seq.rds")## file="the stereo-seq data"
seurat_obj<-readRDS("seurat_preprocessed_stereo-seq.rds") ##if you want to recreate the and review the plots

df <- seurat_obj@meta.data
# ggplot UMAP with default colors
# 7. UMAP plot
p1<-DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = "seurat_clusters") + NoLegend()
ggsave("DimPlotseuratstereo-seq.png/tiff", plot = p1, width = 8, height = 6, dpi = 1200)
##spatial plot
df <- seurat_obj@meta.data
# ggplot spatial with default colors
ggplot(df, aes(x = x, y = y, color = seurat_clusters)) +
  geom_point(size = 1.4) +
  coord_fixed() +
  theme_void() +
  labs(title = "Spatial Plot of Seurat_Clusters")
## UMAp and spatial plot for stereo-seq data can be created with custom color  palate

ggsave("Spatial_plot2seuratclusters_stereo-seq.png/tiff", width = 12, height = 10, dpi = 600, bg = "transparent")
##Top 10 or 50 genes
# 8. Marker gene detection
markers <- FindAllMarkers(seurat_obj,only.pos = TRUE)
write.csv(markers, file = "seurat_positive_allMarkers_stereo-seq.csv")

##Tp10 or top 50 genes  can be extracted from seurat 
# 4. Extract gene list
# Top 10 markers per cluster
top10 <- markers_filtered %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)##in case of 50 genes, n=50
gene_list <- unique(top10$gene)
DotPlot(seurat_obj, features = gene_list, group.by = "seurat_clusters") +
  RotatedAxis() +
  theme(text = element_text(size = 8)) +
  labs(title = "Top Marker Genes per Cluster")
ggsave("Dot_plottop10/50_stereo-seq.png/tiff", width = 12, height = 10, dpi = 600, bg = "transparent")
### Cell type annotation with Reference single cell data
# 1. Load annotated scRNA-seq reference (already processed and labeled)
reference <- readRDS("ScRNA-seq.RDS")  # must contain reference$cell_type
reference 
class(seurat_obj)
class(reference)
# 2. Prepare both datasets with SCTransform
reference<- SCTransform(reference, verbose = FALSE)
query <- SCTransform(seurat_obj, verbose = FALSE)  # your Stereo-seq data
#reference$
# 3. Find anchors
anchors <- FindTransferAnchors(reference = reference, query = query, 
                               normalization.method = "SCT", 
                               dims = 1:20)

# 4. Transfer cell-type labels
query <- TransferData(anchorset = anchors, 
                      refdata = reference$CellType, 
                      dims = 1:20, 
                      verbose = FALSE)
table(query$predicted.id)
# 5. Add predictions to your Stereo-seq Seurat object
seurat_obj$predicted_cell_type <- query$predicted.id
saveRDS(seurat_obj, file = "seurat_obj_with_celltypestereo-seq.rds")

## UMAP and spatial plots are generated following  code as  described earlier and saved
