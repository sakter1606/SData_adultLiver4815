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
DefaultAssay(seurat_obj) <- "originalexp"
seurat_obj[["RNA"]] <- seurat_obj[["originalexp"]]
#seurat_obj[["originalexp"]] <- NULL
DefaultAssay(seurat_obj) <- "RNA"
##Preprocessing/filtration
gene_filter <- Matrix::rowSums(seurat_obj[["RNA"]]@counts > 0) >= 150
seurat_obj <- seurat_obj[gene_filter, ]
cell_filter <- Matrix::colSums(seurat_obj[["RNA"]]@counts > 0) >= 10
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
seurat_obj<- FindNeighbors(seurat_obj,dims = 1:20, k.param = 30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.55)
seurat_obj <- RunUMAP(seurat_obj,dims = 1:20)
#p <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = "seurat_clusters") + NoLegend()
#p ##if you want labels in the graph
df <- seurat_obj@meta.data
# ggplot UMAP with manual colors
ggplot(df, aes(x = x, y = y, color = seurat_clusters)) +
  geom_point(size = 1.4, alpha = 0.8) +
  scale_color_manual(values = my_colors_seurat) +
  theme_void(base_size = 14) +
  ggtitle("Seurat Clusters with Custom Colors")

saveRDS(seurat_obj, file = "seurat_preprocessed_B03022E1_4815Liver_150_10genefilteredParam30res0.55seurat.rds")
# 6. Save outputs
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
# 7. UMAP plot
seurat_obj<-readRDS("seurat_preprocessed_B03022E1_4815Liver_150_10genefilteredParam30res0.55seurat.rds")
my_colors_seurat <- c(
  "0" = "#205444", "1" = "#ff7f0e", "2" = "#2ca02c", "3" = "#c2a515",
  "4" ="#e0d643" , "5" = "#bf454f",
  "6" = "#661f45","7" = "#7d12c9","8"="#1aba7d",
  "9" = "#79b7c7","10" = "#b6e043",
  "11"= "#70580d", "12"="#a38577","13"="#3043bf","14"= "#e0ba7b"
  # Add al cell types here with distinct colors
)

# Show plot
#print(p1)
p1 <- DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE,
  group.by = "seurat_clusters",
  #group.by = "predicted_cell_type",
  cols = my_colors_seurat
) + NoLegend() #with labels in the graph

p1
p1 <- DimPlot(
  seurat_obj,
  reduction = "umap",
  label = FALSE,
  group.by = "seurat_clusters",
  #group.by = "predicted_cell_type",
  cols = my_colors_seurat
) 
p1
ggsave("DimPlotseuratNewCLB03022E1_4815Liver_150_10genefilteredParam30res0.55.png", plot = p1, width = 8, height = 6, dpi = 1200)

# 8. Marker gene detection
markers <- FindAllMarkers(seurat_obj)## both positive and negative expression
write.csv(markers, file = "seurat_allMarkersB03022E1_4815Liver_150_10genefilteredParam30res0.55seurat_09_10_2025.csv")

markers <- FindAllMarkers(seurat_obj,only.pos = TRUE)
write.csv(markers, file = "seurat_positivemarkersB03022E1_4815Liver_150_10genefilteredParam30res0.55_09_10_2025.csv")


df <- seurat_obj@meta.data
# ggplot UMAP with manual colors
ggplot(df, aes(x = x, y = y, color = seurat_clusters)) +
  geom_point(size = 1.4, alpha = 0.8) +
  scale_color_manual(values = my_colors_seurat) +
  theme_void(base_size = 14) +
  ggtitle("Seurat Clusters with Custom Colors")

ggsave("Spatial_plot2seuratclusters_B03022E1_4815Liver_150_10genefilteredParam30res0.55.png", width = 12, height = 10, dpi = 600, bg = "transparent")
##Top 10 or 50 genes
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.1,logfc.threshold = 0.25)
markers_filtered <- markers %>% filter(p_val_adj <= 0.05)
top50 <- markers_filtered %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)
write.csv(top50, file = "top50_B03022E1_4815Liver_150_10genefilteredParam30res0.55.csv")
top10 <- markers_filtered %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
write.csv(top10, file = "top10_B03022E1_4815Liver_150_10genefilteredParam30res0.55.csv")
# 4. Extract gene list
gene_list <- unique(top10$gene)

DotPlot(seurat_obj, features = gene_list, group.by = "seurat_clusters") +
  RotatedAxis() +
  theme(text = element_text(size = 8)) +
  labs(title = "Top Marker Genes per Cluster")
ggsave("Dot_plottop10seuratclusters_B03022E1_4815Liver_150_10genefilteredParam30res0.55.png", width = 12, height = 10, dpi = 600, bg = "transparent")
### Cell type annotation with Reference single cell data
# 1. Load annotated scRNA-seq reference (already processed and labeled)
reference <- readRDS("ScRef_Liver_IRIS.RDS")  # must contain reference$cell_type
#seurat_obj <- readRDS("seurat_preprocessed_B03022E1_4815Liver_150_10genefilteredParam30res0.55seurat.rds")
#assayNames(reference)
#obj<-readRDS("seurat_obj_with_celltypeB03022E1_4815Liver_150_10genefilteredParam30res0.55seurat.rds")
reference 
#logcounts(reference) <- log1p(assay(reference, "counts"))
#reference <- as.Seurat(reference, counts = "X", data = "logcounts")
#DefaultAssay(reference_seurat) <- "RNA"
class(seurat_obj)
class(reference)
# 2. Prepare both datasets with SCTransform
reference<- SCTransform(reference, verbose = FALSE)
#seurat_obj$predicted_cell_type
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
saveRDS(seurat_obj, file = "seurat_obj_with_celltypeB03022E1_4815Liver_150_10genefilteredParam30res0.55seurat.rds")
obj<-readRDS("seurat_obj_with_celltypeB03022E1_4815Liver_150_10genefilteredParam30res0.55seurat.rds")

my_colors_dim <- c(
  "Kupffer cells" = "#2ca02c","Cholangiocytes" = "#7d12c9",
  "Hepatocytes" = "#79b7c7",
  "Liver sinusoidal endothelial cells"="#ff7f0e","Neutrophils"="#c2a515"
  # Add al cell types here with distinct colors only for 5 cell types
)
my_colors_dim <- c(
  "B cells" = "#205444", "T cells" = "#ff7f0e", "NK cells" = "#2ca02c", "NK T cells" = "#c2a515",
  "Dendritic cells" ="#e0d643" , "Myofibroblasts" = "#bf454f",
  "Kupffer cells" = "#661f45","Cholangiocytes" = "#7d12c9","Plasma cells"="#cc72c8",
  "Hepatocytes" = "#79b7c7","Unknown cells" = "#b6e043",
  "Blood vascular endothelial cells"= "#70580d", "Liver sinusoidal endothelial cells"="#a38577","Neutrophils"="#3043bf","Proliferative dendritic cells"= "#e0ba7b"
  # Add al cell types here with distinct colors
)
p2 <- DimPlot(
  obj,
  reduction = "umap",
  label = FALSE,
  #group.by = "seurat_clusters",
  group.by = "predicted_cell_type",
  cols = my_colors_dim
) 
p2
p2 <- DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE,
  #group.by = "seurat_clusters",
  group.by = "predicted_cell_type",
  cols = my_colors
) + NoLegend()
p2
#ggsave("DimPlot_5celltypeannotatedB03022E1_4815Liver_150_10genefilteredParam30res0.55.tif", plot = p2, width = 8, height = 6, dpi = 1200)
ggsave("DimPlot_NewcelltypeannotatedB03022E1_4815Liver_150_10genefilteredParam30res0.55.tif", plot = p2, width = 8, height = 6, dpi = 1200)
ggplot(seurat_obj@meta.data, aes(x = x, y = y, color = predicted_cell_type)) +
  geom_point(size = 1.2, alpha = 1.0) +
  coord_fixed() +  # 👈 Coordinate setting goes here
  theme_void(base_size = 14) +  # 👈 Then apply the theme
  scale_color_manual(values = my_colors) +
  labs(color = "Cell Type") +
  ggtitle("Spatial Distribution of Cell Types") +
  
  theme(
    plot.background = element_rect(fill = "White", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
ggsave("Spatialplot1Finalseurat_celltypeannotated_B03022E1_4815Liver_150_10genefilteredParam30res0.55seurat.tiff", width = 12, height = 10, dpi = 1200, bg = "transparent")
