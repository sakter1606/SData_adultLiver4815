
# %% [markdown]
# # Stereopy Pipeline (Single Sample, Square Bin or Cell Bin)
# - Reads a SAW output GEF/GEM
# - QC → filtering → normalization
# - HVGs → PCA → neighbors (incl. spatial) → UMAP
# - Clustering (Leiden) → marker genes
# - Optional: interactive high-resolution export, AnnData conversion
#
# Requirements:
#   pip install stereopy
# 
# %%bash
##follow  https://stereopy.readthedocs.io/en/latest/content/00_Installation.html
#it is possible to stereopy with anaconda navigator in windows pc
# pip install --quiet --upgrade stereopy
#Tested with: stereopy >= 1.6.0

# %% [markdown]
# ## Imports & runtime configuration

# %%
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import os
import matplotlib.pyplot as pl
import matplotlib.pyplot as plt
import anndata 
from scipy import io
from scipy.sparse import coo_matrix,csr_matrix

import warnings
warnings.filterwarnings("ignore")

import stereo as st
import os
from pathlib import Path

# --- User paths: update these for your dataset ---
DATA_DIR   = Path("/path/to/data")
# Example files (choose one that exists in your outputs):
TISSUE_GEF = DATA_DIR / "B03022E1.tissue.gef"    # Square-bin tissue GEF from SAW
RAW_GEF    = DATA_DIR / "SS200000135TL_D1.raw.gef"       # Raw square-bin GEF from SAW
CELL_GEF   = DATA_DIR / "SS200000135TL_D1.cellbin.gef"   # Cell-bin GEF from SAW
GEM_FILE   = DATA_DIR / "SS200000135TL_D1.cellbin.gem"   # (optional) GEM

OUT_DIR    = Path("./stpy_out")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Choose input type:
USE_CELLBIN   = False      # set True if you want to read a .cellbin.gef
USE_GEM       = False      # set True if your input is GEM instead of GEF
BIN_SIZE      = 50         # only used for square-bin (bins) data
N_THREADS_IO  = -1         # -1 means all cores for read_gef (bins)
##3) Read data

# %% [markdown]
# ## Read data into StereoExpData

# %%
if USE_GEM and GEM_FILE.exists():
    data = st.io.read_gem(
        file_path=str(GEM_FILE),
        sep="\t",
        bin_type=('cell_bins' if USE_CELLBIN else 'bins'),
        is_sparse=True,
    )
elif USE_CELLBIN and CELL_GEF.exists():
    data = st.io.read_gef(
        file_path=str(CELL_GEF),
        bin_type='cell_bins',
        is_sparse=True,
    )
else:
    # default: square-bin tissue/ raw GEF
    src = str(TISSUE_GEF if TISSUE_GEF.exists() else RAW_GEF)
    data = st.io.read_gef(
        file_path=src,
        bin_type='bins',
        bin_size=BIN_SIZE,
        is_sparse=True,
        num_threads=N_THREADS_IO,
    )
data
##  for one file only
data_path = 'location of the stereo-seq.tissue.gef file'##/project/methylation/liu.yang/data2024/MirxesPackage/....../.tissue.gef
data = st.io.read_gef(file_path=data_path, bin_size=25)
##4) Quality control (QC) & visualization
##QC computes standard metrics (total_counts, n_genes_by_counts, % mito) and provides violin & spatial scatter plots
#
# %% [markdown]
# ## Quality control & overview plots

# %%
data.tl.cal_qc()                # compute QC metrics
data.plt.violin()               # violin plot of QC metrics
data.plt.spatial_scatter()      # spatial scatter of QC metrics
## the code for save each plot 
plt.savefig("plot or figure's name", dpi=600, bbox_inches='tight')
##5) Filtering (cells/bins/genes)
## Normalization & HVGs
##Basic normalization pipeline combines normalize_total and log1p; you can also use scTransform (which auto-finds HVGs). Otherwise, call highly_variable_genes to select HVGs.
##
# %% [markdown]
# ## Normalization & highly variable genes (HVGs)
data.tl.raw_checkpoint()
data.tl.raw
# %%
# Recommended basic normalization:
data.tl.normalize_total(target_sum=10000)
data.tl.log1p()

# Identify HVGs (skip if you used scTransform)
data.tl.highly_variable_genes(
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
    n_top_genes=2000,
    res_key="highly_variable_genes",
)
data.plt.highly_variable_genes(res_key="highly_variable_genes")
##
# %% [markdown]
# ## Normalization & highly variable genes (HVGs)

# %%
# Recommended basic normalization:
data.tl.normalize_total(target_sum=10_000)
data.tl.log1p()

# Identify HVGs (skip if you used scTransform)
data.tl.highly_variable_genes(
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
    n_top_genes=2000,
    res_key="highly_variable_genes",
)
data.plt.highly_variable_genes(res_key="highly_variable_genes")
##
data.tl.key_record
# %% [markdown]
# ## PCA → neighbors → spatial neighbors → UMAP

# %%
data.tl.pca(use_highly_genes=False, n_pcs=30, res_key="pca")
data.plt.elbow(pca_res_key="pca")                      # elbow plot for PCs

data.tl.neighbors(pca_res_key="pca", n_pcs=30, res_key="neighbors")
data.tl.spatial_neighbors(neighbors_res_key="neighbors", res_key="spatial_neighbors")

data.tl.umap(pca_res_key="pca", neighbors_res_key="neighbors", res_key="umap")
data.plt.umap(res_key="umap")
#
# %% [markdown]
# ## Clustering & visualization

# %%
# Leiden using kNN neighbors:
data.tl.leiden(neighbors_res_key="neighbors", res_key="leiden")
data.plt.cluster_scatter(res_key="leiden")          # spatial distribution
data.plt.umap(res_key="umap", cluster_key="leiden") # UMAP colored by clusters

# Optional: also compute spatial-leiden using spatial_neighbors
data.tl.leiden(neighbors_res_key="spatial_neighbors", res_key="spatial_leiden")
data.plt.cluster_scatter(res_key="spatial_leiden")
##9) Marker genes

# %% [markdown]
# ## Marker genes (DE testing across clusters)

# %%
print(data.raw) 
data.tl.raw_checkpoint()
data.tl.find_marker_genes(
    cluster_res_key="leiden",
    method="t_test",
    use_highly_genes=False,
    use_raw=True,
)
data.plt.marker_genes_text(res_key="marker_genes", markers_num=10, sort_key="scores")
data.plt.marker_genes_scatter(res_key="marker_genes", markers_num=5)
data_sub = st.io.stereo_to_anndata(sub, flavor="scanpy", output=str(sub_file))
# Compute QC metrics
adata.obs['nCount_RNA'] = adata.X.sum(axis=1).A1 if hasattr(adata.X, "A1") else adata.X.sum(axis=1)
adata.obs['nFeature_RNA'] = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, "A1") else (adata.X > 0).sum(axis=1)
## SCipy  pipeline
# Compute mitochondrial gene percentage (genes starting with 'MT-')
## read data
adata = sc.read_h5ad("4815_Refannotation.h5ad")
sc.pp.filter_cells(adata, min_genes=10)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
sc.pl.violin(
    adata,
    ["nFeature_RNA", "nCount_RNA", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,

qc_metrics = ["nFeature_RNA", "nCount_RNA", "pct_counts_mt"]

plt.figure(figsize=(10,4))
for i, metric in enumerate(qc_metrics):
    plt.subplot(1, 3, i+1)
    sns.violinplot(y=adata.obs[metric], color="#00AFC1", inner="point", linewidth=0.8)
    plt.title(metric)
    plt.xlabel("")
    plt.ylabel("")
plt.tight_layout()
plt.show()
adata = adata[
    (adata.obs["nFeature_RNA"] > 500) &
    (adata.obs["nFeature_RNA"] < 7500) &
    (adata.obs["pct_counts_mt"] < 25),
    :
]
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)
## creating violin plot  of marker genes for 5 cell types
cell_types_to_keep = [
    'Hepatocytes',
    'Cholangiocytes',
    'Kupffer cells',
    'Liver sinusoidal endothelial cells',
    'Neutrophils'
]

adata1 = adata[adata.obs['cell_type'].isin(cell_types_to_keep)].copy()
adata1.obs['cell_type']
marker_genes = [
   "APOA2"
]
# Check that the genes exist in your data
existing_genes = [g for g in marker_genes if g in adata.var_names]
print(f"✅ Found {len(existing_genes)} marker genes: {existing_genes}")

# Create violin plot by cell type
sc.pl.violin(
    adata1,
    keys=existing_genes,
    groupby="cell_type",        # change to your annotation column
    stripplot=False,           # removes jittered points
    multi_panel=False,          # plot each gene in its own row
    rotation=90,
    jitter=0,
    show=True,
    figsize=(12,10)
  ##  cell type annotation in spatial location
  data.cells['cellType']=adata.obs['cell_type']
data.plt.cluster_scatter(res_key='cellType',dot_size=8,marker='o')
plt.gcf().set_size_inches(16, 10)  # width x height in inches
plt.tight_layout()
plt.savefig("spatial_plot_4815Liver.png", dpi=600, bbox_inches='tight')

plt.show()

##  all code is run in Jupyter notebook,each plot is saved as .jpg or .tif  from the output
