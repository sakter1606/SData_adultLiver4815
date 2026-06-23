## SData_adultLiver4815

Spatial transcriptomics preserves both molecular and morphological context, enabling the direct mapping of gene expression across tissue architecture. In this project, we analyzed a high-resolution spatial transcriptomics dataset of bovine adult liver generated using Stereo-seq, also known as SpaTial Enhanced REsolution Omics-sequencing.

## Cattle Liver Stereo-seq Data

Raw FASTQ files were processed and analyzed using the SAW pipeline, version 7.0.0. The processed outputs were generated as `.h5ad` or `.tissue.gef` files, which contain gene expression data with spatial coordinates.

In Stereo-seq, DNA nanoballs or spatial spots are organized into different bin sizes, such as bin20, bin50, bin100, and bin200, based on the spatial coordinates of the liver tissue. These bins capture gene expression signals across the tissue and allow visualization of spatial gene expression patterns. In this project, the bovine liver Stereo-seq data were analyzed at bin50 resolution.
Analysis Workflow:
FASTQ → SAW v7.0.0 → GEF/H5AD → Stereo-py QC/clustering i)→ SingleR annotation 
                                                        ii)→ Seurat UMAP/markers → Seurat cell type annotation 
 

## Steps in the Stereo-seq Liver Data Analysis        
1. **Quality control and raw data processing using the SAW pipeline**
   Script: `Saw_pipeline_stereo-seq.sh`

2. **Quality control, preprocessing, filtering, clustering, and marker gene extraction using Stereo-py**
   Scripts/notebooks: `Stereo-py_cattle4815.py` and `stereo_seq_Liver_4815_50binFinal.ipynb`

3. **Seurat clustering, UMAP visualization, marker gene identification, and cell type annotation in R**
   Script: `StereoSeq_4815Liver-suratclusterMarkergenes_annotattion_Sdata.R`

4. **Cell type annotation using the SingleR package in R**
   Script: `Analysis_SingleRannotation_cells_4815liver.R`
