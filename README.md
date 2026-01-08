# SData_adultLiver4815
Spatial transcriptomics preserves molecular and morphological context, enabling direct mapping of gene expression across tissue architecture. Here, we present a high-resolution spatial transcriptomic dataset of bovine liver generated with Stereo-seq (SpaTial Enhanced REsolution Omics-sequencing). 
#Cattle Liver Stereo-seq data 
Raw data: fastq files were processed and analyzed by SAW pipeline; SAW version:v7.0.0 and output as .h5ad or .tissue.gef as gene expression data with spatial co-ordinates. DNA nano ball or spots are organized in different bin; 20,50,100,200 in spatial co-rodinates of liver tissue where cells and gene expression of cells are spotted.

Several steps in Stereo-seq liver data analysis:
# Quality control and Raw data (Fastq files) analysis with SAW pipe;code: Rawdata_saw.sh
# QC, preproceesing, filtering,and clustering, marker genes extrection with Stereo-py pipeline;stereo_seq_Liver_4815_50binFinal.ipynb
# Seural clustering,UMAP clustering , cell type annotation with Seurat pipeline in R
# Cell type annotation with SingleR package in R;
