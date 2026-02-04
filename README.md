# SData_adultLiver4815
Spatial transcriptomics preserves molecular and morphological context, enabling direct mapping of gene expression across tissue architecture. Here, we present a high-resolution spatial transcriptomic dataset of bovine liver generated with Stereo-seq (SpaTial Enhanced REsolution Omics-sequencing). 
#Cattle Liver Stereo-seq data 
Raw data: fastq files were processed and analyzed by SAW pipeline; SAW version:v7.0.0 and output as .h5ad or .tissue.gef as gene expression data with spatial co-ordinates. DNA nano ball or spots are organized in different bins; 20,50,100,200 in spatial co-rodinates of liver tissue where cells and gene expression of cells are spotted.In this project stereo-seq liver data with bin50 was analyzed and generated results. 

Several steps in Stereo-seq liver data analysis:
# Quality control and Raw data (Fastq files) analysis with SAW pipeline:Saw_pipeline_stereo-seq.sh
# QC, preproceesing, filtering,and clustering, marker genes extraction with Stereo-py pipeline; Stereo-py_cattle4815.py & stereo_seq_Liver_4815_50binFinal.ipynb
# Seural clustering,UMAP clustering , cell type annotation with Seurat pipeline in R:StereoSeq_4815Liver-suratclusterMarkergenes_annotattion_Sdata.R
# Cell type annotation with SingleR package in R:Analysis_SingleRannotation_cells_4815liver.R

