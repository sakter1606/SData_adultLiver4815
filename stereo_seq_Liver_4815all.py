#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import os
import matplotlib.pyplot as pl
import igraph
#import scvelo as scv
#import loonpy as lnp
import anndata 
from scipy import io
from scipy.sparse import coo_matrix,csr_matrix


# In[2]:


import stereo as st
import warnings
warnings.filterwarnings('ignore')


# In[21]:


import matplotlib.pyplot as plt
import stereo as st

# --- Publication Standard Global Settings ---
plt.rcParams.update({
    'figure.dpi': 300,              # High resolution for viewing in notebooks
    'savefig.dpi': 600,             # Ultra-high resolution for saving
    'font.size': 24,                # Base font size
    'font.weight': 'bold',          # Make text bold
    'axes.labelsize': 24,           # X and Y axis label size
    'axes.labelweight': 'bold',     # Make X and Y axis labels bold
    'axes.titlesize': 24,           # Title size
    'axes.titleweight': 'bold',     # Title bold
    'xtick.labelsize': 22,          # X tick numbers
    #'ytick.fontstyle : 'italic',      # Y tick numbers
    'legend.fontsize': 22,          # Legend text size
    'savefig.bbox': 'tight',        # Removes empty white space around the image
    'savefig.format': 'tif'         # Default to saving as vector format
})


# In[ ]:


# Make border/axis box thicker
for spine in ax.spines.values():
    spine.set_linewidth(2.0)
    spine.set_color("black")


# In[3]:


adata = sc.read_h5ad("B03022E1_bin50out.h5ad")


# In[4]:


sc.pl.highest_expr_genes(adata, n_top=25)


# In[238]:


adata = sc.read_h5ad("CattleLiver4815_seurat2_bin150_08_17_2026.h5ad")


# In[240]:


sc.pl.highest_expr_genes(adata, n_top=25)


# In[13]:


import matplotlib.pyplot as plt
import stereo as st
plt.rcParams['axes.linewidth'] = 2.5  
plt.rcParams['axes.edgecolor'] = 'black'


# In[15]:


#data = st.io.read_gef('B03022E1.tissue.gef')


# In[7]:


data_path = './B03022E1.tissue.gef'
st.io.read_gef_info(data_path)


# In[23]:


data = st.io.read_gef(file_path=data_path, bin_size=100)


# In[ ]:


data


# In[39]:


data = st.io.read_gef(file_path=data_path, bin_size=50)


# In[28]:


data = st.io.read_gef(file_path=data_path, bin_size=150)


# In[33]:


data = st.io.read_gef(file_path=data_path, bin_size=200)


# In[42]:


data


# In[43]:


data.tl.cal_qc()


# In[44]:


data.plt.violin()


# In[45]:


data.plt.spatial_scatter()


# In[ ]:


data.plt.violin()


# In[ ]:





# In[ ]:





# In[ ]:





# In[37]:


data.plt.violin()


# In[31]:


data.plt.violin()


# In[32]:


data.plt.spatial_scatter()


# In[26]:


data.plt.violin()


# In[27]:


data.plt.spatial_scatter()


# In[11]:


data.plt.violin()


# In[38]:


data.plt.spatial_scatter()


# In[61]:


data.plt.spatial_scatter()


# In[249]:


data.plt.genes_count()


# In[13]:


data.plt.genes_count()


# In[ ]:


data.tl.filter_cells(
    min_gene=3,  # Minimum number of genes per cell
    max_gene=None,  # Optional: Maximum number of genes per cell
    mincounts=20,  # Minimum UMI counts per cell
    maxcounts=None,  # Optional: Maximum UMI counts per cell
    pctcountsmt=0  # Mitochondrial gene percentage threshold
)


# In[21]:


data.tl.filter_cells(
        min_counts=20,
        min_genes=3,
        pct_counts_mt=5,
        inplace=True
        )
data


# In[250]:


data.tl.filter_cells(
        min_counts=25,
        min_genes=10,
        pct_counts_mt=5,
        inplace=True
        )
data


# In[251]:


data


# In[252]:


data.tl.raw_checkpoint()


# In[253]:


data.tl.raw


# In[254]:


##Normalization


# In[ ]:





# In[255]:


data.tl.normalize_total(target_sum=10000)
data.tl.log1p()


# In[256]:


data.tl.highly_variable_genes(
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_top_genes=2000,
        res_key='highly_variable_genes'
        )


# In[257]:


data.plt.highly_variable_genes(res_key='highly_variable_genes')


# In[260]:


data.tl.pca(
        use_highly_genes=False,
        n_pcs=30,
        res_key='pca'
        )


# In[261]:


data.tl.key_record


# In[262]:


data.plt.elbow(pca_res_key='pca')


# In[278]:


##Neighborhood graph
data.tl.neighbors(
        pca_res_key='pca',
        n_pcs=30,
        res_key='neighbors'
        )


# In[274]:


# compute spatial neighbors
data.tl.spatial_neighbors(
        neighbors_res_key='neighbors',
        res_key='spatial_neighbors'
        )


# In[279]:


#
data.tl.umap(pca_res_key='pca', neighbors_res_key='neighbors', res_key='umap')


# In[ ]:


import numpy as np
array = np.memmap('large_array.dat', dtype='int64', mode='w+', shape=(2602689504, 2))


# In[280]:


##Clustering
data.tl.leiden(neighbors_res_key='neighbors', res_key='leiden')


# In[281]:


data.plt.cluster_scatter(res_key='leiden',dot_size=20,marker='o')


# In[220]:


data.plt.cluster_scatter(res_key='leiden',dot_size=20,marker='o')


# In[167]:


data.plt.cluster_scatter(res_key='leiden',dot_size=20,marker='o')


# In[123]:


data.plt.cluster_scatter(res_key='leiden',dot_size=20,marker='o')


# In[78]:


data.plt.cluster_scatter(res_key='leiden',dot_size=20,marker='o')


# In[29]:


data.plt.cluster_scatter(res_key='leiden',dot_size=20,marker='o')


# In[35]:


data.plt.cluster_scatter(res_key='leiden',dot_size=20,marker='o')


# In[31]:


data.plt.cluster_scatter(res_key='leiden',dot_size=20,marker='o')


# In[37]:


data.plt.cluster_scatter(res_key='leiden', groups=['1', '2'],dot_size=20,marker='o')


# In[168]:


data.plt.umap(res_key='umap', cluster_key='leiden',dot_size=20,marker='o')


# In[221]:


data.plt.umap(res_key='umap', cluster_key='leiden',dot_size=20,marker='o')


# In[79]:


data.plt.umap(res_key='umap', cluster_key='leiden',dot_size=20,marker='o')


# In[39]:


data.plt.umap(res_key='umap', cluster_key='leiden',dot_size=20,marker='o')


# In[125]:


data.plt.umap(gene_names=['ALB','ORM1','ORM1','APOA2','FGA','APOE','C3','AMBP','SELENOP','TF','FGB','CFB','CAT','VTN','FGG','CYTB'], res_key='umap') 


# In[33]:


data.plt.umap(gene_names=['ALB','ORM1','ORM1','APOA2','FGA','APOE','C3','AMBP','SELENOP','TF','FGB','CFB','CAT','VTN','FGG','CYTB'], res_key='umap') 


# In[223]:


data.tl.leiden(neighbors_res_key='spatial_neighbors', res_key='spatial_leiden')


# In[170]:


data.plt.cluster_scatter(res_key='spatial_leiden',dot_size=20,marker='o')


# In[82]:


data.plt.cluster_scatter(res_key='spatial_leiden',dot_size=20,marker='o')


# In[224]:


data.plt.cluster_scatter(res_key='spatial_leiden',dot_size=20,marker='o')


# In[225]:


data.tl.louvain(neighbors_res_key='neighbors', res_key='louvain')


# In[226]:


data.plt.cluster_scatter(res_key='louvain',dot_size=20,marker='o')


# In[130]:


data.plt.cluster_scatter(res_key='louvain',dot_size=20,marker='o')


# In[172]:


data.plt.cluster_scatter(res_key='louvain',dot_size=20,marker='o')


# In[41]:


data.plt.cluster_scatter(res_key='louvain',dot_size=20,marker='o')


# In[131]:


data.plt.interact_cluster(res_key='leiden')


# In[227]:


print(data.raw) 


# In[228]:


data.tl.find_marker_genes(
        cluster_res_key='leiden',
        method='t_test',
        use_highly_genes=False,
        use_raw=True
        )


# In[229]:


data.plt.marker_genes_text(
        res_key='marker_genes',
        markers_num=10,
        sort_key='scores'
        )


# In[ ]:





# In[230]:


data.plt.marker_genes_scatter(res_key='marker_genes', markers_num=5)


# In[ ]:


data.plt.marker_genes_scatter(res_key='marker_genes', markers_num=5)


# In[ ]:





# In[89]:


data.plt.marker_genes_scatter(res_key='marker_genes', markers_num=5)


# In[177]:


data.plt.marker_gene_volcano(group_name='2.vs.rest', vlines=False,dot_size=30)


# In[90]:


data.plt.marker_genes_scatter(res_key='marker_genes', markers_num=5)


# In[60]:


data.plt.marker_genes_scatter(res_key='marker_genes', markers_num=10)


# In[ ]:





# In[47]:


data.plt.marker_gene_volcano(group_name='3.vs.rest', vlines=False)


# In[231]:


data.tl.filter_marker_genes(
    marker_genes_res_key='marker_genes',
    min_fold_change=1,
    min_in_group_fraction=0.25,
    max_out_group_fraction=0.5,
    res_key='marker_genes_filtered'
)


# In[232]:


##Annotation
data.plt.interact_annotation_cluster(
            res_cluster_key='leiden',
            res_marker_gene_key='marker_genes',
            res_key='leiden_annotation'
            )


# In[233]:


##Annotation
data.plt.interact_annotation_cluster(
            res_cluster_key='leiden',
            res_marker_gene_key='marker_genes',
            res_key='leiden_annotation'
            )


# In[96]:


data.plt.marker_genes_scatter(res_key='marker_genes', markers_num=5)


# In[234]:


deg_results = data.tl.result['marker_genes']


# In[235]:


deg_results


# In[184]:


deg_results = data.tl.result['marker_genes']


# In[ ]:


deg_results


# In[236]:


import pandas as pd

# Assuming you already created this
df = pd.DataFrame({
    'sample_id': list(range(1,9)),
    'data': [
        deg_results["1.vs.rest"], deg_results["2.vs.rest"], deg_results["3.vs.rest"],
        deg_results["4.vs.rest"], deg_results["5.vs.rest"], deg_results["6.vs.rest"],
        deg_results["7.vs.rest"], deg_results["8.vs.rest"]
    ]
})
# ✅ Loop through each row and save as individual CSVs
for index, row in df.iterrows():
    sample_id = row['sample_id']
    sample_df = row['data']
    
    # Optional: Add a sample_id column to each sub-DataFrame
    sample_df = sample_df.copy()
    sample_df['sample_id'] = sample_id

    # Save to CSV
    filename = f"Cattleliver4815_{sample_id}_deg_resultsLeiden_bin150.csv"
    sample_df.to_csv(filename, index=False)

    print(f"Saved: {filename}")


# In[ ]:


import pandas as pd

# Assuming you already created this
df = pd.DataFrame({
    'sample_id': list(range(1,21)),
    'data': [
        deg_results["1.vs.rest"], deg_results["2.vs.rest"], deg_results["3.vs.rest"],
        deg_results["4.vs.rest"], deg_results["5.vs.rest"], deg_results["6.vs.rest"],
        deg_results["7.vs.rest"], deg_results["8.vs.rest"], deg_results["9.vs.rest"],
        deg_results["10.vs.rest"], deg_results["11.vs.rest"], deg_results["12.vs.rest"],
        deg_results["13.vs.rest"], deg_results["14.vs.rest"],deg_results["15.vs.rest"],
        deg_results["16.vs.rest"],deg_results["17.vs.rest"],deg_results["18.vs.rest"],
        deg_results["19.vs.rest"],deg_results["20.vs.rest"]
    ]
})
# ✅ Loop through each row and save as individual CSVs
for index, row in df.iterrows():
    sample_id = row['sample_id']
    sample_df = row['data']
    
    # Optional: Add a sample_id column to each sub-DataFrame
    sample_df = sample_df.copy()
    sample_df['sample_id'] = sample_id

    # Save to CSV
    filename = f"Cattleliver4815_{sample_id}_deg_results.csv"
    sample_df.to_csv(filename, index=False)

    print(f"Saved: {filename}")


# In[80]:


def classify_gene(log2fc, pval):
    if log2fc > 1 and pval < 0.05:
        return "Upregulated"
    elif log2fc < -1 and pval < 0.05:
        return "Downregulated"
    else:
        return "Not Significant"

# Example usage
print(classify_gene(1.5, 0.01))  # Upregulated
print(classify_gene(-2.0, 0.04)) # Downregulated
print(classify_gene(0.5, 0.2))   # Not Significant


# In[143]:


adata = st.io.stereo_to_anndata(data, output='CattleLiver4815_seurat2_bin200_08_17_2026.h5ad')


# In[9]:


import pandas as pd
import glob


# In[10]:


all_files = glob.glob("Cattleliver4815_*_deg_results.csv")
df_list = [pd.read_csv(file) for file in all_files]


# In[11]:


df_list


# In[96]:


df = pd.DataFrame({
    'Sample_id': list(range(1,21))})


# In[ ]:


obj={}
def classify_gene(log2fc, pval):
    if log2fc > 1 and pval < 0.05:
        return "Upregulated"
    elif log2fc < -1 and pval < 0.05:
        return "Downregulated"
    else:
        return "Not Significant"

# Example usage
print(classify_gene(1.5, 0.01))  # Upregulated
print(classify_gene(-2.0, 0.04)) # Downregulated
print(classify_gene(0.5, 0.2))   # Not Significant
for i in range(0,12):
    d1=df['Sample_id'][i]
    print(d1)
    d=df_list[i]
    obj[d1]=[]
# Extract upregulated gene names
    upregulated_genes = d.loc[(d["log2fc"] > 1) & (d["pvalues_adj"] < 0.05), "genes"]

# Extract downregulated gene names
    downregulated_genes = d.loc[(d["log2fc"] < -1) & (d["pvalues_adj"] < 0.05), "genes"]

# Print first few gene names
    print("Upregulated Genes:", upregulated_genes.head())
    print("Downregulated Genes:", downregulated_genes.head())
       
    s={'DEG_up':upregulated_genes,'DEG_down':downregulated_genes}
    df2=pd.DataFrame(data=s)
    with open(f"{d1}_sampleDEGS.csv", 'a') as f:
        (df2).to_csv(f, header=True) 
            


# In[ ]:





# In[ ]:


obj={}
def classify_gene(log2fc, pval):
    if log2fc > 1 and pval < 0.05:
        return "Upregulated"
    elif log2fc < -1 and pval < 0.05:
        return "Downregulated"
    else:
        return "Not Significant"

# Example usage
print(classify_gene(1.5, 0.01))  # Upregulated
print(classify_gene(-2.0, 0.04)) # Downregulated
print(classify_gene(0.5, 0.2))   # Not Significant
for i in range(0,20):
    d1=df['Sample_id'][i]
    print(d1)
    d=df_list[i]
    obj[d1]=[]
# Extract upregulated gene names
    d["class"] = d.apply(lambda row: classify_gene(row["log2fc"], row["pvalues_adj"]), axis=1)

    print(d)
# Extract downregulated gene names
    #downregulated_genes = d.loc[(d["log2fc"] < -1) & (d["pvalues_adj"] < 0.05), "genes"]
    counts = d["class"].value_counts()
    print(counts)

# Print first few gene names
    #print("Upregulated Genes:", upregulated_genes.head())
    #print("Downregulated Genes:", downregulated_genes.head())
       
    #s={'DEG_up':upregulated_genes,'DEG_down':downregulated_genes}
    #df2=pd.DataFrame(data=d)
    #with open(f"{d1}_4815LiverstereseqDEGS.csv", 'a') as f:
        #(df2).to_csv(f, header=True) 
            


# In[ ]:


df1["class"] = df1.apply(lambda row: classify_gene(row["log2fc"], row["pvalues_adj"]), axis=1)

print(df1)


# In[ ]:





# In[75]:


deg_key = '4.vs.rest'  # Replace with the key you need
deg_data = deg_results[deg_key]  # Access its value

print(type(deg_data))  # Check its type (likely dict or DataFrame)
print(deg_data)  # Print content


# In[64]:


print(type(deg_results))  # Should be <class 'dict'>
print(deg_results.keys())


# In[125]:


d1=deg_results["2.vs.rest"]


# In[72]:


d2=deg_results["1.vs.rest"]


# In[101]:


d3=deg_results["3.vs.rest"]


# In[73]:


d2


# In[74]:


import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import os
import matplotlib.pyplot as pl
import igraph
#import scvelo as scv
#import loonpy as lnp
import anndata
from scipy import io
from scipy.sparse import coo_matrix,csr_matrix


# In[126]:


deg_df = pd.DataFrame(d1)  # Convert to DataFrame
deg_df.to_csv(f"2.vs.rest_DEGs_4815_liverbin_200.csv", index=False)  # Save as CSV


# In[97]:


deg_df = pd.DataFrame(d2)  # Convert to DataFrame
deg_df.to_csv(f"1.vs.rest_DEGs_4815_liverbin200.csv", index=False)  # Save as CSV


# In[102]:


deg_df = pd.DataFrame(d3)  # Convert to DataFrame
deg_df.to_csv(f"3.vs.rest_DEGs_4815_liverbin200.csv", index=False)  # Save as CSV


# In[116]:


df1=pd.read_csv("1.vs.rest_DEGs_4815_liverbin200.csv")


# In[117]:


def classify_gene(log2fc, pval):
    if log2fc > 1 and pval < 0.05:
        return "Upregulated"
    elif log2fc < -1 and pval < 0.05:
        return "Downregulated"
    else:
        return "Not Significant"


# In[118]:


# Example usage
print(classify_gene(1.5, 0.01))  # Upregulated
print(classify_gene(-2.0, 0.04)) # Downregulated
print(classify_gene(0.5, 0.2))   # Not Significant


# In[119]:


# Apply classification
df1["class"] = df1.apply(lambda row: classify_gene(row["log2fc"], row["pvalues_adj"]), axis=1)

print(df1)


# In[ ]:


#import numpy as np

#conditions = [
    #(df1["log2FoldChange"] > 1) & (df1["pvals_adj"] < 0.05),  # Upregulated
    #(df1["log2FoldChange"] < -1) & (df1["pvals_adj"] < 0.05)  # Downregulated

#choices = ["Upregulated", "Downregulated"]

#df1["class"] = np.select(conditions, choices, default="Not Significant")
#print(df)


# In[120]:


counts = df1["class"].value_counts()
print(counts)


# In[122]:


df1.to_csv(f"degs_4815Liver_1vsrestbin200.csv", index=False)  # Save as CSV


# In[128]:


import pandas as pd

# Load DEGs file (modify path accordingly)
df = pd.read_csv("2.vs.rest_DEGs_4815_liverbin_200.csv")


# In[129]:


df


# In[130]:


from gprofiler import GProfiler
gp = GProfiler(return_dataframe=True)
results = gp.profile(organism="btaurus", query=df["genes"].tolist())
# Display results
print(results.head())


# In[131]:


df2= pd.DataFrame(results)  # Convert to DataFrame
df2.to_csv(f"funation_ana_4815Liver_2vsrestbin200.csv", index=False)  # Save as CSV


# In[ ]:


##https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html#


# In[103]:


import pandas as pd
import scanpy as sc


# In[99]:


adata = sc.read_h5ad("B03022E1_bin50out.h5ad")


# In[ ]:





# In[100]:


adata


# In[102]:


sc.pl.highest_expr_genes(adata, n_top=50)


# In[143]:


results_file = "liver_B03022E1_bin50out.h5ad"


# In[106]:


sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)


# In[107]:


adata


# In[109]:


##annotate the group of mitochondrial genes as "mt"
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
   adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)


# In[110]:


sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)


# In[111]:


sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")


# In[112]:


adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()


# In[113]:


sc.pp.log1p(adata)


# In[114]:


sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


# In[115]:


sc.pl.highly_variable_genes(adata)


# In[116]:


adata.raw = adata.copy()


# In[117]:


adata = adata[:, adata.var.highly_variable]


# In[118]:


sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])


# In[119]:


sc.pp.scale(adata, max_value=10)


# In[120]:


sc.tl.pca(adata, svd_solver="arpack")


# In[129]:


sc.pl.pca(adata, color="ALB")


# In[125]:


sc.pl.pca(adata, color="APOA2")


# In[133]:


sc.pl.pca(adata, color="FGA")


# In[134]:


sc.pl.pca_variance_ratio(adata, log=True)


# In[144]:


adata.write(results_file)


# In[145]:


adata


# In[141]:


sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)


# In[146]:


sc.tl.louvain(adata) 


# In[148]:


sc.tl.tsne(adata)


# In[147]:


sc.tl.leiden(adata)


# In[149]:


sc.tl.paga(adata)


# In[150]:


sc.tl.paga(adata)
sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata, init_pos='paga')


# In[153]:


sc.tl.umap(adata)


# In[163]:


sc.pl.umap(adata, color=['ALB','ORM1','ORM1','APOA2','FGA','APOE','C3','AMBP','SELENOP','TF','FGB','CFB','CAT','VTN','FGG','CYTB'])


# In[ ]:


sc.pl.umap(adata, color=['ALB','ORM1','FGA','AMBP','SELENOP','TF','FGB','CFB','CAT','VTN','FGG','CYTB'], use_raw=False)


# In[179]:


adata.write(results_file)

