#!/usr/bin/env python
# coding: utf-8

# <font size="8">Dataset filtering and annotations</font>

# # Environment

# In[ ]:


import ipynbname
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy.external as sce
from datetime import datetime
from gprofiler import GProfiler
import decoupler as dc
import triku as tk
from gprofiler import GProfiler
import os


# In[ ]:


sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100)
sc.logging.print_header()


# In[ ]:


def selectMarkers(adataObj, mList):
    """  From a list of gene names select only the genes that are present in adata.var
    """
    #Select markers present in adata
    p = adataObj.var_names[adataObj.var_names.isin(mList) == True]
    #Keep the same order as input list
    p = [x for x in mList if x in p]   
    
    #Select missing genes
    ab = set(mList).difference(set(adataObj.var_names))
    
    #Print message SISTEMA
    if len(ab) == len(mList):
        print('\nAll markers are missing')
    else:
        print('\nThe following marker genes are missing: ', ab)
        
    return(p)


# In[ ]:


def CustomGO(adata, cluster, rank, n_markers=40,  show=10):
    """  
        GO analysis with GProfiler for cluster top-marker genes. Adapted for toxo gene names.
    """
    
    GroupMarkers = pd.DataFrame(adata.uns[rank]['names']).head(n_markers)   
    q = GroupMarkers[cluster].str.replace('cellranger_gex-GRCh38-2020-A_', '').tolist()
    u = adata.var_names.str.replace('cellranger_gex-GRCh38-2020-A_', '').tolist()
    return gp.profile(organism='hsapiens', sources=['GO:BP', 'GO:CC'], query=q, 
           background=u, no_iea=True).head(show)


# In[ ]:


def CustomUmap(adata, genes):
    #genes = ['cellranger_gex-GRCh38-2020-A_' + gene for gene in genes]
    genes = selectMarkers(adata, genes)
    sc.pl.umap(adata, color=genes, size=10, frameon=False,
               vmin='p1',  vmax='p99')


# In[ ]:


def CustomDA(adata, genes):
    #genes = ['cellranger_gex-GRCh38-2020-A_' + gene for gene in genes]
    genes = selectMarkers(adata, genes)
    sc.pl.draw_graph(adata, color=genes, size=10, frameon=False,
               vmin='p1',  vmax='p99')


# In[ ]:


data_folder = "./data"
result_folder = './results'


# In[ ]:


print(datetime.now())


# # Data load and concatenation

# In[ ]:


adata1 = sc.read_10x_mtx(
    os.path.join(data_folder, "ctl_raw_feature_bc_matrix"),
    var_names='gene_symbols',
    cache=True)   
adata1.obs['batch_id'] = "CNT"


# In[ ]:


adata2 = sc.read_10x_mtx(
    os.path.join(data_folder, "toxoHA_raw_feature_bc_matrix"),
    var_names='gene_symbols',
    cache=True)   
adata2.obs['batch_id'] = "HA"


# In[ ]:


adata3 = sc.read_10x_mtx(
    os.path.join(data_folder, "toxoMeCP2_raw_feature_bc_matrix"),
    var_names='gene_symbols',
    cache=True)   
adata3.obs['batch_id'] = "MECP2"


# In[ ]:


print(adata1.shape) 
print(adata2.shape) 
print(adata3.shape) 


# In[ ]:


anndata_list = []
anndata_list.append(adata1)
anndata_list.append(adata2)
anndata_list.append(adata3)
batch_list = ["CNT", "HA", "MECP2"]

adata = anndata_list[0].concatenate(anndata_list[1:],
                                                     join='outer',
                                                     batch_key='batch_id',
                                                     batch_categories=batch_list,
                                                     uns_merge=None,
                                                     index_unique='-',
                                                     fill_value=0.0)
    

adata


# In[ ]:


#Clean up environment
del adata1
del adata2
del adata3


# ## read CMO info and cellranger filtered barcodes

# In[ ]:


cmo_sample_id = pd.read_csv(os.path.join(data_folder, "MultiSeq_annotation_prefiltered_mtx.csv"), index_col=0)
cmo_sample_id


# In[ ]:


adata = adata[cmo_sample_id.index]
adata


# In[ ]:


adata.obs["sample_id"] = cmo_sample_id["Sample"]
adata.obs["sample_id"]


# In[ ]:


adata.obs["sample_id"].unique()


# ## Initial numbers

# In[ ]:


#adata.obs 
print('Initial number of cells:', adata.n_obs) 
 
# To see the row names: 
print('Cell names: ', adata.obs_names[:5].tolist()) 
 
# To see the columns of the metadata (information available for each cell)  
print('Available metadata for each cell: ', adata.obs.columns) 


# In[ ]:


print('Initial number of genes:', adata.n_vars) 
 
# To see the columns names: 
print('Gene names: ', adata.var_names[:5].tolist()) 
 
# To see the gene metadata (information available for each gene)  
print('Available metadata for each gene: ', adata.var.columns) 


# In[ ]:


adata.obs['batch_id'].value_counts()


# In[ ]:


adata.obs['batch_id'].value_counts().plot.bar(color=['orange', 'magenta', 'limegreen'])


# In[ ]:


adata.obs['sample_id'].value_counts()


# In[ ]:


adata.obs['sample_id'].value_counts()/adata.obs['sample_id'].value_counts().sum()*100


# In[ ]:


adata.obs['sample_id'].value_counts().plot.bar()


# # Quality checks

# ## Top-Expressed genes

# In[ ]:


sc.pl.highest_expr_genes(adata, n_top=20)


# ## Automated QC metrics
# 

# **The string based selection is not completely accurate** because it can also include genes that are not mitocondrial/ribosomal but share the string used for the selection. It should be anyway a good enough approximation for our purposes.
# *CAREFUL: str.contains must be used with | for the 2 alternatives, cannot be used with the same synthax as startswith with 2 strings.*
# 

# In[ ]:


#qc_vars wants a column of adata.var containing T/F or 1/0 indicating the genes to be selected for sub-statistics
adata.var['ribo']= adata.var_names.str.contains('^cellranger_gex-GRCh38-2020-A_RPS|^cellranger_gex-GRCh38-2020-A_RPL')
# CAREFUL: str.contains must be used with | for the 2 alternatives, cannot be used with the same synthax as startswith with 2 strings 
#adata.var['rb'] = adata.var_names.str.startswith(("RPS","RPL"))
adata.var['mito']= adata.var_names.str.startswith('cellranger_gex-GRCh38-2020-A_MT-')
adata.var['toxo'] = adata.var_names.str.startswith('ToxoDB_tgondii_ME49_')
adata.var['human'] = adata.var_names.str.startswith('cellranger_gex-GRCh38-2020-A')
adata.var['HA'] = adata.var_names.str.startswith('ToxoDB_tgondii_ME49_mod______HAstop')
adata.var['Mecp2'] = adata.var_names.str.startswith('ToxoDB_tgondii_ME49_mod______MeCP2opt')


# In[ ]:


#Compute metrics (inplace=True to append to adata)
sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo', 'mito', 'toxo','human','HA','Mecp2'], inplace=True,
                           log1p=False, percent_top=None)
#adata.obs


# In[ ]:


adata.obs


# ## Inspect quality-related parameters

# In[ ]:


fig, ax = plt.subplots(figsize=(17,5))
ax.set_ylim(0, 60000)
sc.pl.violin(adata, ['total_counts'], groupby='sample_id', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()


# In[ ]:


fig, ax = plt.subplots(figsize=(17,5))
ax.set_ylim(0, 10000)
sc.pl.violin(adata, ['n_genes_by_counts'], groupby='sample_id', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)


# In[ ]:


fig, ax = plt.subplots(figsize=(17,5))
ax.set_ylim(0, 20)
sc.pl.violin(adata, ['pct_counts_mito'], groupby='sample_id', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)


# In[ ]:


fig, ax = plt.subplots(figsize=(17,5))
ax.set_ylim(0, 30)
sc.pl.violin(adata, ['pct_counts_ribo'], groupby='sample_id', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)


# In[ ]:


fig, ax = plt.subplots(figsize=(17,5))
ax.set_ylim(0, 10)
sc.pl.violin(adata, ['pct_counts_toxo'], groupby='sample_id', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)


# In[ ]:


fig, ax = plt.subplots(figsize=(17,5))
ax.set_ylim(90, 100)
sc.pl.violin(adata, ['pct_counts_human'], groupby='sample_id', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)


# In[ ]:


fig, ax = plt.subplots(figsize=(17,5))
ax.set_ylim(0, 0.3)
sc.pl.violin(adata, ['pct_counts_HA'], groupby='sample_id', multi_panel=True, jitter=False, log=False, rotation=90., scale="count", ax=ax)


# In[ ]:


fig, ax = plt.subplots(figsize=(17,5))
ax.set_ylim(0, 0.3)
sc.pl.violin(adata, ['pct_counts_Mecp2'], groupby='sample_id', multi_panel=True, jitter=False, log=False, rotation=90., scale="count", ax=ax)


# # Filters

# ## Filter "doublets" and "negative" cells

# In[ ]:


adata = adata[~adata.obs["sample_id"].isin(["Doublet", "Negative"])]
adata


# ## Thresholds

# In[ ]:


MIN_CELLS = 200
#
MIN_GENES = 1000
MAX_GENES = 8000
MIN_COUNTS= 2000
MAX_COUNTS= 30000
#
MT_PERCENTAGE = 15
RIBO_PERCENTAGE = 35
HUMAN_PERCENTAGE = 50


# In[ ]:


plt.figure(figsize=(13,5), tight_layout=True)

ax1 = plt.subplot(1, 2, 1)
sns.histplot(adata.obs['total_counts'], stat="count", bins=500, color='chocolate', kde=True, ax=ax1)
plt.axvline(MIN_COUNTS, 0, 1, c='red')
plt.axvline(MAX_COUNTS, 0, 1, c='red')
ax1.set_xlim([0., 35000.])

ax2 = plt.subplot(1, 2, 2)
sns.histplot(adata.obs['n_genes_by_counts'], stat="count", bins=100, color='orange', kde=True, ax=ax2)
plt.axvline(MIN_GENES, 0, 1, c='red')
plt.axvline(MAX_GENES, 0, 1, c='red')
ax2.set_xlim([0., 10000.])

plt.show()


# In[ ]:


plt.figure(figsize=(13,5), tight_layout=True)

ax1 = plt.subplot(1, 2, 1)
sns.histplot(adata.obs['pct_counts_mito'], stat="count", bins=100, kde=True, color='limegreen', ax=ax1)
plt.axvline(MT_PERCENTAGE, 0, 1, c='red')
ax1.set_xlim([0., 20.])

ax2 = plt.subplot(1, 2, 2)
sns.histplot(adata.obs['pct_counts_ribo'], stat="count", bins=100, kde=True, color='deepskyblue', ax=ax2)
plt.axvline(RIBO_PERCENTAGE, 0, 1, c='red')
ax2.set_xlim([0., 60.])

plt.show()


# In[ ]:


plt.figure(figsize=(13,5), tight_layout=True)

ax1 = plt.subplot(1, 2, 1)
sns.histplot(adata.obs['pct_counts_human'], stat="count", bins=100, kde=True, color='limegreen', ax=ax1)
plt.axvline(HUMAN_PERCENTAGE, 0, 1, c='red')
ax1.set_xlim([0., 100.])
ax1.set_ylim([0., 60.])

ax2 = plt.subplot(1, 2, 2)
sns.histplot(adata.obs['pct_counts_toxo'], stat="count", bins=100, kde=True, color='deepskyblue', ax=ax2)
ax2.set_xlim([0., 10.])

plt.show()


# ## Filtering cells

# ### Detected genes

# In[ ]:


sc.pp.filter_cells(adata, min_genes=MIN_GENES)
print('After filtering on min detected genes:number of cells:', adata.n_obs)
print('After filtering on min detected genes:number of genes:', adata.n_vars)
print()


# In[ ]:


sc.pp.filter_cells(adata, max_genes=MAX_GENES)
print('After filtering on max detected genes:number of cells:', adata.n_obs)
print('After filtering on max detected genes:number of genes:', adata.n_vars)


# ### UMI counts

# In[ ]:


sc.pp.filter_cells(adata, min_counts=MIN_COUNTS)
print('After filtering on min UMI counts:number of cells:', adata.n_obs)
print('After filtering on min UMI counts:number of genes:', adata.n_vars)
print()


# In[ ]:


sc.pp.filter_cells(adata, max_counts=MAX_COUNTS)
print('After filtering on min UMI counts:number of cells:', adata.n_obs)
print('After filtering on min UMI counts:number of genes:', adata.n_vars)
print()


# ### Mitochondrial RNA

# In[ ]:


adata = adata[adata.obs['pct_counts_mito'] < MT_PERCENTAGE, :]

print('After filtering on mitochondrial RNA: number of cells:', adata.n_obs)


# ### Ribosomal RNA

# In[ ]:


adata = adata[adata.obs['pct_counts_ribo'] < RIBO_PERCENTAGE, :]

print('After filtering on ribosomal protein RNA: number of cells:', adata.n_obs)


# ###  Human genes percentage 

# In[ ]:


adata = adata[adata.obs['pct_counts_human'] > HUMAN_PERCENTAGE, :]

print('After filtering on only human genes: number of cells:', adata.n_obs)


# ## Filtering genes
# 
# We implement a manual filtering so that beyond the genes respecting the threhsold also our genes of interest from the constructs are kept. 
# 

# In[ ]:


print('Before gene filtering: number of genes:', adata.n_vars)
print('Before gene filtering: number of cells:', adata.n_obs)


# In[ ]:


keep = sc.pp.filter_genes(adata, min_cells=MIN_CELLS, inplace=False)

keep[0][-1] = True
keep[0][-2] = True
keep[0][-3] = True

adata = adata[:, keep[0]]


# In[ ]:


adata.var_names


# In[ ]:


print('After gene filtering: number of genes:', adata.n_vars)
print('After filtering: number of cells:', adata.n_obs)


# # Numbers and Visualization after filtering

# In[ ]:


print('After applied filtering: number of cells:', adata.n_obs)
print('After applied filtering: number of genes:', adata.n_vars)


# ## Violin plots
# 

# In[ ]:


fig, ax = plt.subplots(figsize=(7,5), tight_layout=True)

sc.pl.violin(adata, keys=['n_genes_by_counts'],
              multi_panel=True, groupby='sample_id', rotation=45, jitter=False, log=True, ax=ax)


# In[ ]:


fig, ax = plt.subplots(figsize=(7,5), tight_layout=True)

sc.pl.violin(adata, keys=['total_counts'],
              multi_panel=True, groupby='sample_id', rotation=45, jitter=False, log=True, ax=ax)


# In[ ]:


fig, ax = plt.subplots(figsize=(7,5), tight_layout=True)

sc.pl.violin(adata, keys=['pct_counts_mito'],
              multi_panel=True, groupby='sample_id', rotation=45, jitter=False, log=False, ax=ax)


# In[ ]:


fig, ax = plt.subplots(figsize=(7,5), tight_layout=True)

sc.pl.violin(adata, keys=['pct_counts_ribo'],
              multi_panel=True, groupby='sample_id', rotation=45, jitter=False, log=False, ax=ax)


# In[ ]:


fig, ax = plt.subplots(figsize=(7,5), tight_layout=True)

sc.pl.violin(adata, keys=['pct_counts_toxo'],
              multi_panel=True, groupby='sample_id', rotation=45, jitter=False, log=False, ax=ax)


# In[ ]:


fig, ax = plt.subplots(figsize=(7,5), tight_layout=True)

sc.pl.violin(adata, keys=['pct_counts_human'],
              multi_panel=True, groupby='sample_id', rotation=45, jitter=False, log=False, ax=ax)


# ## Density plots

# In[ ]:


#Plot them in line so they take up less space
plt.figure(figsize=(20,5))

plt.subplot(1, 5, 1)
sns.kdeplot(np.log10(adata.obs['n_genes_by_counts']), shade=True, color='cornflowerblue')
plt.axvline(np.log10(MIN_GENES), 0, 1, c='red')  #set manually for chosen threshold

plt.subplot(1, 5, 2)
sns.kdeplot(np.log10(adata.obs['total_counts']), shade=True, color='forestgreen')
plt.axvline(np.log10(MIN_COUNTS), 0, 1, c='red')  #set manually for chosen threshold

plt.subplot(1, 5, 3)
sns.kdeplot(adata.obs['pct_counts_mito'], shade=True, color='coral')
plt.axvline(MT_PERCENTAGE, 0, 1, c='red')  #set manually for chosen threshold

plt.subplot(1, 5, 4)
sns.kdeplot(adata.obs['pct_counts_ribo'], shade=True, color='orchid')
plt.axvline(RIBO_PERCENTAGE, 0, 1, c='red')  #set manually for chosen threshold


# ## batches and IDs

# In[ ]:


adata.obs['batch_id'].value_counts().plot.bar(color=['orange', 'magenta', 'limegreen'])


# In[ ]:


adata.obs['sample_id'].value_counts().plot.bar()


# # Dimensionality reduction and annotation of the full dataset

# ## Normalize and Log-transform
# For this quick exploration I adopt the basic normalization from scanpy

# In[ ]:


adata.layers["counts"] = adata.X.copy()


# In[ ]:


sc.pp.normalize_total(adata, exclude_highly_expressed=True)
sc.pp.log1p(adata)


# ## Triku gene selection

# In[ ]:


sc.pp.pca(adata)
sc.pp.neighbors(adata, metric='cosine', n_neighbors=int(0.5 * len(adata) ** 0.5))


# In[ ]:


tk.tl.triku(adata, use_raw=False)


# In[ ]:


Top20Triku = adata.var.sort_values(by=['triku_distance'], ascending=False).head(20).index
Top20Triku


# In[ ]:


print('Number of Higly Variable Genes', len(adata.var_names[adata.var['highly_variable'] == True])) ##numero HVG


# ## PCA

# In[ ]:


sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')


# In[ ]:


sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)


# In[ ]:


sc.pl.pca(adata, color=['n_genes_by_counts', 'total_counts', 'pct_counts_mito', 'pct_counts_ribo', 'pct_counts_toxo', 'pct_counts_human'])


# In[ ]:


sc.pl.pca(adata, color=['batch_id', 'sample_id'])


# ## UMAP

# In[ ]:


sc.pp.neighbors(adata, n_neighbors=80, n_pcs=50)


# In[ ]:


sc.tl.umap(adata, random_state=1)


# In[ ]:


sc.pl.umap(adata, color=['n_genes_by_counts', 'total_counts', 'pct_counts_mito', 'pct_counts_ribo', 'pct_counts_toxo', 'pct_counts_human'], vmin='p1', vmax='p99')


# In[ ]:


sc.pl.umap(adata, color=['batch_id', 'sample_id'], vmin='p1', vmax='p99')


# In[ ]:


sc.pl.umap(adata, color=['ToxoDB_tgondii_ME49_mod______GRA16',
                         'ToxoDB_tgondii_ME49_mod______HAstop',
                         'ToxoDB_tgondii_ME49_mod______MeCP2opt'],
           size=50, cmap="Reds", add_outline=True, vmax=1)


# In[ ]:


sc.tl.embedding_density(adata, groupby='batch_id')
sc.pl.embedding_density(adata, groupby='batch_id')


# In[ ]:


Top10Triku = adata.var.sort_values(by=['triku_distance'], ascending=False).head(10).index
sc.pl.umap(adata, color=Top10Triku, vmin='p1', vmax='p99')


# ## Cluster identification

# In[ ]:


res = [0.3, 0.4, 0.5, 0.6]
leiden_labels = []

for x in res:
    label = "Leiden_" + str(x).replace('.', '')
    leiden_labels.append(label) 
    sc.tl.leiden(adata, resolution = x, key_added= label)


# In[ ]:


sc.pl.umap(adata, color=leiden_labels)


# ## *bona-fide* infected cells

# In[ ]:


adata.obs["infected"] = (adata.obs['Leiden_03'] == "6").astype(str)


# In[ ]:


adata.obs["batch_infect"] = adata.obs["batch_id"].astype(str) + "_" + adata.obs["infected"].map({"True": "infect", "False": "not_infect"}).astype(str)


# In[ ]:


sc.pl.umap(adata, color="infected")


# In[ ]:


cell_numbers = pd.crosstab(adata.obs["batch_id"], adata.obs["infected"])
cell_numbers


# In[ ]:


cell_frac = pd.crosstab(adata.obs["batch_id"], adata.obs["infected"], normalize="index")
cell_frac


# In[ ]:


fig1, [ax1, ax2] = plt.subplots(ncols=1, nrows=2, figsize=(5, 8), sharex=True)
cell_numbers.plot.bar(stacked=True, ax=ax1).legend(loc='lower right')
cell_frac.plot.bar(stacked=True, ax=ax2).legend(loc='lower right')
fig1.tight_layout()
fig1.show()


# In[ ]:


cell_numbers = pd.crosstab(adata.obs["sample_id"], adata.obs["infected"])
cell_numbers


# In[ ]:


cell_frac = pd.crosstab(adata.obs["sample_id"], adata.obs["infected"], normalize="index")
cell_frac


# In[ ]:


fig1, [ax1, ax2] = plt.subplots(ncols=1, nrows=2, figsize=(5, 8), sharex=True)
cell_numbers.plot.bar(stacked=True, ax=ax1).legend(loc='lower right')
cell_frac.plot.bar(stacked=True, ax=ax2).legend(loc='lower right')
fig1.tight_layout()
fig1.show()


# ## Diffusion map

# In[ ]:


sc.tl.diffmap(adata)


# In[ ]:


sc.pl.diffmap(adata, color=['n_genes_by_counts', 'total_counts', 'pct_counts_ribo', 'pct_counts_mito'])


# In[ ]:


sc.pl.diffmap(adata, color=['batch_id', 'sample_id', 'infected'], size=2)


# ## Draw graph

# In[ ]:


sc.tl.draw_graph(adata)


# In[ ]:


sc.pl.draw_graph(adata, color=['n_genes_by_counts', 'total_counts', 'pct_counts_ribo', 'pct_counts_mito'])


# In[ ]:


sc.pl.draw_graph(adata, color=['batch_id', 'sample_id', 'infected'], size=2)


# ## Marker plots

# **the marker plots are collected in the notebook "02-markers.ipynb"**

# # Dimensionality reduction and annotation of the human-only dataset

# ## human genes only - recover raw counts

# In[ ]:


onlyhuman = sc.AnnData(
    obs=adata[:, adata.var_names.str.startswith('cellranger_gex-GRCh38-2020-A')].obs.copy(),
    var=adata[:, adata.var_names.str.startswith('cellranger_gex-GRCh38-2020-A')].var.copy(),
    X=adata[:, adata.var_names.str.startswith('cellranger_gex-GRCh38-2020-A')].layers["counts"]
)


# In[ ]:


onlyhuman.var.index = onlyhuman.var.index.str.replace("cellranger_gex-GRCh38-2020-A_", "")


# In[ ]:


onlyhuman


# ## Normalize and Log-transform
# For this quick exploration I adopt the basic normalization from scanpy

# In[ ]:


sc.pp.normalize_total(onlyhuman, exclude_highly_expressed=True)
sc.pp.log1p(onlyhuman)


# ## Triku gene selection

# In[ ]:


sc.pp.pca(onlyhuman)
sc.pp.neighbors(onlyhuman, metric='cosine', n_neighbors=int(0.5 * len(adata) ** 0.5))


# In[ ]:


tk.tl.triku(onlyhuman, use_raw=False)


# In[ ]:


Top20Triku = onlyhuman.var.sort_values(by=['triku_distance'], ascending=False).head(20).index
Top20Triku


# In[ ]:


print('Number of Higly Variable Genes', len(onlyhuman.var_names[onlyhuman.var['highly_variable'] == True])) ##numero HVG


# ## PCA

# In[ ]:


sc.pp.pca(onlyhuman, n_comps=50, use_highly_variable=True, svd_solver='arpack')


# In[ ]:


sc.pl.pca_variance_ratio(onlyhuman, log=True, n_pcs=50)


# In[ ]:


sc.pl.pca(onlyhuman, color=['n_genes_by_counts', 'total_counts', 'pct_counts_toxo', 'pct_counts_human'])


# In[ ]:


sc.pl.pca(onlyhuman, color=['batch_id', 'sample_id', 'infected'])


# ## UMAP

# In[ ]:


sc.pp.neighbors(onlyhuman, n_neighbors=80, n_pcs=18)


# In[ ]:


sc.tl.umap(onlyhuman, random_state=1)


# In[ ]:


sc.pl.umap(onlyhuman, color=['n_genes_by_counts', 'total_counts', 'pct_counts_mito', 'pct_counts_ribo', 'pct_counts_toxo', 'pct_counts_human'], vmin='p1', vmax='p99')


# In[ ]:


sc.pl.umap(onlyhuman, color=['batch_id', 'sample_id',], vmin='p1', vmax='p99')


# In[ ]:


sc.pl.umap(onlyhuman, color="infected", groups="True", size=50, add_outline=False,)


# In[ ]:


sc.pl.umap(onlyhuman[onlyhuman.obs["batch_id"] == "CNT"], color=['batch_id', 'sample_id'], vmin='p1', vmax='p99')
sc.pl.umap(onlyhuman[onlyhuman.obs["batch_id"] == "HA"], color=['batch_id', 'sample_id'], vmin='p1', vmax='p99')
sc.pl.umap(onlyhuman[onlyhuman.obs["batch_id"] == "MECP2"], color=['batch_id', 'sample_id'], vmin='p1', vmax='p99')


# In[ ]:


sc.tl.embedding_density(onlyhuman, groupby='batch_id')
sc.pl.embedding_density(onlyhuman, groupby='batch_id')


# ## Diffusion maps

# In[ ]:


sc.tl.diffmap(onlyhuman)


# In[ ]:


sc.pl.diffmap(onlyhuman, color=['batch_id', 'pct_counts_ribo', 'sample_id'], size=2)


# In[ ]:


sc.pl.diffmap(onlyhuman[onlyhuman.obs["batch_id"] == "CNT"], color=['batch_id', 'sample_id'], vmin='p1', vmax='p99')
sc.pl.diffmap(onlyhuman[onlyhuman.obs["batch_id"] == "HA"], color=['batch_id', 'sample_id'], vmin='p1', vmax='p99')
sc.pl.diffmap(onlyhuman[onlyhuman.obs["batch_id"] == "MECP2"], color=['batch_id', 'sample_id'], vmin='p1', vmax='p99')


# ## Draw graph

# In[ ]:


sc.tl.draw_graph(onlyhuman)


# In[ ]:


sc.pl.draw_graph(onlyhuman, color=['n_genes_by_counts', 'total_counts', 'pct_counts_ribo', 'pct_counts_mito'])


# In[ ]:


sc.pl.draw_graph(onlyhuman, color=['batch_id', 'sample_id', 'infected'], size=2)


# ## Marker plots

# **the marker plots are collected in the notebook "02-markers.ipynb"**

# ## Cluster identification

# In[ ]:


res = [0.3, 0.4, 0.5, 0.6]
leiden_labels = []

for x in res:
    label = "Leiden_" + str(x).replace('.', '')
    leiden_labels.append(label) 
    sc.tl.leiden(onlyhuman, resolution = x, key_added= label)


# In[ ]:


sc.pl.umap(onlyhuman, color=leiden_labels)


# In[ ]:


sc.pl.diffmap(onlyhuman, color=leiden_labels)


# In[ ]:


chosen_leiden = 'Leiden_06'
key_leiden = 'rank_L' + chosen_leiden[-2:]


# In[ ]:


key_leiden


# In[ ]:


sc.tl.rank_genes_groups(onlyhuman, chosen_leiden, method='wilcoxon', key_added=key_leiden)


# In[ ]:


GroupMarkers = pd.DataFrame(onlyhuman.uns[key_leiden]['names']).head(41)
GroupMarkers.columns = 'Cl_' + GroupMarkers.columns
GroupMarkers.head(26)


# ## Cluster functional analysis

# In[ ]:


gp = GProfiler(return_dataframe=True)


# **Cluster 0**

# In[ ]:


CustomGO(onlyhuman, cluster='0', rank=key_leiden, n_markers=50,  show=8)


# **Cluster 1**

# In[ ]:


CustomGO(onlyhuman, cluster='1', rank=key_leiden, n_markers=50,  show=8)


# **Cluster 2**

# In[ ]:


CustomGO(onlyhuman, cluster='2', rank=key_leiden, n_markers=50,  show=8)


# **Cluster 3**

# In[ ]:


CustomGO(onlyhuman, cluster='3', rank=key_leiden, n_markers=50,  show=8)


# **Cluster 4**

# In[ ]:


CustomGO(onlyhuman, cluster='4', rank=key_leiden, n_markers=50,  show=8)


# **Cluster 5**

# In[ ]:


CustomGO(onlyhuman, cluster='5', rank=key_leiden, n_markers=50,  show=8)


# **Cluster 6**

# In[ ]:


CustomGO(onlyhuman, cluster='6', rank=key_leiden, n_markers=50,  show=8)


# **Cluster 7**

# In[ ]:


CustomGO(onlyhuman, cluster='7', rank=key_leiden, n_markers=50,  show=8)


# **Cluster 8**

# In[ ]:


CustomGO(onlyhuman, cluster='8', rank=key_leiden, n_markers=50,  show=8)


# **Cluster 9**

# In[ ]:


CustomGO(onlyhuman, cluster='9', rank=key_leiden, n_markers=50,  show=8)


# **Cluster 10**

# In[ ]:


CustomGO(onlyhuman, cluster='10', rank=key_leiden, n_markers=50,  show=8)


# **Cluster 11**

# In[ ]:


CustomGO(onlyhuman, cluster='11', rank=key_leiden, n_markers=50,  show=8)


# ## Cluster annotation

# In[ ]:


sc.pl.umap(onlyhuman, color=[chosen_leiden], legend_loc='on data')


# In[ ]:


onlyhuman.obs["celltype"] = onlyhuman.obs[chosen_leiden].map(
 {"0": "N_IP",
  "1": "N_UPR",
  "2": "N_metabolism",
  "3": "N1",
  "4": "N2",
  "5": "vRG_oRG",
  "6": "N3",
  "7": "N_UPR2",
  "8": "N_UPR3", 
  "9": "vRG_oRG2",
  "10": "CyclingProg",
  "11": "N_Proj",}   
)


# In[ ]:


sc.pl.umap(onlyhuman, color=["celltype"], legend_loc='on data')


# In[ ]:


adata.obs["celltype"] = onlyhuman.obs["celltype"]


# In[ ]:


fig, ax = plt.subplots(figsize=(8,3))
sc.pl.violin(onlyhuman, ['pct_counts_ribo'], groupby='celltype', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(8,3))
sc.pl.violin(onlyhuman, ['pct_counts_toxo'], groupby='celltype', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(8,3))
sc.pl.violin(onlyhuman, ['total_counts'], groupby='celltype', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()


# # Metrics of infected vs infected neurons

# In[ ]:


onlyhuman.layers['scaled'] = sc.pp.scale(onlyhuman, copy=True).X


# In[ ]:


onlyhuman.obs["sample_infect"] = onlyhuman.obs["sample_id"].astype(str) + "_" +  onlyhuman.obs["infected"].astype(str)


# In[ ]:


N_clusters = onlyhuman[onlyhuman.obs["celltype"].isin(['N1', 'N2', 'N3'])].copy()
N_clusters


# In[ ]:


fig, ax = plt.subplots(figsize=(4,3))
sc.pl.violin(N_clusters, ['n_genes_by_counts'], groupby='batch_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(4,3))
sc.pl.violin(N_clusters, ['total_counts'], groupby='batch_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(4,3))
sc.pl.violin(N_clusters, ['total_counts_ribo'], groupby='batch_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(4,3))
sc.pl.violin(N_clusters, ['pct_counts_ribo'], groupby='batch_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(4,3))
sc.pl.violin(N_clusters, ['total_counts_mito'], groupby='batch_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(4,3))
sc.pl.violin(N_clusters, ['pct_counts_mito'], groupby='batch_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(4,3))
sc.pl.violin(N_clusters, ['pct_counts_toxo'], groupby='batch_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(4,3))
sc.pl.violin(N_clusters, ['pct_counts_human'], groupby='batch_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()


# In[ ]:


fig, ax = plt.subplots(figsize=(6,3))
sc.pl.violin(N_clusters, ['n_genes_by_counts'], groupby='sample_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(6,3))
sc.pl.violin(N_clusters, ['total_counts'], groupby='sample_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(6,3))
sc.pl.violin(N_clusters, ['total_counts_ribo'], groupby='sample_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(6,3))
sc.pl.violin(N_clusters, ['pct_counts_ribo'], groupby='sample_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(6,3))
sc.pl.violin(N_clusters, ['total_counts_mito'], groupby='sample_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(6,3))
sc.pl.violin(N_clusters, ['pct_counts_mito'], groupby='sample_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(6,3))
sc.pl.violin(N_clusters, ['pct_counts_toxo'], groupby='sample_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(6,3))
sc.pl.violin(N_clusters, ['pct_counts_human'], groupby='sample_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()


# # Saving

# ## Save adata

# In[ ]:


results_file = os.path.join(data_folder, "annotated_dataset_human.h5ad")
onlyhuman.write(results_file)


# In[ ]:


results_file = os.path.join(data_folder, "annotated_dataset_both.h5ad")
adata.write(results_file)


# ## Timestamp finished computations

# In[ ]:


print(datetime.now())


# ## Save python and html version of notebook

# In[ ]:


nb_fname = ipynbname.name()
nb_fname


# In[ ]:


get_ipython().run_cell_magic('bash', '-s "$nb_fname"', 'sleep 120\njupyter nbconvert "$1".ipynb --to="python" --ClearOutputPreprocessor.enabled=True\njupyter nbconvert "$1".ipynb --to="html"')


# In[ ]:




