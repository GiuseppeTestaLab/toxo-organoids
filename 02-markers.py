#!/usr/bin/env python
# coding: utf-8

# <font size="8">Marker plots</font>

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
sc.settings.set_figure_params(dpi=60)
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


data_folder = "./data"
result_folder = './results'


# In[ ]:


print(datetime.now())


# # Human and toxo dataset

# In[ ]:


def CustomUmap(adata, genes):
    genes = ['cellranger_gex-GRCh38-2020-A_' + gene for gene in genes]
    genes = selectMarkers(adata, genes)
    sc.pl.umap(adata, color=genes, size=10, frameon=False,
               vmin='p1',  vmax='p99')


# In[ ]:


def CustomDA(adata, genes):
    genes = ['cellranger_gex-GRCh38-2020-A_' + gene for gene in genes]
    genes = selectMarkers(adata, genes)
    sc.pl.draw_graph(adata, color=genes, size=10, frameon=False,
               vmin='p1',  vmax='p99')


# ## Data Load

# In[ ]:


adata = sc.read(os.path.join(data_folder, "annotated_dataset_both.h5ad"))


# In[ ]:


adata.shape 


# In[ ]:


adata.obs['batch_id'].value_counts()


# In[ ]:


#adata.obs 
print('Number of cells:', adata.n_obs) 
print('Number of genes:', adata.n_vars) 


# ## UMAP: marker genes 

# In[ ]:


sc.pl.umap(adata, color=['batch_id'], vmin='p1', vmax='p99', frameon=False)


# ### Top-20 Triku

# In[ ]:


Top10Triku = adata.var.sort_values(by=['triku_distance'], ascending=False).head(20).index
sc.pl.umap(adata, color=Top10Triku, vmin='p1', vmax='p99', frameon=False)
#Top20Triku


# ### Proliferation

# In[ ]:


CustomUmap(adata, genes=['MKI67', 'CDC20', 'HMGB1', 'HMGB2', 'CCNB1', 'CCNB2'])


# ### Apical Progenitors

# In[ ]:


CustomUmap(adata, genes=['FABP7', 'GLI3', 'PAX6', 'NES', 'SOX1', 'SOX2', 'SOX9', 'VIM'])


# ### Intermediate Progenitors

# In[ ]:


CustomUmap(adata, genes=['ELAVL4', 'NHLH1', 'KCNQ3'])


# ### Pan-neuronal

# In[ ]:


CustomUmap(adata, genes=['DCX', 'STMN2', 'SYT1', 'SYP', 'GAP43', 'MEF2C',  'FOXG1'])


# ### Excitatory

# In[ ]:


CustomUmap(adata, genes=['NEUROD1', 'NEUROD2', 'NEUROD6', 'NEUROG2', #excitatory progenitors  
                        'SLC17A7', 'TBR1', 'BCL11A', 'CUX1', 'CUX2',
                         'SATB2', 'POU3F2', 'POU3F3'])


# ### Inhibitory

# In[ ]:


CustomUmap(adata, genes=['DLX5', 'DLX6', 'DLX1', 'DLX2', 'DLX6-AS1', 'CALB2', 'GAD1', 'GAD2', 'PROX1', 'NR2F2', 'EGFR', 'NKX2-1', 
       'SST', 'NPY', 'CALB1', 'CALB2', 'VIP', 'PVALB'])


# ### oRG

# In[ ]:


CustomUmap(adata, genes=['FAM107A', 'HOPX',  'PTPRZ1', 'TNC', 'ITGB5'])


# ### Astrocytes

# In[ ]:


CustomUmap(adata, genes=['GFAP', 'SLC1A3',  'S100B', 'AQP4', 'FGFR3', 'ALDH1L1'])


# ### Off-target markers

# In[ ]:


CustomUmap(adata, genes=['HOXA3', 'HOXB3', 'VSX2', 'OTX2', 'RSPO2', 'TTR', 'MYH3', 'MYH8', 'MYL1', 'MYLPF'] )


# ### Stress-related markers

# In[ ]:


CustomUmap(adata, genes=['ALDOA', 'BNIP3', 'PGK1', 'ARCN1', 'GORASP2'] )


# ### Interferon-sensitive genes

# In[ ]:


CustomUmap(adata, genes=['MX1', 'IFIT1', 'ISG15', 'STAT1', 'IRF7', 'IRF3', 'IRF1', 'IRF9'])


# ### Cytokines

# In[ ]:


CustomUmap(adata, genes=['GFAP', 'CXCL1', 'CCL3', 'CCL4', 'TGFB1', 'IL1A', 'IL1B', 'IL6', 'CXCL8',
        'TNF', 'IFNG', 'IL27', 'MIF', 'SERPINE1'])


# ## Force Graph: marker genes

# In[ ]:


sc.pl.draw_graph(adata, color=['batch_id'], frameon=False)


# ### Top-20 Triku

# In[ ]:


sc.pl.draw_graph(adata, color=Top10Triku, vmin='p1', vmax='p99', frameon=False)


# ### Proliferation

# In[ ]:


CustomDA(adata, genes=['MKI67', 'CDC20', 'HMGB1', 'HMGB2', 'CCNB1', 'CCNB2'])


# ### Apical Progenitors

# In[ ]:


CustomDA(adata, genes=['FABP7', 'GLI3', 'PAX6', 'NES', 'SOX1', 'SOX2', 'SOX9', 'VIM'])


# ### Intermediate Progenitors

# In[ ]:


CustomDA(adata, genes=['EOMES', 'ELAVL4', 'NHLH1', 'KCNQ3'])


# ### Pan-neuronal

# In[ ]:


CustomDA(adata, genes=['DCX', 'STMN2', 'SYT1', 'SYP', 'GAP43', 'MEF2C',  'FOXG1'])


# ### Excitatory

# In[ ]:


CustomDA(adata, genes=['NEUROD1', 'NEUROD2', 'NEUROD6', 'NEUROG2', #excitatory progenitors  
                        'SLC17A7', 'TBR1', 'BCL11A', 'CUX1', 'CUX2',
                         'SATB2', 'POU3F2', 'POU3F3'])


# ### Inhibitory

# In[ ]:


CustomDA(adata, genes=['DLX5', 'DLX6', 'DLX1', 'DLX2', 'DLX6-AS1', 'CALB2', 'GAD1', 'GAD2', 'PROX1', 'NR2F2', 'EGFR', 'NKX2-1', 
       'SST', 'NPY', 'CALB1', 'CALB2', 'VIP', 'PVALB'])


# ### oRG

# In[ ]:


CustomDA(adata, genes=['FAM107A', 'HOPX',  'PTPRZ1', 'TNC', 'ITGB5'])


# ### Astrocytes

# In[ ]:


CustomDA(adata, genes=['GFAP', 'SLC1A3',  'S100B', 'AQP4', 'FGFR3', 'ALDH1L1'])


# ### Off-target markers

# In[ ]:


CustomDA(adata, genes=['HOXA3', 'HOXB3', 'VSX2', 'OTX2', 'RSPO2', 'TTR', 'MYH3', 'MYH8', 'MYL1', 'MYLPF'])


# ### Stress markers

# In[ ]:


CustomDA(adata, genes=['ALDOA', 'BNIP3', 'PGK1', 'ARCN1', 'GORASP2'] )


# ### Interferon-sensitive genes

# In[ ]:


CustomDA(adata, genes=['MX1', 'IFIT1', 'ISG15', 'STAT1', 'IRF7', 'IRF3', 'IRF1', 'IRF9'])


# ### Cytokines

# In[ ]:


CustomDA(adata, genes=['GFAP', 'CXCL1', 'CCL3', 'CCL4', 'TGFB1', 'IL1A', 'IL1B', 'IL6', 'CXCL8',
        'TNF', 'IFNG', 'IL27', 'MIF', 'SERPINE1'])


# # "Only-human" dataset

# In[ ]:


def CustomUmap(adata, genes):
    genes = selectMarkers(adata, genes)
    sc.pl.umap(adata, color=genes, size=10, frameon=False,
               vmin='p1',  vmax='p99')


# In[ ]:


def CustomDA(adata, genes):
    genes = selectMarkers(adata, genes)
    sc.pl.draw_graph(adata, color=genes, size=10, frameon=False,
               vmin='p1',  vmax='p99')


# ## Data Load

# In[ ]:


adata = sc.read(os.path.join(data_folder, "annotated_dataset_human.h5ad"))


# In[ ]:


adata.shape 


# In[ ]:


adata.obs['batch_id'].value_counts()


# In[ ]:


#adata.obs 
print('Number of cells:', adata.n_obs) 
print('Number of genes:', adata.n_vars) 


# ## UMAP: marker genes 

# In[ ]:


sc.pl.umap(adata, color=['batch_id'], vmin='p1', vmax='p99', frameon=False)


# ### Top-20 Triku

# In[ ]:


Top10Triku = adata.var.sort_values(by=['triku_distance'], ascending=False).head(20).index
sc.pl.umap(adata, color=Top10Triku, vmin='p1', vmax='p99', frameon=False)
#Top20Triku


# ### Proliferation

# In[ ]:


CustomUmap(adata, genes=['MKI67', 'CDC20', 'HMGB1', 'HMGB2', 'CCNB1', 'CCNB2'])


# ### Apical Progenitors

# In[ ]:


CustomUmap(adata, genes=['FABP7', 'GLI3', 'PAX6', 'NES', 'SOX1', 'SOX2', 'SOX9', 'VIM'])


# ### Intermediate Progenitors

# In[ ]:


CustomUmap(adata, genes=['EOMES', 'ELAVL4', 'NHLH1', 'KCNQ3'])


# ### Pan-neuronal

# In[ ]:


CustomUmap(adata, genes=['DCX', 'STMN2', 'SYT1', 'SYP', 'GAP43', 'MEF2C',  'FOXG1'])


# ### Excitatory

# In[ ]:


CustomUmap(adata, genes=['NEUROD1', 'NEUROD2', 'NEUROD6', 'NEUROG2', #excitatory progenitors  
                        'SLC17A7', 'TBR1', 'BCL11A', 'CUX1', 'CUX2',
                         'SATB2', 'POU3F2', 'POU3F3'])


# ### Inhibitory

# In[ ]:


CustomUmap(adata, genes=['DLX5', 'DLX6', 'DLX1', 'DLX2', 'DLX6-AS1', 'CALB2', 'GAD1', 'GAD2', 'PROX1', 'NR2F2', 'EGFR', 'NKX2-1', 
       'SST', 'NPY', 'CALB1', 'CALB2', 'VIP', 'PVALB'])


# ### oRG

# In[ ]:


CustomUmap(adata, genes=['FAM107A', 'HOPX',  'PTPRZ1', 'TNC', 'ITGB5'])


# ### Astrocytes

# In[ ]:


CustomUmap(adata, genes=['GFAP', 'SLC1A3',  'S100B', 'AQP4', 'FGFR3', 'ALDH1L1'])


# ### Off-target markers

# In[ ]:


CustomUmap(adata, genes=['HOXA3', 'HOXB3', 'VSX2', 'OTX2', 'RSPO2', 'TTR', 'MYH3', 'MYH8', 'MYL1', 'MYLPF'] )


# ### Stress-related markers

# In[ ]:


CustomUmap(adata, genes=['ALDOA', 'BNIP3', 'PGK1', 'ARCN1', 'GORASP2'] )


# ### Interferon-sensitive genes

# In[ ]:


CustomUmap(adata, genes=['MX1', 'IFIT1', 'ISG15', 'STAT1', 'IRF7', 'IRF3', 'IRF1', 'IRF9'])


# ### Cytokines

# In[ ]:


CustomUmap(adata, genes=['GFAP', 'CXCL1', 'CCL3', 'CCL4', 'TGFB1', 'IL1A', 'IL1B', 'IL6', 'CXCL8',
        'TNF', 'IFNG', 'IL27', 'MIF', 'SERPINE1'])


# ## Force Graph: marker genes

# In[ ]:


sc.pl.draw_graph(adata, color=['batch_id'], frameon=False)


# ### Top-20 Triku

# In[ ]:


sc.pl.draw_graph(adata, color=Top10Triku, vmin='p1', vmax='p99', frameon=False)


# ### Proliferation

# In[ ]:


CustomDA(adata, genes=['MKI67', 'CDC20', 'HMGB1', 'HMGB2', 'CCNB1', 'CCNB2'])


# ### Apical Progenitors

# In[ ]:


CustomDA(adata, genes=['FABP7', 'GLI3', 'PAX6', 'NES', 'SOX1', 'SOX2', 'SOX9', 'VIM'])


# ### Intermediate Progenitors

# In[ ]:


CustomDA(adata, genes=['EOMES', 'ELAVL4', 'NHLH1', 'KCNQ3'])


# ### Pan-neuronal

# In[ ]:


CustomDA(adata, genes=['DCX', 'STMN2', 'SYT1', 'SYP', 'GAP43', 'MEF2C',  'FOXG1'])


# ### Excitatory

# In[ ]:


CustomDA(adata, genes=['NEUROD1', 'NEUROD2', 'NEUROD6', 'NEUROG2', #excitatory progenitors  
                        'SLC17A7', 'TBR1', 'BCL11A', 'CUX1', 'CUX2',
                         'SATB2', 'POU3F2', 'POU3F3'])


# ### Inhibitory

# In[ ]:


CustomDA(adata, genes=['DLX5', 'DLX6', 'DLX1', 'DLX2', 'DLX6-AS1', 'CALB2', 'GAD1', 'GAD2', 'PROX1', 'NR2F2', 'EGFR', 'NKX2-1', 
       'SST', 'NPY', 'CALB1', 'CALB2', 'VIP', 'PVALB'])


# ### oRG

# In[ ]:


CustomDA(adata, genes=['FAM107A', 'HOPX',  'PTPRZ1', 'TNC', 'ITGB5'])


# ### Astrocytes

# In[ ]:


CustomDA(adata, genes=['GFAP', 'SLC1A3',  'S100B', 'AQP4', 'FGFR3', 'ALDH1L1'])


# ### Off-target markers

# In[ ]:


CustomDA(adata, genes=['HOXA3', 'HOXB3', 'VSX2', 'OTX2', 'RSPO2', 'TTR', 'MYH3', 'MYH8', 'MYL1', 'MYLPF'])


# ### Stress markers

# In[ ]:


CustomDA(adata, genes=['ALDOA', 'BNIP3', 'PGK1', 'ARCN1', 'GORASP2'] )


# ### Interferon-sensitive genes

# In[ ]:


CustomDA(adata, genes=['MX1', 'IFIT1', 'ISG15', 'STAT1', 'IRF7', 'IRF3', 'IRF1', 'IRF9'])


# ### Cytokines

# In[ ]:


CustomDA(adata, genes=['GFAP', 'CXCL1', 'CCL3', 'CCL4', 'TGFB1', 'IL1A', 'IL1B', 'IL6', 'CXCL8',
        'TNF', 'IFNG', 'IL27', 'MIF', 'SERPINE1'])


# # Saving

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




