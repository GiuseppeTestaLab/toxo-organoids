#!/usr/bin/env python
# coding: utf-8

# <font size="8">Expression levels of MECP2-related genes</font>

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
sc.settings.set_figure_params(dpi=120)
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


def CustomUmap(adata, genes):
    genes = selectMarkers(adata, genes)
    sc.pl.umap(adata, color=genes, size=10, frameon=False,
               vmin='p1',  vmax='p99')


# In[ ]:


def CustomDA(adata, genes):
    genes = selectMarkers(adata, genes)
    sc.pl.draw_graph(adata, color=genes, size=10, frameon=False,
               vmin='p1',  vmax='p99')


# In[ ]:


data_folder = "./data"
result_folder = './results'


# In[ ]:


print(datetime.now())


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


# In[ ]:


N_clusters = adata[adata.obs["celltype"].isin(['N1', 'N2', 'N3'])].copy()
N_clusters


# # Supervised exploration of MECP2 targets in neurons
# 

# ## repressed by MECP2
# 
# https://www.sciencedirect.com/science/article/pii/S0022283619305959
# 
# https://ars.els-cdn.com/content/image/1-s2.0-S0022283619305959-gr3_lrg.jpg
# 
# or upregulated in Rett syndrome

# In[ ]:


MECP2_markers = ['TBL1X', 'SIN3A', 'SKI', 'YY1', 'GTF2B', 'RBPJ', 'PRMT6', 'SP3', 'SOX2', 'SMARCA2',"BDNF", "FKBP5", "IGF2", "DLX5", "DLX6", "SGK1", "MPP1", "GAMT", "FXYD1"]

#for ISG15: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3587786/


#Select present markers for plotting
present_MECP2 = N_clusters.var_names[N_clusters.var_names.isin(MECP2_markers) == True]
present_MECP2 = [x for x in MECP2_markers if x in present_MECP2] #Keep the order of the markers, better way??

#Print missing markers
print('\nThe following marker genes are missing: ', set(MECP2_markers).difference(set(N_clusters.var_names)))


# In[ ]:


sc.pl.umap(N_clusters, color=present_MECP2, vmin='p1', vmax='p99')


# In[ ]:


sc.pl.violin(N_clusters, keys=present_MECP2[0:3],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(N_clusters, keys=present_MECP2[3:6],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(N_clusters, keys=present_MECP2[6:9],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(N_clusters, keys=present_MECP2[9:12],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(N_clusters, keys=present_MECP2[12:15],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(N_clusters, keys=present_MECP2[15:18],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# ## activated by MECP2
# 
# https://www.sciencedirect.com/science/article/pii/S0022283619305959
# 
# https://ars.els-cdn.com/content/image/1-s2.0-S0022283619305959-gr3_lrg.jpg
# 
# or downregulated in Rett Syndrome
# 
# https://ojrd.biomedcentral.com/articles/10.1186/s13023-016-0545-5/tables/3

# In[ ]:


MECP2_markers = ['CREB1', 'MYCN',"UBE3A", "GRID1", "MECP2"]

#for ISG15: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3587786/


#Select present markers for plotting
present_MECP2 = N_clusters.var_names[N_clusters.var_names.isin(MECP2_markers) == True]
present_MECP2 = [x for x in MECP2_markers if x in present_MECP2] #Keep the order of the markers, better way??

#Print missing markers
print('\nThe following marker genes are missing: ', set(MECP2_markers).difference(set(N_clusters.var_names)))


# In[ ]:


sc.pl.umap(N_clusters, color=present_MECP2)


# In[ ]:


sc.pl.violin(N_clusters, keys=present_MECP2[0:4],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# # Supervised exploration of "Transcriptional Regulation by MECP2" reactome pathway in neurons

# In[ ]:


MECP2_markers = ["MECP2","MOBP","SGK1","RBFOX1","CREB1","MEF2C","PPARG","MIR137","GAMT","PVALB",
                 "MIR132","IRAK1","BDNF","CRH","SST","DLL1","GAD2","GAD1","PTPN1","NOTCH1",
                 "MET","OPRK1","PTPN4","GPRIN1","TRPC3","OPRM1","GRIN2B","SLC2A3","GRIN2A","GRIA2",
                 "FKBP5","PTEN"]

#Select present markers for plotting
present_MECP2 = N_clusters.var_names[N_clusters.var_names.isin(MECP2_markers) == True]
present_MECP2 = [x for x in MECP2_markers if x in present_MECP2] #Keep the order of the markers, better way??

#Print missing markers
print('\nThe following marker genes are missing: ', set(MECP2_markers).difference(set(N_clusters.var_names)))


# In[ ]:


sc.pl.umap(N_clusters, color=present_MECP2, vmin='p1', vmax='p99')


# In[ ]:


sc.pl.violin(N_clusters, keys=present_MECP2[0:3],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(N_clusters, keys=present_MECP2[3:6],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(N_clusters, keys=present_MECP2[6:9],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(N_clusters, keys=present_MECP2[9:12],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(N_clusters, keys=present_MECP2[12:15],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(N_clusters, keys=present_MECP2[15:18],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(N_clusters, keys=present_MECP2[18:21],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(N_clusters, keys=present_MECP2[21:24],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(N_clusters, keys=present_MECP2[24:27],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# # Supervised exploration of "Transcriptional Regulation by MECP2" reactome pathway in infeccted neurons

# In[ ]:


only_infected = N_clusters[N_clusters.obs["infected"].isin(["True"]),:]

MECP2_markers = ["MECP2","MOBP","SGK1","RBFOX1","CREB1","MEF2C","PPARG","MIR137","GAMT","PVALB",
                 "MIR132","IRAK1","BDNF","CRH","SST","DLL1","GAD2","GAD1","PTPN1","NOTCH1",
                 "MET","OPRK1","PTPN4","GPRIN1","TRPC3","OPRM1","GRIN2B","SLC2A3","GRIN2A","GRIA2",
                 "FKBP5","PTEN"]

#Select present markers for plotting
present_MECP2 = only_infected.var_names[only_infected.var_names.isin(MECP2_markers) == True]
present_MECP2 = [x for x in MECP2_markers if x in present_MECP2] #Keep the order of the markers, better way??

#Print missing markers
print('\nThe following marker genes are missing: ', set(MECP2_markers).difference(set(only_infected.var_names)))


# In[ ]:


sc.pl.umap(only_infected, color=present_MECP2, vmin='p1', vmax='p99')


# In[ ]:


sc.pl.violin(only_infected, keys=present_MECP2[0:3],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(only_infected, keys=present_MECP2[3:6],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(only_infected, keys=present_MECP2[6:9],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(only_infected, keys=present_MECP2[9:12],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(only_infected, keys=present_MECP2[12:15],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(only_infected, keys=present_MECP2[15:18],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(only_infected, keys=present_MECP2[18:21],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(only_infected, keys=present_MECP2[21:24],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


# In[ ]:


sc.pl.violin(only_infected, keys=present_MECP2[24:27],
             jitter=False, multi_panel=True, groupby='batch_id', rotation=45)


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




