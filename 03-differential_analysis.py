#!/usr/bin/env python
# coding: utf-8

# <font size="8">Analysis of toxo-infected organoids</font>

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
import decoupler as dc
import os


# In[ ]:


sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100)
sc.logging.print_header()


# In[ ]:


print(datetime.now())


# In[ ]:


msigdb = dc.get_resource('MSigDB')
msigdb


# In[ ]:


curated = msigdb[msigdb['collection'].isin(['reactome_pathways', 'kegg_pathways'])]
curated = curated[~curated.duplicated(['geneset', 'genesymbol'])]
curated.geneset.unique()


# In[ ]:


aggregated = curated[["geneset", "genesymbol"]].groupby("geneset").count()
curated = curated[~curated.geneset.isin(aggregated[aggregated.genesymbol > 100].index.tolist())]


# In[ ]:


data_folder = "./data"
result_folder = './results'


# # Data Load

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


# In[ ]:


N_clusters.layers['scaled'] = sc.pp.scale(N_clusters, copy=True).X


# In[ ]:


N_TOP_GENES = 500


# In[ ]:


P_THRESHOLD = 0.001


# # Differential expression

# ## Differential expression **MECP2** vs  **HA batch** batch neurons

# In[ ]:


sc.tl.rank_genes_groups(N_clusters, 'batch_id', groups=["MECP2"], reference="HA", method='wilcoxon')


# In[ ]:


degs_table = sc.get.rank_genes_groups_df(N_clusters, group="MECP2")
degs_table


# In[ ]:


degs_table.to_csv(os.path.join(result_folder, "neurons_MECP2vsHA_DEGS.csv"))


# In[ ]:


"up-regulated:, ", len(degs_table[(degs_table.pvals_adj < 0.01) & (degs_table.logfoldchanges > 0)])


# In[ ]:


"down-regulated:, ", len(degs_table[(degs_table.pvals_adj < 0.01) & (degs_table.logfoldchanges < 0)])


# In[ ]:


rank = pd.DataFrame()
rank["score"] = -np.sign(degs_table.set_index("names").logfoldchanges) * np.log10(degs_table.set_index("names").pvals_adj + 1.E-300)
rank.sort_values(by="score")


# In[ ]:


up = rank.sort_values(by="score", ascending=False).head(100).index.tolist()
down = rank.sort_values(by="score", ascending=True).head(100).index.tolist()
sns.set()
sc.pl.heatmap(N_clusters, {"up": up, "down": down}, groupby='batch_id', swap_axes=False, layer="scaled", cmap="viridis")


# In[ ]:


pthreshold = P_THRESHOLD
n_up = N_TOP_GENES
n_bottom = 0

estimate, pvals = dc.run_ora(mat=rank.T, net=curated, source='geneset', target='genesymbol', verbose=False, n_up=n_up, n_bottom=n_bottom)
ora = pd.concat([estimate, pvals]).T
ora.columns = ["estimate", "pvals"]
ora["log10_pval"] = -np.log10(ora["pvals"])
ora = ora.reset_index()

sns.set(font_scale=0.8) 
height = 3+len(ora[ora.pvals < pthreshold])*0.2
ax = sns.relplot(data=ora[ora.pvals < pthreshold].sort_values(by="log10_pval", ascending=False), x="log10_pval", y="source", 
                 height=height, aspect=15/height, legend="auto")


# In[ ]:


ora[ora.source.str.contains("MECP2")]


# In[ ]:


ora[ora.pvals < 1].to_csv(os.path.join(result_folder, "neurons_MECP2vsHA_funct_TOP{}.csv".format(N_TOP_GENES)))


# In[ ]:


pthreshold = P_THRESHOLD
n_up = 0
n_bottom = N_TOP_GENES

estimate, pvals = dc.run_ora(mat=rank.T, net=curated, source='geneset', target='genesymbol', verbose=False, n_up=n_up, n_bottom=n_bottom)
ora = pd.concat([estimate, pvals]).T
ora.columns = ["estimate", "pvals"]
ora["log10_pval"] = -np.log10(ora["pvals"])
ora = ora.reset_index()

sns.set(font_scale=0.8)
height = 3+len(ora[ora.pvals < pthreshold])*0.2
ax = sns.relplot(data=ora[ora.pvals < pthreshold].sort_values(by="log10_pval", ascending=False), x="log10_pval", y="source", 
                 height=height, aspect=15/height, legend="auto")


# In[ ]:


ora[ora.source.str.contains("MECP2")]


# In[ ]:


ora[ora.pvals < 1].to_csv(os.path.join(result_folder, "neurons_MECP2vsHA_funct_BOTTOM{}.csv".format(N_TOP_GENES)))


# In[ ]:


pthreshold = P_THRESHOLD

estimate, norm, pvals = dc.run_gsea(mat=rank.T, net=curated, source='geneset', target='genesymbol', verbose=False)
ora = pd.concat([estimate, norm, pvals]).T
ora.columns = ["estimate", "norm", "pvals"]
ora["log10_pval"] = -np.log10(ora["pvals"])
ora = ora.reset_index()

sns.set(font_scale=0.8)
maxv, minv = ora["log10_pval"].max(), ora["log10_pval"].min()
height = 3+len(ora[ora.pvals < pthreshold])*0.2
ax = sns.relplot(data=ora[ora.pvals < pthreshold].sort_values(by="norm", ascending=False), x="norm", y="source", size="log10_pval", 
                 sizes=(10,200), size_norm=(minv, maxv), 
                 hue="log10_pval", hue_norm=(minv, maxv), palette="copper", 
                 height=height, aspect=15/height, legend="auto")


# In[ ]:


ora[ora.source.str.contains("MECP2")]


# In[ ]:


ora[ora.pvals < 1].to_csv(os.path.join(result_folder, "neurons_MECP2vsHA_funct_GSEA.csv"))


# ## Differential expression **MECP2 batch** vs **CNT batch** **neurons**

# In[ ]:


sc.tl.rank_genes_groups(N_clusters, 'batch_id', groups=["MECP2"], reference="CNT", method='wilcoxon')


# In[ ]:


degs_table = sc.get.rank_genes_groups_df(N_clusters, group="MECP2")
degs_table


# In[ ]:


degs_table.to_csv(os.path.join(result_folder, "neurons_MECP2vsCNT_DEGS.csv"))


# In[ ]:


"up-regulated:, ", len(degs_table[(degs_table.pvals_adj < 0.01) & (degs_table.logfoldchanges > 0)])


# In[ ]:


"down-regulated:, ", len(degs_table[(degs_table.pvals_adj < 0.01) & (degs_table.logfoldchanges < 0)])


# In[ ]:


rank = pd.DataFrame()
rank["score"] = -np.sign(degs_table.set_index("names").logfoldchanges) * np.log10(degs_table.set_index("names").pvals_adj + 1.E-300)
rank.sort_values(by="score")


# In[ ]:


up = rank.sort_values(by="score", ascending=False).head(100).index.tolist()
down = rank.sort_values(by="score", ascending=True).head(100).index.tolist()
sns.set()
sc.pl.heatmap(N_clusters, {"up": up, "down": down}, groupby='batch_id', swap_axes=False, layer="scaled", cmap="viridis")


# In[ ]:


pthreshold = P_THRESHOLD
n_up = N_TOP_GENES
n_bottom = 0

estimate, pvals = dc.run_ora(mat=rank.T, net=curated, source='geneset', target='genesymbol', verbose=False, n_up=n_up, n_bottom=n_bottom)
ora = pd.concat([estimate, pvals]).T
ora.columns = ["estimate", "pvals"]
ora["log10_pval"] = -np.log10(ora["pvals"])
ora = ora.reset_index()

sns.set(font_scale=0.8) 
height = 3+len(ora[ora.pvals < pthreshold])*0.2
ax = sns.relplot(data=ora[ora.pvals < pthreshold].sort_values(by="log10_pval", ascending=False), x="log10_pval", y="source", 
                 height=height, aspect=15/height, legend="auto")


# In[ ]:


ora[ora.source.str.contains("MECP2")]


# In[ ]:


ora[ora.pvals < 1].to_csv(os.path.join(result_folder, "neurons_MECP2vsCNT_funct_TOP{}.csv".format(N_TOP_GENES)))


# In[ ]:


pthreshold = P_THRESHOLD
n_up = 0
n_bottom = N_TOP_GENES

estimate, pvals = dc.run_ora(mat=rank.T, net=curated, source='geneset', target='genesymbol', verbose=False, n_up=n_up, n_bottom=n_bottom)
ora = pd.concat([estimate, pvals]).T
ora.columns = ["estimate", "pvals"]
ora["log10_pval"] = -np.log10(ora["pvals"])
ora = ora.reset_index()

sns.set(font_scale=0.8)
height = 3+len(ora[ora.pvals < pthreshold])*0.2
ax = sns.relplot(data=ora[ora.pvals < pthreshold].sort_values(by="log10_pval", ascending=False), x="log10_pval", y="source", 
                 height=height, aspect=15/height, legend="auto")


# In[ ]:


ora[ora.pvals < 1].to_csv(os.path.join(result_folder, "neurons_MECP2vsCNT_funct_BOTTOM{}.csv".format(N_TOP_GENES)))


# In[ ]:


ora[ora.source.str.contains("MECP2")]


# In[ ]:


pthreshold = P_THRESHOLD

estimate, norm, pvals = dc.run_gsea(mat=rank.T, net=curated, source='geneset', target='genesymbol', verbose=False)
ora = pd.concat([estimate, norm, pvals]).T
ora.columns = ["estimate", "norm", "pvals"]
ora["log10_pval"] = -np.log10(ora["pvals"])
ora = ora.reset_index()

sns.set(font_scale=0.8)
maxv, minv = ora["log10_pval"].max(), ora["log10_pval"].min()
height = 3+len(ora[ora.pvals < pthreshold])*0.2
ax = sns.relplot(data=ora[ora.pvals < pthreshold].sort_values(by="norm", ascending=False), x="norm", y="source", size="log10_pval", 
                 sizes=(10,200), size_norm=(minv, maxv), 
                 hue="log10_pval", hue_norm=(minv, maxv), palette="copper", 
                 height=height, aspect=15/height, legend="auto")


# In[ ]:


ora[ora.pvals < 1].to_csv(os.path.join(result_folder, "neurons_MECP2vsCNT_funct_GSEA.csv"))


# In[ ]:


ora[ora.source.str.contains("MECP2")]


# ## Differential expression **HA batch** vs **CNT batch neurons**

# In[ ]:


sc.tl.rank_genes_groups(N_clusters, 'batch_id', groups=["HA"], reference="CNT", method='wilcoxon')


# In[ ]:


degs_table = sc.get.rank_genes_groups_df(N_clusters, group="HA")
degs_table


# In[ ]:


degs_table.to_csv(os.path.join(result_folder, "neurons_HAvsCNT_DEGS.csv"))


# In[ ]:


"up-regulated:, ", len(degs_table[(degs_table.pvals_adj < 0.01) & (degs_table.logfoldchanges > 0)])


# In[ ]:


"down-regulated:, ", len(degs_table[(degs_table.pvals_adj < 0.01) & (degs_table.logfoldchanges < 0)])


# In[ ]:


rank = pd.DataFrame()
rank["score"] = -np.sign(degs_table.set_index("names").logfoldchanges) * np.log10(degs_table.set_index("names").pvals_adj + 1.E-300)
rank.sort_values(by="score")


# In[ ]:


up = rank.sort_values(by="score", ascending=False).head(100).index.tolist()
down = rank.sort_values(by="score", ascending=True).head(100).index.tolist()
sns.set()
sc.pl.heatmap(N_clusters, {"up": up, "down": down}, groupby='batch_id', swap_axes=False, layer="scaled", cmap="viridis")


# In[ ]:


pthreshold = P_THRESHOLD
n_up = N_TOP_GENES
n_bottom = 0

estimate, pvals = dc.run_ora(mat=rank.T, net=curated, source='geneset', target='genesymbol', verbose=False, n_up=n_up, n_bottom=n_bottom)
ora = pd.concat([estimate, pvals]).T
ora.columns = ["estimate", "pvals"]
ora["log10_pval"] = -np.log10(ora["pvals"])
ora = ora.reset_index()

sns.set(font_scale=0.8) 
height = 3+len(ora[ora.pvals < pthreshold])*0.2
ax = sns.relplot(data=ora[ora.pvals < pthreshold].sort_values(by="log10_pval", ascending=False), x="log10_pval", y="source", 
                 height=height, aspect=15/height, legend="auto")


# In[ ]:


ora[ora.source.str.contains("MECP2")]


# In[ ]:


ora[ora.pvals < 1].to_csv(os.path.join(result_folder, "neurons_HAvsCNT_funct_TOP{}.csv".format(N_TOP_GENES)))


# In[ ]:


pthreshold = P_THRESHOLD
n_up = 0
n_bottom = N_TOP_GENES

estimate, pvals = dc.run_ora(mat=rank.T, net=curated, source='geneset', target='genesymbol', verbose=False, n_up=n_up, n_bottom=n_bottom)
ora = pd.concat([estimate, pvals]).T
ora.columns = ["estimate", "pvals"]
ora["log10_pval"] = -np.log10(ora["pvals"])
ora = ora.reset_index()

sns.set(font_scale=0.8)
height = 3+len(ora[ora.pvals < pthreshold])*0.2
ax = sns.relplot(data=ora[ora.pvals < pthreshold].sort_values(by="log10_pval", ascending=False), x="log10_pval", y="source", 
                 height=height, aspect=15/height, legend="auto")


# In[ ]:


ora[ora.source.str.contains("MECP2")]


# In[ ]:


ora[ora.pvals < 1].to_csv(os.path.join(result_folder, "neurons_HAvsCNT_funct_BOTTOM{}.csv".format(N_TOP_GENES)))


# In[ ]:


pthreshold = P_THRESHOLD

estimate, norm, pvals = dc.run_gsea(mat=rank.T, net=curated, source='geneset', target='genesymbol', verbose=False)
ora = pd.concat([estimate, norm, pvals]).T
ora.columns = ["estimate", "norm", "pvals"]
ora["log10_pval"] = -np.log10(ora["pvals"])
ora = ora.reset_index()

sns.set(font_scale=0.8)
maxv, minv = ora["log10_pval"].max(), ora["log10_pval"].min()
height = 3+len(ora[ora.pvals < pthreshold])*0.2
ax = sns.relplot(data=ora[ora.pvals < pthreshold].sort_values(by="norm", ascending=False), x="norm", y="source", size="log10_pval", 
                 sizes=(10,200), size_norm=(minv, maxv), 
                 hue="log10_pval", hue_norm=(minv, maxv), palette="copper", 
                 height=height, aspect=15/height, legend="auto")


# In[ ]:


ora[ora.pvals < 1].to_csv(os.path.join(result_folder, "neurons_HAvsCNT_funct_GSEA.csv"))


# In[ ]:


ora[ora.source.str.contains("MECP2")]


# # REACTOME_TRANSCRIPTIONAL_REGULATION_BY_MECP2 score

# In[ ]:


MECP2_list = curated[curated.geneset == "REACTOME_TRANSCRIPTIONAL_REGULATION_BY_MECP2"].genesymbol.tolist()
len(MECP2_list)


# In[ ]:


sc.tl.score_genes(adata, MECP2_list, score_name='MECP2_ptw')


# In[ ]:


sc.settings.set_figure_params(dpi=100, facecolor='white')
sc.pl.umap(adata, color=['MECP2_ptw', 'celltype'], vmin=0., vmax='p99')


# In[ ]:


N_clusters = adata[adata.obs["celltype"].isin(['N1', 'N2', 'N3'])].copy()
N_clusters


# In[ ]:


N_clusters.layers['scaled'] = sc.pp.scale(N_clusters, copy=True).X


# In[ ]:


fig, ax = plt.subplots(figsize=(4,3))
sc.pl.violin(N_clusters, ['MECP2_ptw'], groupby='batch_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()


# In[ ]:


sc.settings.set_figure_params(dpi=100, facecolor='white')
fig, ax = plt.subplots(figsize=(4,3))
sc.pl.violin(N_clusters, ['MECP2_ptw'], groupby='batch_id', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()


# In[ ]:


fig, ax = plt.subplots(figsize=(4,3))
sc.pl.violin(adata[adata.obs["celltype"].isin(['CyclingProg', 'vRG_oRG', 'vRG_oRG2'])], ['MECP2_ptw'], groupby='batch_infect', multi_panel=True, jitter=False, log=False, rotation=90., ax=ax)
plt.show()


# # 10. Transcription factor activity inference

# In[ ]:


net = dc.get_dorothea(organism='human', levels=['A','B','C'])
net


# In[ ]:


dc.run_mlm(mat=adata, net=net, source='source', target='target', weight='weight', verbose=True, use_raw=False)


# In[ ]:


acts = dc.get_acts(adata, obsm_key='mlm_estimate')
acts


# In[ ]:


mean_acts = dc.summarize_acts(acts, groupby='batch_id', min_std=0.3)
mean_acts


# In[ ]:


mean_acts = dc.summarize_acts(acts, groupby="celltype", min_std=0.5)
mean_acts


# In[ ]:


mean_acts = dc.summarize_acts(acts[acts.obs["celltype"].isin(['N1', 'N2', 'N3'])], groupby="batch_infect", min_std=0.2)
sns.clustermap(mean_acts.T, center=0, vmax=5, cmap='coolwarm', figsize=(4, 15))
plt.show()


# In[ ]:


mean_acts = dc.summarize_acts(acts[acts.obs["celltype"].isin(['N_UPR', 'N_UPR2', 'N_UPR3'])], groupby="batch_infect", min_std=0.2)
sns.clustermap(mean_acts.T, center=0, vmax=5, cmap='coolwarm', figsize=(4, 15))
plt.show()


# In[ ]:


mean_acts = dc.summarize_acts(acts[acts.obs["celltype"].isin(['CyclingProg', 'vRG_oRG', 'vRG_oRG2'])], groupby="batch_infect", min_std=0.2)
sns.clustermap(mean_acts.T, center=0, vmax=5, cmap='coolwarm', figsize=(4, 30))
plt.show()


# In[ ]:


mean_acts = dc.summarize_acts(acts[(acts.obs.infected == "True") & (acts.obs["celltype"].isin(["N1", "N2", "N3"]))], groupby='batch_id', min_std=0.3)
tf = mean_acts.columns
mean_acts = dc.summarize_acts(acts[acts.obs["celltype"].isin(['N1', 'N2', 'N3'])], groupby="batch_infect", min_std=0.0)
mean_acts = mean_acts[tf]
sns.clustermap(mean_acts.T, center=0, vmax=5, cmap='coolwarm', figsize=(4, 12))
plt.show()


# In[ ]:


sc.settings.set_figure_params(dpi=100, facecolor='white')
sc.pl.umap(acts, color=['HSF1', 'celltype'], cmap='coolwarm', vcenter=0)


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




