import os
import json
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from pathlib import Path

DATA_DIR = '../data/'
MISC_DIR = '../misc/'

gene_scores = json.load(open(os.path.join(MISC_DIR, 'gene_list.json'), 'r'))

adata = ad.read_h5ad(os.path.join(DATA_DIR, 'hlca_full.h5ad'))
adata.var.index = adata.var['feature_name']

for celltype, genes in gene_scores.items():
    sc.tl.score_genes(adata, genes, n_bins=24, ctrl_size=100, score_name=celltype, use_raw=False)

# celltype_orig - original markers
adata.obs['celltype_orig'] = 'Other'
adata.obs.loc[adata.obs['Neuroendocrine'] > 1, 'celltype_orig'] = 'Neuroendocrine'
adata.obs.loc[adata.obs['TIP'] > 1.6, 'celltype_orig'] = 'POU2F3+ tuft-like'
adata.obs.loc[adata.obs['Ionocyte'] > 1, 'celltype_orig'] = 'Ionocyte'
adata.obs.loc[adata.obs['Tuft'] > 1, 'celltype_orig'] = 'Mature tuft'

# celltype_custom - custom markers
adata.obs['celltype_custom'] = 'Other'
adata.obs.loc[adata.obs['Neuroendocrine'] > 1, 'celltype_custom'] = 'Neuroendocrine'
adata.obs.loc[adata.obs['TIP'] > 1.6, 'celltype_custom'] = 'POU2F3+ tuft-like'
adata.obs.loc[adata.obs['Ionocyte'] > 1, 'celltype_custom'] = 'Ionocyte'
adata.obs.loc[adata.obs['custom_Tuft'] > 0.5, 'celltype_custom'] = 'Mature tuft'

# subset normal
adata = adata[adata.obs['disease'] == 'normal']
adata.write_h5ad(os.path.join(DATA_DIR, 'hlca_normal_scored.h5ad'))
