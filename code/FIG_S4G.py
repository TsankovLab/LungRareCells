import os
import pandas as pd
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt

DATA_DIR = '../data/'
FIG_DIR = '../figures/'

def update_index_cols(df):
    indices = []
    for idx in df.index:
        indices.append(f'{idx} ({df.loc[idx].sum()})')
    df.index = indices
    colnames = []
    for col in df.columns:
        colnames.append(f'{col} ({df.loc[:,col].sum()})')
    df.columns = colnames
    return df

adata = ad.read_h5ad(os.path.join(DATA_DIR, 'hlca_normal_scored.h5ad'))
final_annotations = pd.read_csv(os.path.join(DATA_DIR, 'hlca_normal_annot_after_doublet_removal.csv'), index_col=0)

adata.obs = pd.concat([adata.obs, final_annotations], axis=1)

adata.obs['ann_finest_level_v2'] = adata.obs['ann_finest_level'].copy()
adata.obs['ann_finest_level_v2'] = adata.obs['ann_finest_level_v2'].astype(str)
adata.obs.loc[~adata.obs['ann_finest_level_v2'].isin(['Neuroendocrine', 'Ionocyte', 'Tuft']), 'ann_finest_level_v2'] = 'Other'
adata.obs['ann_finest_level_v2'] = adata.obs['ann_finest_level_v2'].astype('category')

ann_cm = pd.crosstab(index=adata.obs['custom_annot'], columns=adata.obs['ann_finest_level_v2'])
ann_cm = ann_cm.loc[['Neuroendocrine', 'Ionocyte', 'POU2F3+ tuft-like', 'Mature tuft']]
ann_cm = ann_cm.loc[:, ['Neuroendocrine', 'Ionocyte', 'Tuft', 'Other']]
ann_cm = update_index_cols(ann_cm)

ann_cm_normed = ann_cm.div(ann_cm.sum(axis=1), axis=0)

fig = plt.figure()
sns.heatmap(ann_cm_normed, annot=ann_cm, cmap='Oranges', fmt='d')
plt.xlabel('HLCA cell annotations')
plt.ylabel('Marker scoring cell annotations')
plt.savefig(os.path.join(FIG_DIR, 'FIG_S4G_left.pdf'), bbox_inches='tight')
plt.close()

ann_cm = pd.crosstab(index=adata.obs['custom_annot'], columns=adata.obs['tissue'])
ann_cm = ann_cm.loc[['Neuroendocrine', 'Ionocyte', 'POU2F3+ tuft-like', 'Mature tuft']]
ann_cm = update_index_cols(ann_cm)

ann_cm_normed = ann_cm.div(ann_cm.sum(axis=1), axis=0)

fig = plt.figure()
sns.heatmap(ann_cm_normed, annot=ann_cm, cmap='Oranges', fmt='d')
plt.xlabel('HLCA regions sampled')
plt.ylabel('Marker scoring cell annotations')
plt.savefig(os.path.join(FIG_DIR, 'FIG_S4G_right.pdf'), bbox_inches='tight')
plt.close()
