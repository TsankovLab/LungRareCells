import os
import pandas as pd
import anndata as ad
import scanpy as sc
from pathlib import Path
from scanpy._settings import ScanpyConfig
from matplotlib.colors import LinearSegmentedColormap

colors = ["blue", "red"]
cmap = LinearSegmentedColormap.from_list("BlueRed", colors)

DATA_DIR = '../data/'
FIG_DIR = '../figures/'
ScanpyConfig.figdir = Path(FIG_DIR)

adata = ad.read_h5ad(os.path.join(DATA_DIR, 'hlca_normal_scored.h5ad'))
final_annotations = pd.read_csv(os.path.join(DATA_DIR, 'hlca_normal_annot_after_doublet_removal.csv'), index_col=0)

adata.obs = pd.concat([adata.obs, final_annotations], axis=1)

plt_gene_scores = {
        'Original Tuft Geneset': ['CPB1', 'HCK', 'HOTAIRM1', 'PBXIP1'],
        'Common Genes': ['BMX', 'GNG13', 'RGS13', 'SH2D6', 'TRPM5'],
        'New Tuft Geneset': ['ALKAL1',
                          'POU2AF2',
                          'CHAT',
                          'GNAT3',
                          'HTR3C',
                          'HTR3E',
                          'IL17RB',
                          'IL25',
                          'TAS1R3']
}

subset_ctypes = ['Ionocyte', 'Neuroendocrine', 'Mature tuft', 'POU2F3+ tuft-like']

adata.obs['gene_exp_annot'] = adata.obs['ann_level_2'].copy()
adata.obs['gene_exp_annot'] = adata.obs['gene_exp_annot'].astype(str)
adata.obs.loc[adata.obs['original_annot'].isin(subset_ctypes), 'gene_exp_annot'] = adata.obs['original_annot']

sc.pl.dotplot(adata, plt_gene_scores, groupby='gene_exp_annot', cmap=cmap, use_raw=False, figsize=(10,6), dot_max=1, save='FIG_S4F.pdf')
