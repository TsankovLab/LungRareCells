import os
import anndata as ad
import scanpy as sc
from pathlib import Path
from scanpy._settings import ScanpyConfig

DATA_DIR = '../data/'
FIG_DIR = '../figures/'
ScanpyConfig.figdir = Path(FIG_DIR)

adata = ad.read_h5ad(os.path.join(DATA_DIR, 'hlca_normal_rare_cells_only_custom_genes.h5ad'))

palette = {
    'Ionocyte': '#F8766D',
    'Neuroendocrine': '#0CB702',
    'Mature tuft': 'tab:purple',
    'POU2F3+ tuft-like': '#00B8E7',
}

sc.pl.umap(adata, color=['custom_annot', 'TIP', 'custom_Tuft', 'Neuroendocrine', 'Ionocyte'], palette=palette, cmap='jet', vmin=0, vmax=2, ncols=3, wspace=0.35, show=False, save='_FIG_S4E.pdf')
