# LungRareCells
This repository contains custom code to reproduce rare cell detection in HLCA normal dataset and figures in our manuscript, "A deep lung cell atlas reveals cytokine-mediated lineage switching of a rare cell progenitor of the human airway epithelium"

---- 

Download the integrated Human Lung Cell Atlas (HLCA) v1.0 (full) [here](https://data.humancellatlas.org/hca-bio-networks/lung/atlases/lung-v1-0). (https://data.humancellatlas.org/hca-bio-networks/lung/atlases/lung-v1-0)

----

## Installation

To download this repository, run:

```
git clone https://github.com/TsankovLab/LungRareCells.git
```

and download the HLCA .h5ad object into the data folder.

### Installing packages

```
conda env create -f 'LungRareCells/misc/lung_rare_cells.yml'
```

## Detecting rare cells in HLCA normal dataset (full)

Run script `detect_rare_cells.py` with command:

```
python detect_rare_cells.py
```

resulting object, `hlca_normal_scored.h5ad`, will be saved in the `data` folder.

Scripts utilizing `hlca_normal_scored.h5ad` object will also load `hlca_normal_annot_after_doublet_removal.csv` which contain rare celltype annotations after doublet removal.
