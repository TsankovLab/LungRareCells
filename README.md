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
conda env create -f 'LungRareCells/environment/lung_rare_cell.yml'
```
