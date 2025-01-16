import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['figure.dpi'] = 300

DATA_DIR = '../data/'
FIG_DIR = '../figures/'
os.makedirs(FIG_DIR, exist_ok=True)

sig_thresh = -np.log10(0.05)    # significant threshold

tuft_vs_iono = pd.read_csv(os.path.join(DATA_DIR, 'tuft_vs_ionocyte_gsea.tsv'), sep='\t')
tuft_vs_iono = tuft_vs_iono.iloc[:4]
tuft_vs_iono = tuft_vs_iono.sort_values(['p-value'], ascending=False)
tuft_vs_iono = tuft_vs_iono[['Gene Set Name', 'FDR q-value']]
tuft_vs_iono['enriched_in'] = 'tuft'

iono_vs_tuft = pd.read_csv(os.path.join(DATA_DIR, 'ionocyte_vs_tuft_gsea.tsv'), sep='\t')
iono_vs_tuft = iono_vs_tuft.iloc[:4][['Gene Set Name', 'FDR q-value']]
iono_vs_tuft['enriched_in'] = 'ionocyte'

df = pd.concat([tuft_vs_iono, iono_vs_tuft])
df['log_q_value'] = np.log10(df['FDR q-value'])
df['geneset_names'] = df['Gene Set Name'] + '_' + df['enriched_in']

fig, ax = plt.subplots(1, 1, figsize=(15,15))
ax.barh(df.loc[df['enriched_in'] == 'ionocyte', 'geneset_names'].values, df.loc[df['enriched_in'] == 'ionocyte', 'log_q_value'].values, color = 'c')
ax.axvline(-sig_thresh, color='black', linewidth=0.8, linestyle='--')
plt.savefig(os.path.join(FIG_DIR, 'FIG_3A_left.pdf'), bbox_inches='tight')
plt.close()

fig, ax = plt.subplots(1, 1, figsize=(15,15))
ax.barh(df.loc[df['enriched_in'] == 'tuft', 'geneset_names'].values, -df.loc[df['enriched_in'] == 'tuft', 'log_q_value'].values, color = 'y')
ax.axvline(sig_thresh, color='black', linewidth=0.8, linestyle='--')
ax.set_xticks([0, 1, 2, 3, 4])
plt.savefig(os.path.join(FIG_DIR, 'FIG_3A_right.pdf'), bbox_inches='tight')
plt.close()
