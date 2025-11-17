# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:50:27 2022

@author: trinkya
"""
import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA 
import matplotlib.pyplot as plt
import sys
sys.path.append('C:/Users/trinkya/OneDrive - Bar Ilan University/github/paper2_splicing/src/figure_1/')
from wilms_module_1 import run_pca_on_normalized_counts


def plot_rmats_on_pca_genes2d(
    score,                          
    n_arcs,                          
    arcs,                            
    inc_levels_path,                 
    events_path,                     
    save=False,
    folder='',
    file_type='png'
):
    """
    Plot rMATS inclusion levels on 2D gene-expression PCA space (PC1 vs PC2),
    one figure per event/gene. Last n_arcs rows of `score` are treated as archetypes.

    Parameters
    ----------
    score : array-like or pandas.DataFrame
        Coordinates; only the first two columns are used (PC1, PC2).
        Assumed order: regular samples first, then the `n_arcs` archetype rows at the end.
    n_arcs : int
        Number of archetype rows at the end of `score`.
    arcs : list of str
        Names for the archetypes, length must equal `n_arcs`.
    inc_levels_path : str or Path
        CSV with inclusion-level matrix; rows = events (gene identifiers), columns = samples.
    events_path : str or Path
        CSV with full rMATS event rows; indexed by the same event IDs as inc_levels_path.
    save : bool
        If True, save figures to `folder` with extension `file_type`.
    folder : str
        Output directory.
    file_type : str
        File extension for saved figures, e.g. 'png', 'svg', 'pdf'.
    """
    
    if save:
        os.makedirs(folder, exist_ok=True)

    inc_matrix = pd.read_csv(inc_levels_path, index_col=0)
    events_all = pd.read_csv(events_path, index_col=0)

  
    score = pd.DataFrame(score)
    if score.shape[1] < 2:
        raise ValueError("`score` must have at least 2 columns (PC1, PC2).")

    score.columns = [f"PC{i+1}" for i in range(score.shape[1])]
    if n_arcs <= 0 or n_arcs > len(score):
        raise ValueError("`n_arcs` must be in 1..n_rows(score).")

    if len(arcs) != n_arcs:
        raise ValueError(f"Length of `arcs` ({len(arcs)}) must equal n_arcs ({n_arcs}).")

    score_reg = score.iloc[:-n_arcs, :2]   # PC1, PC2 for regular samples
    score_arcs = score.iloc[-n_arcs:, :2]  # PC1, PC2 for archetypes (at the end)

 
    for g in inc_matrix.index:
        
        if g not in events_all.index:
            continue

        event = events_all.loc[g, :]
        exon_start = event.get('exonStart_0base', 'NA')

        
        sz = inc_matrix.loc[g, :].astype(float).to_numpy()

       
        minv = np.nanmin(sz)
        sz = sz - minv
        maxv = np.nanmax(sz)
        if maxv > 0:
            sz = (sz / maxv * 50.0) + 5.0
        else:
            sz = np.full_like(sz, 5.0, dtype=float)  # all equal, avoid div by zero

        # Plot
        plt.figure()
        plt.scatter(
            score_reg.iloc[:, 0], score_reg.iloc[:, 1],
            s=sz, c=sz, cmap='bwr', linewidths=0.7, edgecolors='black'
        )
        plt.scatter(
            score_arcs.iloc[:, 0], score_arcs.iloc[:, 1],
            color='black', marker='*', s=20
        )

       
        for idx, name in enumerate(arcs):
            plt.text(
                score_arcs.iloc[idx, 0], score_arcs.iloc[idx, 1],
                s=name, ha='right'
            )

        plt.xlabel('PC1', fontsize=15)
        plt.ylabel('PC2', fontsize=15)
        plt.title(f"{g.split('-')[0]} exon start: {exon_start}", fontsize=20)

        if save:
        
            fname = f"{g.split('(')[0]}exon_start_{exon_start}"
            fname = fname.replace(' ', '_')
            fname = fname.replace(':', '').replace('-', '_') + f".{file_type}"
            out_path = os.path.join(folder, fname)
            plt.savefig(out_path, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
            
            





# %%


counts_path = '../data/tcga_counts_all_threeArcs.csv'
expr = pd.read_csv(counts_path, index_col=0)



pca_res = run_pca_on_normalized_counts(
    expr=expr,
    n_components=3,
    arc_cols=None,
    n_arcs=3,
    exclude_archetypes_from_pca=False,
    center=True,
    scale=True,     
    log1p=False,    
    random_state=0,
    svd_solver="auto",
)




PC = pca_res["PC"]                     
arcs = pca_res["archetype_cols"]       
samples_all = list(PC.index)


non_arcs = [s for s in samples_all if s not in arcs]
score_ordered = PC.loc[non_arcs + arcs]  


inc_levels_path = '../data_rmats_and_rmaps/inclusion_level_matrices_merged/inc_level_df_se.csv'
events_path     = '../data_rmats_and_rmaps/inclusion_level_matrices_merged/inc_level_df_se_whole_events.csv'

out_dir = './rmats_pcaGE'
os.makedirs(out_dir, exist_ok=True)

plot_rmats_on_pca_genes2d(
    score=score_ordered,            
    n_arcs=len(arcs),               
    arcs=arcs,                      
    inc_levels_path=inc_levels_path,
    events_path=events_path,
    save=True,
    folder=out_dir,
    file_type='svg'                 
)
















