#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot deconvolution results for Figure 2
Author: Yaron Trink
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------
# Setup
# --------------------------------------------------------------------------
sys.path.append("../figure_1")
from _wilms_nD import wilms_Narcs

basepath = "./"
deconvolution_path = os.path.join(basepath, "deconvolution_pca")
os.makedirs(deconvolution_path, mode=0o777, exist_ok=True)

# --------------------------------------------------------------------------
# Load data
# --------------------------------------------------------------------------
wilms = wilms_Narcs(filename="../data/tcga_counts_all_threeArcs.csv")
metadata = pd.read_csv("../data/metadataAll.csv", index_col=0)

# --------------------------------------------------------------------------
# Define plotting function
# --------------------------------------------------------------------------
def plot_deconvolution(score, arcs, nArcs, labels,
                       all_samples=True,
                       save_to_file=False,
                       folder="Desktop",
                       file_end="basic",
                       file_type="svg",
                       title="Hello"):
    """
    Plot deconvolution results in PCA space.
    """
    filename = os.path.join(folder, f"{file_end}.{file_type}")

    score_normal = score.iloc[:-nArcs]
    score_arcs = score.tail(nArcs)

    fig, ax = plt.subplots()

    sz = np.array(labels, dtype=float)
    sz = sz - sz.min()
    sz = (sz / sz.max() * 50) + 5

    if all_samples:
        ax.scatter(score.iloc[:, 0], score.iloc[:, 1],
                   c=sz, s=sz, cmap="bwr", marker="o",
                   linewidths=0.7, edgecolors="black")
    else:
        sz = np.array(labels, dtype=float)
        sz = sz - sz.min()
        sz = (sz / sz.max() * 100) + 5

        ax.scatter(score_normal.iloc[:, 0], score_normal.iloc[:, 1],
                   c=sz, s=sz, cmap="bwr", marker="o")
        ax.scatter(score_arcs.iloc[:, 0], score_arcs.iloc[:, 1],
                   color="red", marker="*")
        ax.legend()

    # Label archetypes
    for idx, arc_name in enumerate(arcs):
        ax.text(score_arcs.iloc[idx, 0], score_arcs.iloc[idx, 1], arc_name)

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(title, fontsize=20)

    if save_to_file:
        os.makedirs(folder, exist_ok=True)
        plt.savefig(filename, bbox_inches="tight")
        print(f"Saved plot: {filename}")

    plt.show()
    plt.close()

# --------------------------------------------------------------------------
# Run deconvolution plots
# --------------------------------------------------------------------------
# Get data from wilms_Narcs instance
score = wilms.get_score()
arcs = wilms.get_arcs()
nArcs = wilms.get_nArcs()

# Cell-type labels to iterate
labels_to_plot = [
    "Cap mesenchyme", "Endothelium", "Loop of Henle", "Podocyte",
    "Proliferating cap mesenchyme", "Proximal tubule", "UB",
    "Macrophage", "Fibroblast", "Myofibroblast", "Renal Vesicle",
    "S-shaped body"
]

for label_name in labels_to_plot:
    print(f"Plotting {label_name}...")
    label_values = metadata[label_name].values
    plot_deconvolution(
        score=score,
        arcs=arcs,
        nArcs=nArcs,
        labels=label_values,
        save_to_file=True,
        folder=deconvolution_path,
        file_end=label_name.replace(" ", "_"),
        file_type="svg",
        title=label_name
    )
