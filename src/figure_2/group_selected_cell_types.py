#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Group fetal kidney single-cell labels into major categories.
Author: Yaron Trink
"""

import pandas as pd
from pathlib import Path

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------
INPUT_PATH = Path("../deconvolution_files/cell_labels_selected_only_not_grouped.csv")
OUTPUT_PATH = Path("../deconvolution_files/cell_labels_selected_only_grouped.csv")


celltypes_selected = pd.read_csv(INPUT_PATH)["x"]

# -------------------------------------------------------------------------
# Define grouping dictionary
# -------------------------------------------------------------------------
group_dict = {
    # Ureteric bud lineage
    "CNT/PC - proximal UB": "UB",
    "Proximal UB": "UB",
    "Pelvic epithelium - distal UB": "UB",

    # Macrophages
    "Macrophage 1": "Macrophage",
    "Macrophage 2": "Macrophage",
    "Proliferating macrophage": "Macrophage",

    # Fibroblasts
    "Fibroblast 1": "Fibroblast",
    "Fibroblast 2": "Fibroblast",
    "Proliferating fibroblast": "Fibroblast",

    # Myofibroblasts
    "Myofibroblast 1": "Myofibroblast",
    "Myofibroblast 2": "Myofibroblast",
    "Proliferating myofibroblast": "Myofibroblast",

    # Renal vesicle
    "Proliferating distal renal vesicle": "Renal Vesicle",
    "Distal renal vesicle": "Renal Vesicle",
    "Proximal renal vesicle": "Renal Vesicle",

    # Stroma progenitor
    "Stroma progenitor": "Stroma progenitor",
    "Proliferating stroma progenitor": "Stroma progenitor",

    # S-shaped body
    "Medial S shaped body": "S-shaped body",
    "Proximal S shaped body": "S-shaped body",
    "Distal S shaped body": "S-shaped body",
}


new_celltypes = [group_dict.get(cell, cell) for cell in celltypes_selected]

pd.DataFrame(new_celltypes, columns=["celltype"]).to_csv(OUTPUT_PATH, index=False)


print(f"Saved grouped cell types to: {OUTPUT_PATH}")
print(f"Unique groups: {len(set(new_celltypes))}")
print(f"Groups: {sorted(set(new_celltypes))}")




  
    
    
