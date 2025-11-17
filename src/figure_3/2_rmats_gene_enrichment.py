
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get significant genes/events (SE) for ToppGene from rMATS outputs.
Creates per-comparison gene lists (.txt) and event tables (.csv).
Author: yaron trink
"""

import os
import pandas as pd
from all_inclusion_levels import get_events

# --- Config ---
FDR = 0.0
INC_LEVEL = 0.1
OUT_DIR = "./significant_events"
BASE = "../data_rmats_and_rmaps/output_all_three"

os.makedirs(OUT_DIR, mode=0o777, exist_ok=True)

# archetype comparisons
COMPS = [
    ("arc1_vs_arc3", "arc1_arc3"),
    ("arc2_vs_arc3", "arc2_arc3"),
    ("arc1_vs_arc2", "arc1_arc2"),
]

def run_one(comp_dir: str, tag: str) -> None:
    # Read SE results for this comparison
    in_path = os.path.join(BASE, comp_dir, "SE.MATS.JC.txt")
    rmats_output = pd.read_table(in_path, sep="\t")

    # Filter significant events
    sig_events = get_events(rmats_output, fdr_threshold=FDR, inclevel_threshold=INC_LEVEL)

   
    genes = (
        sig_events["GeneID"]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .sort_values()
        .tolist()
    )

    # Save outputs
    pd.DataFrame(genes, columns=["Gene"]).to_csv(
        os.path.join(OUT_DIR, f"{tag}_sig_genes.txt"), index=False
    )
    sig_events.to_csv(
        os.path.join(OUT_DIR, f"{tag}_sig_events.csv"), index=False
    )

    print(f"[{tag}] genes: {len(genes)}  events: {sig_events.shape[0]}")

if __name__ == "__main__":
    for comp_dir, tag in COMPS:
        run_one(comp_dir, tag)
