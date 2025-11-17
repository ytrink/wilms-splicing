
# -*- coding: utf-8 -*-
"""
Run inclusion-level extraction from rMATS outputs (SE & MXE) for all archetype
pairwise comparisons, then merge and deduplicate results.

@author: trinkya
"""

from __future__ import annotations

import os
import numpy as np
import pandas as pd



# import custom functions for rMATS results analysis
from all_inclusion_levels import (
    run_inc_level_analysis_all_samples_SE,
    run_inc_level_analysis_all_samples_MXE,
    inc_level_merge,
)


FDR_THRESHOLD = 0.0
INCLEVEL_THRESHOLD = 0.2


colnames_target = list(pd.read_csv("../data/tcga_counts_all_threeArcs.csv", index_col=0).columns)

BASE = "../data_rmats_and_rmaps"
OUT_INDIV = f"{BASE}/inclusion_level_matrices"
OUT_MERGED = f"{BASE}/inclusion_level_matrices_merged"
os.makedirs(OUT_INDIV, mode=0o777, exist_ok=True)
os.makedirs(OUT_MERGED, mode=0o777, exist_ok=True)


ALL_EVENTS = {
    "SE": pd.read_table(f"{BASE}/output_all_samples_comparison/SE.MATS.JC.txt", sep="\t"),
    "MXE": pd.read_table(f"{BASE}/output_all_samples_comparison/MXE.MATS.JC.txt", sep="\t"),
}


COMPARISONS = {
    "Arc1-Arc3": "arc1_vs_arc3",
    "Arc2-Arc3": "arc2_vs_arc3",
    "Arc1-Arc2": "arc1_vs_arc2",
}


def run_and_save_one(comparison: str, comp_dir: str, etype: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run inclusion-level extraction for a single comparison and event type, then save CSVs.
    Returns (incLevels, wholeEvents).
    """
    assert etype in ("SE", "MXE")

    in_path = f"{BASE}/output_all_three/{comp_dir}/{etype}.MATS.JC.txt"
    rmats_output = pd.read_table(in_path, sep="\t")

    if etype == "SE":
        inc_df, whole_df = run_inc_level_analysis_all_samples_SE(
            rmats_output,
            ALL_EVENTS["SE"],
            comparison=comparison,
            fdr_threshold=FDR_THRESHOLD,
            inclevel_threshold=INCLEVEL_THRESHOLD,
        )
    else:  # MXE
        inc_df, whole_df = run_inc_level_analysis_all_samples_MXE(
            rmats_output,
            ALL_EVENTS["MXE"],
            comparison=comparison,
            fdr_threshold=FDR_THRESHOLD,
            inclevel_threshold=INCLEVEL_THRESHOLD,
        )

    inc_path = f"{OUT_INDIV}/{etype}_{comp_dir.replace('_vs_', '_')}.csv"
    whole_path = f"{OUT_INDIV}/{etype}_{comp_dir.replace('_vs_', '_')}_whole_events.csv"
    inc_df.to_csv(inc_path)
    whole_df.to_csv(whole_path)

    print(f"[{etype} | {comparison}] events: {inc_df.shape[0]}  -> {inc_path}")
    return inc_df, whole_df


def merge_three(etype: str, comps: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load three per-comparison results for a given event type and merge
    keeping at most one event per gene (via inc_level_merge).
    Returns (incLevels_merged, wholeEvents_merged).
    """
    # Inclusion levels
    d1 = pd.read_table(f"{OUT_INDIV}/{etype}_arc2_arc3.csv", sep=",", index_col=0)
    d2 = pd.read_table(f"{OUT_INDIV}/{etype}_arc1_arc2.csv", sep=",", index_col=0)
    d3 = pd.read_table(f"{OUT_INDIV}/{etype}_arc1_arc3.csv", sep=",", index_col=0)

    d_int = inc_level_merge(d1, d2)
    d_merged = inc_level_merge(d_int, d3)

    # Whole event rows
    w1 = pd.read_table(f"{OUT_INDIV}/{etype}_arc2_arc3_whole_events.csv", sep=",", index_col=0)
    w2 = pd.read_table(f"{OUT_INDIV}/{etype}_arc1_arc2_whole_events.csv", sep=",", index_col=0)
    w3 = pd.read_table(f"{OUT_INDIV}/{etype}_arc1_arc3_whole_events.csv", sep=",", index_col=0)

    w_int = inc_level_merge(w1, w2)
    w_merged = inc_level_merge(w_int, w3)

    # Deduplicate (keep first) and align whole-events to inc rows
    dup_mask = d_merged.duplicated()
    keep_idx = np.where(~dup_mask)[0]
    d_merged = d_merged.drop_duplicates()
    w_merged = w_merged.iloc[keep_idx, :]

    # Sanity: align index
    w_merged = w_merged.reindex(d_merged.index)

    # Save merged
    d_out = f"{OUT_MERGED}/inc_level_df_{etype.lower()}.csv"
    w_out = f"{OUT_MERGED}/inc_level_df_{etype.lower()}_whole_events.csv"
    d_merged.to_csv(d_out)
    w_merged.to_csv(w_out)
    print(f"[{etype}] merged: {d_merged.shape[0]} -> {d_out}")

    return d_merged, w_merged


def safe_assign_target_colnames(df: pd.DataFrame, target: list[str]) -> pd.DataFrame:
    """
    Keep your original behavior of renaming columns to the first 136 sample names,
    but only if sizes match. Otherwise, keep as-is and warn.
    """
    n = min(len(target), df.shape[1])
    if n == df.shape[1]:
        df.columns = target[:n]
    else:
        print(f"[warn] Column count ({df.shape[1]}) != len(target subset) ({n}); leaving columns unchanged.")
    return df


# run for SE and MXE

for comp_name, comp_dir in COMPARISONS.items():
    _se, _se_w = run_and_save_one(comp_name, comp_dir, "SE")
    _mxe, _mxe_w = run_and_save_one(comp_name, comp_dir, "MXE")


# SE merged
dfSE, dfSE_w = merge_three("SE", COMPARISONS)

# MXE merged
dfMXE, dfMXE_w = merge_three("MXE", COMPARISONS)

# Prefer MXE over SE when a gene has both
print(f"MXE + SE before preference: {dfMXE.shape[0] + dfSE.shape[0]}")
dfFinal = inc_level_merge(dfMXE, dfSE)           # keeps MXE rows; drops SE if same gene appears
dfFinal = safe_assign_target_colnames(dfFinal, colnames_target[:136])
print(f"Final merged (unique genes): {dfFinal.shape[0]}")


dfFinal_w = inc_level_merge(dfSE_w, dfMXE_w).reindex(dfFinal.index)


dup_mask = dfFinal.duplicated()
keep_idx = np.where(~dup_mask)[0]
dfFinal = dfFinal.drop_duplicates()
dfFinal_w = dfFinal_w.iloc[keep_idx, :]
missing = len(set(dfFinal.index).difference(set(dfFinal_w.index)))
print(f"Missing whole-event rows after align: {missing}")

# Save final merged (MXE+SE)
dfFinal.to_csv(f"{OUT_MERGED}/inc_level_df_mxe_and_se.csv")
dfFinal2 = dfFinal_w.join(dfFinal)  # main results table with event metadata
dfFinal2.to_csv(f"{OUT_MERGED}/inc_level_df_mxe_and_se_whole_events.csv")

# ----------------------------- SE-only (for PCA) ---------------------------- #


dfSE_only = dfSE.drop_duplicates()
dup_mask = dfSE.duplicated()
keep_idx = np.where(~dup_mask)[0]
dfSE_w_only = dfSE_w.iloc[keep_idx, :]

dfSE_only.to_csv(f"{OUT_MERGED}/inc_level_df_se.csv")
dfSE_w_only.to_csv(f"{OUT_MERGED}/inc_level_df_se_whole_events.csv")

print("All inclusion-level matrices generated and merged successfully.")
