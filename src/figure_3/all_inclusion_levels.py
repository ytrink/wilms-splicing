# -*- coding: utf-8 -*-


"""
Utilities to extract and merge rMATS splicing event inclusion levels
across all samples, and to look up events discovered in pairwise
comparisons (SE and MXE).

@author: trinkya
"""

from __future__ import annotations

import datetime
from typing import List, Tuple, Dict, Iterable, Union, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ----------------------------- Helpers ------------------------------------- #

def _parse_inc_levels(cell: Union[str, float, int]) -> List[float]:
    """
    Parse a comma-separated inclusion-level cell from rMATS into floats.
    Converts 'NA' to 0, and gracefully handles empty strings / NaNs.
    """
    if isinstance(cell, float) and np.isnan(cell):
        return []
    if not isinstance(cell, str):
        cell = str(cell)
    cell = cell.strip()
    if not cell:
        return []
    vals = []
    for x in cell.split(','):
        x = x.strip()
        if x.upper() == 'NA' or x == '':
            vals.append(0.0)
        else:
            try:
                vals.append(float(x))
            except ValueError:
                vals.append(0.0)
    return vals


def return_inc_levels(mats: pd.DataFrame, rownumber: Union[int, np.integer, np.ndarray]) -> Tuple[List[float], List[float]]:
    """
    Get inclusion levels for all samples from a particular splicing event row.

    Parameters
    ----------
    mats : DataFrame with rMATS output columns 'IncLevel1' and 'IncLevel2'
    rownumber : row index (int or array with a single element)

    Returns
    -------
    inc1 : list of floats (Inclusion Level 1 across samples)
    inc2 : list of floats (Inclusion Level 2 across samples)
    """
    # robust index handling
    if isinstance(rownumber, (np.ndarray, list, tuple)) and len(rownumber) == 1:
        rownumber = int(rownumber[0])
    elif not isinstance(rownumber, (int, np.integer)):
        rownumber = int(rownumber)

    inc1 = _parse_inc_levels(mats.loc[rownumber, 'IncLevel1'])
    inc2 = _parse_inc_levels(mats.loc[rownumber, 'IncLevel2'])
    return inc1, inc2


def _match_event_row(mats_output: pd.DataFrame, event: pd.Series, cols_to_match: Iterable[str]) -> np.ndarray:
    """
    Row-wise exact match on a subset of columns, treating NaNs as equal.
    Returns the matching index values as a numpy array.
    """
    A = mats_output.loc[:, cols_to_match].copy()
    b = event.loc[cols_to_match].copy()

    # Treat NaNs as equal by filling with a sentinel string
    sentinel = '__NA__'
    A_f = A.fillna(sentinel)
    b_f = b.fillna(sentinel)

    mask = A_f.eq(b_f).all(axis=1)
    return A.index[mask].values


def get_splicing_event_index_SE(mats_output: pd.DataFrame, event: pd.Series) -> np.ndarray:
    """
    Find index of an SE event in the large rMATS table.
    """
    cols_to_match = [
        'GeneID', 'geneSymbol', 'chr', 'strand',
        'exonStart_0base', 'exonEnd', 'upstreamES', 'upstreamEE',
        'downstreamES', 'downstreamEE'
    ]
    return _match_event_row(mats_output, event, cols_to_match)


def get_splicing_event_index_MXE(mats_output: pd.DataFrame, event: pd.Series) -> np.ndarray:
    """
    Find index of an MXE event in the large rMATS table.
    """
    cols_to_match = [
        'GeneID', 'geneSymbol', 'chr', 'strand',
        '1stExonStart_0base', '1stExonEnd', '2ndExonStart_0base', '2ndExonEnd',
        'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE'
    ]
    return _match_event_row(mats_output, event, cols_to_match)


def get_splicing_event_inc_levels(mats: pd.DataFrame, idx: Union[int, np.ndarray]) -> Tuple[List[float], List[float]]:
    """
    Convenience wrapper returning inclusion levels for the (first) matched row.
    """
    if isinstance(idx, (np.ndarray, list, tuple)):
        if len(idx) == 0:
            raise IndexError("No matching event row found.")
        idx = int(idx[0])
    elif not isinstance(idx, (int, np.integer)):
        idx = int(idx)
    return return_inc_levels(mats, rownumber=idx)


def get_events(rmats_output: pd.DataFrame, fdr_threshold: float = 0.0, inclevel_threshold: float = 0.2) -> pd.DataFrame:
    """
    Filter rMATS comparison results to events passing FDR and |IncLevelDifference|.
    """
    return rmats_output.loc[
        (rmats_output['FDR'] <= fdr_threshold) &
        (rmats_output['IncLevelDifference'].abs() >= inclevel_threshold),
        :
    ]


def get_splicing_event_exon_locations(mats: pd.DataFrame, idx: Union[int, np.ndarray]) -> pd.Series:
    """
    Extract the full row for the (first) matched splicing event.
    """
    if isinstance(idx, (np.ndarray, list, tuple)):
        if len(idx) == 0:
            raise IndexError("No matching event row found.")
        idx = int(idx[0])
    elif not isinstance(idx, (int, np.integer)):
        idx = int(idx)
    return mats.iloc[idx, :]


# ------------------------- Main extraction APIs ----------------------------- #

def run_inc_level_analysis_all_samples_SE(
    rmats_output: pd.DataFrame,
    all_eventsDF: pd.DataFrame,
    comparison: str,
    fdr_threshold: float = 0.0,
    inclevel_threshold: float = 0.2
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Look up significant SE events from a comparison in the full rMATS table and
    return inclusion levels across all samples + the full event rows.

    Returns
    -------
    incLevelDF : DataFrame (rows = gene/Event_ID; columns = samples)
    wholeEventsDF : DataFrame (rows aligned to incLevelDF; full rMATS event rows)
    """
    incLevelDict: Dict[str, List[float]] = {}
    wholeEventDict: Dict[str, pd.Series] = {}
    genes_used: List[str] = []

    table_of_events = get_events(rmats_output, fdr_threshold, inclevel_threshold)

    for _, event in table_of_events.iterrows():
        geneName = event['GeneID']
        if isinstance(geneName, datetime.datetime):
            continue  # ignore datetime objects

        if geneName in genes_used:
            # Keep first event per gene, skip duplicates
            # print(f"repeat gene: {geneName}")
            continue

        idx = get_splicing_event_index_SE(all_eventsDF, event)
        if len(idx) == 0:
            # Not found in the big table; skip cleanly
            continue

        inc1, inc2 = get_splicing_event_inc_levels(all_eventsDF, idx)
        inc = inc1 + inc2

        event_id = event.get('ID', 'NA')
        row_name = f"{geneName}--Event_ID: {event_id} ({comparison}:SE)"

        incLevelDict[row_name] = inc
        wholeEventDict[row_name] = get_splicing_event_exon_locations(all_eventsDF, idx)
        genes_used.append(geneName)

    incLevelDF = pd.DataFrame.from_dict(incLevelDict, orient='index')
    wholeEventsDF = pd.DataFrame.from_dict(wholeEventDict, orient='index')

    # Ensure numeric (just in case)
    incLevelDF = incLevelDF.apply(pd.to_numeric, errors='coerce')

    # print('--Finished SE--')
    return incLevelDF, wholeEventsDF


def run_inc_level_analysis_all_samples_MXE(
    rmats_output: pd.DataFrame,
    all_eventsDF: pd.DataFrame,
    comparison: str,
    fdr_threshold: float = 0.0,
    inclevel_threshold: float = 0.2
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Same as SE version, for MXE events.
    """
    incLevelDict: Dict[str, List[float]] = {}
    wholeEventDict: Dict[str, pd.Series] = {}
    genes_used: List[str] = []

    table_of_events = get_events(rmats_output, fdr_threshold, inclevel_threshold)

    for _, event in table_of_events.iterrows():
        geneName = event['GeneID']
        if isinstance(geneName, datetime.datetime):
            continue  # ignore datetime objects

        if geneName in genes_used:
            continue

        idx = get_splicing_event_index_MXE(all_eventsDF, event)
        if len(idx) == 0:
            continue

        inc1, inc2 = get_splicing_event_inc_levels(all_eventsDF, idx)
        inc = inc1 + inc2

        event_id = event.get('ID', 'NA')
        row_name = f"{geneName}--Event_ID: {event_id} ({comparison}:MXE)"

        incLevelDict[row_name] = inc
        wholeEventDict[row_name] = get_splicing_event_exon_locations(all_eventsDF, idx)
        genes_used.append(geneName)

    incLevelDF = pd.DataFrame.from_dict(incLevelDict, orient='index')
    wholeEventsDF = pd.DataFrame.from_dict(wholeEventDict, orient='index')

    incLevelDF = incLevelDF.apply(pd.to_numeric, errors='coerce')

    # print('finished MXE')
    return incLevelDF, wholeEventsDF


# --------------------------- Utilities / QC --------------------------------- #

def inc_level_merge(df1: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
    """
    Merge inclusion level matrices keeping at most one event per gene.
    Right-hand events with genes already present in df1 are discarded.

    Assumes row names follow the pattern: 'GeneName--Event_ID: <id> (<comp>:TYPE)'
    """
    def _gene_from_rowname(idx: str) -> str:
        return idx.split('--Event_ID:')[0] if '--Event_ID:' in idx else idx

    genes_left = {_gene_from_rowname(i) for i in df1.index}
    keep_rows = []
    for i in df2.index:
        gene = _gene_from_rowname(i)
        if gene not in genes_left:
            keep_rows.append(i)

    df2_unique = df2.loc[keep_rows, :]
    out = pd.concat([df1, df2_unique], axis=0)
    return out