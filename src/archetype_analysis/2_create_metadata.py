# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 15:46:05 2023

@author: trinkya
"""

import pandas as pd

# %% load
counts = pd.read_csv('../data/tcga_counts_all_threeArcs.csv', index_col=0)
gdc_sample_sheet = pd.read_csv('../data/gdc_sample_sheet.2020-12-01.csv')

# %% reorder gdc_sample_sheet so that it is in the same order as the gene expression matrix

fileID = gdc_sample_sheet['File ID']
colnames = list(counts.columns)

# remove X/x from Matlab/Excel exports
colnames = [i[1:] if isinstance(i, str) and i and i[0] in ('x', 'X') else i for i in colnames]  # FIX: handle both cases

# drop archetype columns if they exist
for arc in ('Arc1_Blast', 'Arc2_Stromal', 'Arc3_Normal'):  # FIX: safe removal
    if arc in colnames:
        colnames.remove(arc)

# TCGA sample IDs typically use '-', not '_'
colnames = [i.replace('_', '-') if isinstance(i, str) else i for i in colnames]

gdc_sample_sheet.set_index('File ID', inplace=True)

new_sample_sheet = gdc_sample_sheet.reindex(colnames)

print(set(colnames).difference(set(fileID)))

new_sample_sheet.to_csv('../data/gdc_sample_sheet.2020-12-01_reordered.csv')



metadata = pd.read_json('../data/metadata.cart.2020-12-01.json')


a = pd.DataFrame(metadata['associated_entities'][0][0], index=[0])

associated_entities = metadata['associated_entities']
for i in range(1, metadata.shape[0]):
    # defensive: if list empty, use {}
    ent = associated_entities[i][0] if (isinstance(associated_entities[i], list) and len(associated_entities[i]) > 0) else {}
    newCols = pd.DataFrame(ent, index=[i])
    a = pd.concat([a, newCols])

metadata.drop('associated_entities', inplace=True, axis=1)
metadata = metadata.join(a)

a = pd.DataFrame(metadata['index_files'][0], index=[0])
for i in range(1, metadata.shape[0]):
    idx_entry = metadata['index_files'][i] if isinstance(metadata['index_files'][i], dict) else {}
    newCols = pd.DataFrame(idx_entry, index=[i])
    a = pd.concat([a, newCols])

metadata.drop('index_files', inplace=True, axis=1)

metadata = pd.merge(metadata, a, left_index=True, right_index=True)
metadata.to_csv('../data/metadata.cart.split.csv')

# %% reorder TARGET metadata.cart
# link is through "file_name" of the bam file

gdc_sample_sheet = pd.read_csv('../data/gdc_sample_sheet.2020-12-01_reordered.csv', index_col=0)
file_names = gdc_sample_sheet['File Name']

metadata.set_index('file_name', inplace=True)

# check if there are any differences between the file names
print(set(metadata.index).difference(set(file_names)))

metadata = metadata.reindex(file_names)
gdc_sample_sheet.set_index('File Name', inplace=True)

metadata2 = pd.merge(gdc_sample_sheet, metadata, left_index=True, right_index=True).reset_index()

# %% merge metadata.cart/gdc_sample_sheet with TARGET_WT_ClinicalData_Discovery using the Case ID

clinical = pd.read_excel('../data/Copy of TARGET_WT_ClinicalData_Discovery_20211111-1.xlsx', index_col=0)
caseID = metadata2['Case ID']

# build clinical data row-aligned; handle missing Case IDs safely  # FIX
rows = []
for i in range(len(caseID)):
    case = caseID.iloc[i]
    if pd.notna(case) and case in clinical.index:
        rows.append(clinical.loc[case, :])
    else:
        rows.append(pd.Series(index=clinical.columns, dtype=object))
clinicalDF = pd.DataFrame(rows).reset_index(drop=True)

metadataFinal = pd.merge(metadata2.reset_index(drop=True),
                         clinicalDF.reset_index(drop=True),
                         left_index=True, right_index=True)

metadataFinal.to_csv('../data/metadata_merged_for_plots.csv')

# %% combine deconvolution and distance from archetypes

metadata_merged_for_plots = pd.read_csv('../data/metadata_merged_for_plots.csv', index_col='file_id_x')

# distance from archetypes
import sys
sys.path.append('../figure_1')

from _wilms_nD import wilms_Narcs
from sklearn.metrics.pairwise import euclidean_distances

# %% first merge with euclidean distance matrix --

wilms = wilms_Narcs('../data/tcga_counts_all_threeArcs.csv')

samples = wilms.get_column_names()
X = wilms.get_score()

dist = euclidean_distances(X, X)
distDF = pd.DataFrame(dist)
distDF.index = samples

# rename columns to "dist_from"
distDF.columns = ['dist_from_' + i for i in samples]

# create rows for archetypes
x = metadata_merged_for_plots.shape[1] - 1

blastRow = ['Arc1_Blast'] + ['Arc'] * x
stromalRow = ['Arc2_Stromal'] + ['Arc'] * x
NormalRow = ['Arc3_Normal'] + ['Arc'] * x

arcsFrame = pd.DataFrame([blastRow, stromalRow, NormalRow])
arcsFrame.columns = metadata_merged_for_plots.columns
arcsFrame.index = ['Arc1_Blast', 'Arc2_Stromal', 'Arc3_Normal']

metadataAll = pd.concat([metadata_merged_for_plots, arcsFrame], axis=0)

# use only distances to archetypes (if present)
distDF2 = distDF[['dist_from_Arc1_Blast', 'dist_from_Arc2_Stromal', 'dist_from_Arc3_Normal']].copy()

# align indices by your underscore convention WITHOUT overwriting metadata indexes  # FIX
distDF2.index = [s.replace('-', '_') for s in distDF2.index]  # make distance index like your metadata rows
metadataAll = metadataAll.join(distDF2, how='left')           # safe join by index (no destructive reindex)

# %% now merge deconvolution

deconvolution_output = pd.read_csv('../data/deconvolution_output/cpm_cell_type_predictions.csv', index_col=0)
sample_names = list(deconvolution_output.index)
sample_names = [i[1:] if isinstance(i, str) and i and i[0] in ('X', 'x') else i for i in sample_names]  # FIX
deconvolution_output.index = sample_names

metadataAll = pd.concat([metadataAll, deconvolution_output], axis=1)

# %% insert normal tags into columns

# avoid chained assignment; write into a new column using .loc  # FIX
if 'Histologic Classification of Primary Tumor' in metadataAll.columns and 'Sample Type' in metadataAll.columns:
    new_col = 'Histologic Classification of Primary Tumor Normals Labelled'
    metadataAll[new_col] = metadataAll['Histologic Classification of Primary Tumor']
    idx = metadataAll['Sample Type'].eq('Solid Tissue Normal')
    metadataAll.loc[idx, new_col] = 'Solid Tissue Normal'

metadataAll.to_csv('../data/metadataAll.csv')  # use this for analysis
