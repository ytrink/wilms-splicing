## Single Cell Deconvolution (CPM) on count matrix with three archetypes
# Using human kidney cell atlas
# --------------------------------------------------------------
# Performs single-cell deconvolution for kidney Wilms tumor bulk expression.
# The function requires four inputs:
# 1) Bulk expression set
# 2) Single-cell expression set
# 3) Labels for each single cell (cell type)
# 4) Embedded single-cell coordinates (PCA/TSNE/UMAP)
# --------------------------------------------------------------

# Install dependencies if missing


# CPM algorithm by Frishberg et al - see details 
# https://pubmed.ncbi.nlm.nih.gov/30886410/



if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if (!requireNamespace("scBio", quietly = TRUE)) install.packages("scBio")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")

# Load libraries
library(stringr)
library(readxl)
library(scBio)
library(limma)

# ---------------------- Load Data ----------------------

# Bulk expression data
bulkExpression <- read.csv('../data/tcga_counts_all_threeArcs.csv', row.names = 1)

# Single-cell counts (from Weinberg et al.)
load('../deconvolution_files/Fetal_full_v3_scnormalized_counts.R')  # contains 'counts'

# Selected cell indices
idx.cells <- read.csv('../deconvolution_files/index_selected_cells.csv')$x

# Cell-type labels
SCLabels.selected.grouped <- read.csv('../deconvolution_files/cell_labels_selected_only_grouped.csv', row.names = 1)
SCLabels.selected.grouped <- SCLabels.selected.grouped$X0

# Prepare single-cell data
SCData <- counts
rm(counts)
SCData.selected <- SCData[, idx.cells]
rownames(SCData) <- toupper(rownames(SCData))

# ---------------------- Filter Genes ----------------------


# 728 relevant genes to reduce noise
genes.728 <- read.csv('../deconvolution_files/gene_names4tSNE (728 genes).csv')
genes.728 <- toupper(genes.728[[1]])

commonGenes <- intersect(rownames(bulkExpression), rownames(SCData))
commonGenes <- intersect(commonGenes, genes.728)

SCData.selected <- SCData.selected[commonGenes, ]
bulkExpression <- bulkExpression[commonGenes, ]

# ---------------------- Embedding ----------------------

cellSpace <- read.csv('../deconvolution_files/fetal_kidney_atlas_umap_embedding.csv', row.names = 1)
cellSpace <- cellSpace[idx.cells, ]

cellSpaceArray <- as.matrix(cellSpace[, 1:2])
rownames(cellSpaceArray) <- rownames(cellSpace)

# ---------------------- CPM Deconvolution ----------------------

# Convert bulk back to absolute scale if log-transformed
bulkExpression.abs <- exp(bulkExpression) - 1

SCData.selected <- as.data.frame(as.matrix(SCData.selected))

set.seed(1)
resAbs.removed <- CPM(
  SCData.selected,
  SCLabels.selected.grouped,
  bulkExpression.abs,
  cellSpaceArray,
  typeTransformation = TRUE,
  quantifyTypes = TRUE
)

# ---------------------- Save Outputs ----------------------

if (!dir.exists('./deconvolution_output')) dir.create('./deconvolution_output')

write.csv(resAbs.removed$predicted, './deconvolution_output/cpm_predicted.csv')
write.csv(resAbs.removed$cellTypePredictions, './deconvolution_output/cpm_cell_type_predictions.csv')

message('CPM deconvolution completed successfully.')
