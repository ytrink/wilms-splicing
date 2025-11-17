# wilms-splicing

For more information, please see our paper:  Characterization of Alternative Splicing in High-Risk Wilms’ Tumors 

All scripts are fully reproducible and modular, focusing on:

1. Principal component and archetype analysis of bulk RNAseq samples.

2. Deconvolution of bulk tumor expression using single-cell references (CPM method).

3. Visualization of inclusion levels from rMATS splicing analysis in PCA space.

4. ComplexHeatmap plots summarizing sample and gene/splicing event relationships.

## Requirements:
Python 3.10+:
Required packages:

pandas
numpy
matplotlib
scikit-learn

R 4.2+:
Required packages:

ComplexHeatmap
dplyr
scBio
stringr
readxl

Matlab 2019+:
ParTI method from Alon Lab: https://github.com/AlonLabWIS/ParTI

## Repository Structure

```bash
src/
├── figure_1/     # Archetype analysis
├── figure_2/     # CPM deconvolution
├── figure_3/     # rMATS splicing inclusion analysis
└── data/         # Expression matrices and metadata
```



For access to gene expression datasets, splicing results, or any other questions, please email trinkya@biu.ac.il









