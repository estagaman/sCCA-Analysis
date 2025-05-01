# Sparse CCA Analysis

This repository contains code and analysis for replicating and extending the work from Witten and Tibshirani's paper "Extensions of Sparse Canonical Correlation Analysis with Applications to Genomic Data" (2009).
https://pubmed.ncbi.nlm.nih.gov/19572827/

## Overview
The project analyzes gene expression and DNA copy number data from lymphoma patients using various methods:
- Sparse Canonical Correlation Analysis (SCCA) using the PMA package 
- Supervised Sparse CCA (sSCCA) using the PMA package
- Supervised PCA using mixOmics

## Data
We used the GSE11318 dataset and performed the following data collapsing steps:
- For gene expression data: averaged expression values from multiple probes targeting the same gene
- For copy number data: averaged every 10 adjacent copy number locations to reduce dimensionality

### Data Files
- `metadata.csv`: Contains patient metadata including disease subtype and survival information
- `collapsed_expression.csv`: collapsed gene expression data
- `collapsed_CGH.csv`: Collapsed copy number data 

## Key Findings
- Identified significant associations between gene expression and copy number variations
- Found stronger correlations with disease subtype compared to survival status
- Chromosomes 3, 6, 8, 9, 12, 13, and 15 showed significant associations with disease subtype
- Chromosomes 8, 12, 13, and 15 showed significant associations with survival outcome

## Dependencies
- R packages: mixOmics, PMA, GEOquery, biomaRt, tidyverse

