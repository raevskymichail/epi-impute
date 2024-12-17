# Epi-Impute

[![](https://img.shields.io/github/languages/code-size/raevskymichail/epi-impute)](https://img.shields.io/github/languages/code-size/raevskymichail/epi-impute)
[![](https://img.shields.io/github/languages/top/raevskymichail/epi-impute)](https://img.shields.io/github/languages/top/raevskymichail/epi-impute)
[![](https://img.shields.io/github/issues/raevskymichail/epi-impute)](https://img.shields.io/github/issues/raevskymichail/epi-impute)
[![](https://img.shields.io/github/license/raevskymichail/epi-impute)](https://img.shields.io/github/license/raevskymichail/epi-impute)

## Introduction

This repository contains primary source code for *"Epi-Impute: single-cell RNA-seq imputation via integration with single-cell ATAC-seq"* that is a computational tool for imputing scRNA-seq data from DNA accessibility data (scATAC-seq) from consistent cell-type populations.

Epi-Impute exploits the idea of open chromatin in active *cis*-regulatory elements of the genes and estimates average accessibility of gene regulatory elements, e.g., promoters and enhancers, in order to add a binarized pseudo-count to the gene expression values reflecting the potential for transcription activation observed at the epigenetic level

## Preprocessing

For preprocessing of scRNA-seq data, please follow the standard processing pipeline to get the expression count matrix, where each row represents a gene and each column represents a cell. Epi-Impute supports both raw and normalized data.

For scATAC-seq data, please, obtain a count matrix and annotations with preprocessing pipeline you are using.

## Requirements

Epi-Impute requires R version 4.0.0 or above and following packages:

* readr
* reshape2
* feather
* rlist
* data.table
* RCurl
* MACS
* UMAP (https://github.com/lmcinnes/umap)

## Installation

### Install from Github
```r
library(devtools)
install_github("raevskymichail/epi-impute/epi.impute")
```

### Install from source codes

Download source codes and then type in R session:

```r
install.packages(path_to_archive, type = "source", rep = NULL)
```

Where `path_to_archive` would represent the full path and file name:
- On Windows it will look something like this: `C:\\Downloads\epi-impute.tar.gz`.
- On UNIX machines it will look like this: `~/Downloads/epi-impute.tar.gz`.

## Quick start

```r
library("epi.impute")

data <- load_example_data()

data_imputed <- epi_impute(sc_exp_data = data[["sc_exp_data"]],
                           sc_atac_data = data[["sc_atac_data"]],
                           sc_atac_cell_names = data[["sc_atac_cell_names"]],
                           sc_atac_peaks_ann = data[["sc_atac_peaks_ann"]],
                           cell_types = c("HSC", "CMP", "GMP"))
```

Where

**sc_exp_data** – scRNA-seq count matrix, where rownames are HGNC genes (HUGO) and colnames are *cell ids*

**sc_atac_data** – scATAC-seq count matrix, where rownames are *cell ids* and colnames are ids for euchromatine peaks (obtained from peak caller, for ex. MACS2)

**sc_atac_cell_names** – matrix, containing description and annotation for cell types observed in scATAC-seq count matrix. It should have rownames (*cell ids*) that match rownames of `sc_atac_data`

**sc_atac_peaks_ann** – matrix with annotations (coordinates) for euchromatine peaks, presented in the scATAC-seq count matrix. It should have rownames (*peak ids*) that match colnames of `sc_atac_data`

**cell_types** – vector, containing names for cell types, presented in count matrix.

**atac_bin_thrld** – numeric value for accessibility threshold used for primary binarization of peaks in scATAC-seq matrix.

## Help

Please feel free to contact Mikhail Raevskiy (raevskii.mm@phystech.edu) if you have any questions about the software.
