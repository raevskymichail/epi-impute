# Epi-Impute

[![](https://img.shields.io/github/languages/code-size/raevskymichail/epi-impute)](https://img.shields.io/github/languages/code-size/raevskymichail/epi-impute)
[![](https://img.shields.io/github/languages/top/raevskymichail/epi-impute)](https://img.shields.io/github/languages/top/raevskymichail/epi-impute)
[![](https://img.shields.io/github/issues/raevskymichail/epi-impute)](https://img.shields.io/github/issues/raevskymichail/epi-impute)
[![](https://img.shields.io/github/license/raevskymichail/epi-impute)](https://img.shields.io/github/license/raevskymichail/epi-impute)

## Introduction

This repository contains primary source code for *"Epi-Impute: single-cell RNA-seq imputation via integration with single-cell ATAC-seq"* that is a computational tool for imputing scRNA-seq data from DNA accessibility data (scATAC-seq) from consistent cell-type populations.

Epi-Impute exploits the idea of open chromatin in active *cis*-regulatory elements of the genes and estimates average accessibility of gene regulatory elements, e.g., promoters and enhancers, in order to add a binarized pseudo-count to the gene expression values reflecting the potential for transcription activation observed at the epigenetic level

## Preprocessing

For preprocessing of scRNA-seq data, please follow the standard processing pipeline to get the expression count matrix, where each row represents a gene and each column represents a cell.

For preprocessing of scATAC-seq data, please first put all the `.bam` files for each cell into a dedicated folder. Then run the preprocessing script we provided to get a count matrix and annotations files.

## Requirements

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
install.packages(path_to_archive, type = 'source', rep = NULL)
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

## Help

Please feel free to contact Mikhail Raevskiy (raevskii.mm@phystech.edu) if you have any questions about the software.
