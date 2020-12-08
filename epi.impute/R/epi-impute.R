#!/usr/bin/env Rscript

#' Epi-Impute: single-cell RNA-seq imputation via integration with single-cell ATAC-seq
#' 
#' Epi-Impute, a computational tool for imputing scRNA-seq data from DNA accessibility data (scATAC-seq) from consistent cell-type populations.
#' @param sc_exp_data a single-cell RNA-seq count matrix, where rownames are
#' HGNC genes (HUGO) and colnames are cell ids
#' @param sc_atac_data a single-cell ATAC-seq count matrix, where rownames are
#' cell ids and colnames are ids for euchromatine peaks (obtained from peak caller,
#' for ex. MACS2)
#' @param sc_atac_cell_names a matrix, containing description and annotation for
#' cell types observed in scATAC-seq count matrix. It should have rownames (cell
#' ids) that
#' match rownames of \code{sc_atac_data}
#' @param sc_atac_peaks_ann a matrix, containing description for euchromatine
#' peaks, presented in count matrix. It should have rownames (peak ids) that
#' match colnames of \code{sc_atac_data}
#' @param cell_types a vector, containing names for cell types, presented in
#' count matrix.
#' @param atac_bin_thrld a numeric value for accessibility threshold used for
#' primary binirization of peaks in scATAC-seq matrix.
#' @return a data frame containg an imputed single-cell RNA-seq matrix
#' 
#' @export
epi_impute <- function(sc_exp_data, sc_atac_data, sc_atac_cell_names,
                                    sc_atac_peaks_ann, cell_types,
                                    atac_bin_thrld = 100, ...) {
  sc_exp_data <- as.data.frame(t(as.matrix(sc_exp_data)))
  sc_exp_data$gene <- rownames(sc_exp_data)

  promoters <- sc_atac_peaks_ann[grepl('promoter|TSS', sc_atac_peaks_ann$X7),]
  promoters_genes <- promoters[order(promoters$X13),]

  sc_atac_data <- merge(sc_atac_data, promoters_genes[c('X13', 'PeakFile_Peak_ID')], by = 'PeakFile_Peak_ID')

  sc_atac_data <- setDT(sc_atac_data)
  sc_atac_data <- sc_atac_data[, lapply(.SD, sum), by = "X13", .SDcols = -"PeakFile_Peak_ID"]
  sc_atac_genes <- sc_atac_data$X13
  sc_atac_data$X13 <- NULL

  genes_common <- intersect(sc_atac_genes, sc_exp_data$gene)
  sc_atac_data <- sc_atac_data[sc_atac_genes %in% genes_common,]
  sc_atac_data <- sc_atac_data[order(genes_common),]

  sc_exp_data_genes <- sc_exp_data[sc_exp_data$gene %in% genes_common,]
  sc_exp_data_genes <- sc_exp_data_genes[order(genes_common),]

  sc_atac_cell_names <- sc_atac_cell_names[grepl(paste(cell_types, collapse = "|"), sc_atac_cell_names$cell_types),]
  sc_atac_data <- as.data.frame(sc_atac_data) # next line will generate bool if comment this because of data.frame or use ", with=FALSE"
  sc_atac_data <- sc_atac_data[, colnames(sc_atac_data) %in% sc_atac_cell_names$cell_id]
  sc_atac_cell_names <- sc_atac_cell_names[sc_atac_cell_names$cell_id %in% colnames(sc_atac_data),]

  celltype_list <- split(sc_atac_cell_names, sc_atac_cell_names$cell_types)
  celltype_list <- celltype_list[cell_types] # there is no MLP

  bulk_peaks_list <- list()
  for (celltype in celltype_list) {
    bulk_peaks_df <- rowSums(sc_atac_data[, celltype$cell_id])
    bulk_peaks_df <- as.data.frame(as.matrix((bulk_peaks_df > atac_bin_thrld) + 0))
    bulk_peaks_list <- list.append(bulk_peaks_list, bulk_peaks_df)
  }
  names(bulk_peaks_list) <- cell_types

  print('Binirize the matrix')
  sc_exp_data_genes_values <- subset(sc_exp_data_genes, select = -gene)
  sc_exp_data_genes_values <- as.data.frame(as.matrix((sc_exp_data_genes_values > 0) + 0))
  sc_exp_data_genes <- cbind(sc_exp_data_genes$gene, sc_exp_data_genes_values)
  colnames(sc_exp_data_genes)[1] <- "gene"

  sc_exp_data_genes <- dcast(melt(sc_exp_data_genes, id.vars = "gene"), variable ~ gene)
  rownames(sc_exp_data_genes) <- sc_exp_data_genes$variable
  sc_exp_data_genes$cells <- gsub('_.*', '', sc_exp_data_genes$variable)
  sc_exp_data_genes <- sc_exp_data_genes[order(sc_exp_data_genes$cells),]
  sc_exp_data_genes$variable <- NULL

  sc_data_list <- split(sc_exp_data_genes, f = sc_exp_data_genes$cells)
  sc_data_list <- lapply(sc_data_list, function(x) { x$cells = NULL; as.data.frame(t(x)) })

  print('Imputation...')
  runtime <- system.time({
    celltypes_to_impute <- cell_types
    scRNAseq_imputed_list <- list()

    for (i in celltypes_to_impute) {
      imputed_i <- as.data.frame(sapply(as.data.frame(sc_data_list[[i]]), function(x) { x + bulk_peaks_list[[i]] }))
      scRNAseq_imputed_list <- list.append(scRNAseq_imputed_list, as.data.frame(imputed_i))
    }

    scRNAseq_imputed_list <- as.data.frame(scRNAseq_imputed_list)
    rownames(scRNAseq_imputed_list) <- genes_common
  })

  return(as.data.frame(t(scRNAseq_imputed_list)))
}

#' Loads example data for imputation
#'
#' @export
load_example_data <- function() {
  library(RCurl)
  library(readr)
  library(feather)

  sc_atac_cell_names <- read_csv(url("http://himorna.fbras.ru/~mraevsky/Epi-Impute/GSE96769_cell_names_matrix.csv"))
  sc_atac_peaks_ann <- read_csv(url("http://himorna.fbras.ru/~mraevsky/Epi-Impute/GSE96769_PeakFile.csv"))
  sc_exp_data <- get(load(url("http://himorna.fbras.ru/~mraevsky/Epi-Impute/GSE117498_scRNAseq_genes.Rdata")))
  sc_atac_data <- bdown("http://himorna.fbras.ru/~mraevsky/Epi-Impute/GSE96769_scATACseq_matrix.feather", "GSE96769_scATACseq_matrix.feather")

  return(list(sc_atac_cell_names = sc_atac_cell_names,
              sc_atac_data = sc_atac_data,
              sc_atac_peaks_ann = sc_atac_peaks_ann,
              sc_exp_data = sc_exp_data))
}

#' Download file from remote FTP server with specified name
#'
#' @export
bdown <- function(url, file) {
  library("RCurl")
  f <- CFILE(file, mode = "wb")
  a <- curlPerform(url = url, writedata = f@ref, noprogress = FALSE)
  close(f)
  return(a)
}