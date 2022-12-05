#!/usr/bin/Rscript
# -*- coding: utf-8 -*

if (!require("pacman")) install.packages("pacman")
pacman::p_load("readr", "feather", "rlist", "data.table")

sc_exp_data = get(load("/home/mraevsky/MIPT/epi-impute/data/GSE117498_scRNAseq_genes.Rdata"))
sc_atac_cell_names_ <- read_csv("/home/mraevsky/MIPT/data/GSE96769_cell_names_matrix.csv")
sc_atac_data <- read_feather("/home/mraevsky/MIPT/data/GSE96769_scATACseq_matrix.feather")
sc_atac_peaks_ann_ <- read_csv("/home/mraevsky/MIPT/data/GSE96769_PeakFile.csv")

#' Compare chromatin accessibility and expression of a particular gene
#'
#' @import rlist
#' @import data.table
#'
#' @export
compare_access_exp_of_gene <- function(sc_exp_data, sc_atac_data,
									   celltypes, gene, gene_regions_in_atac,
									   sc_atac_cell_names = sc_atac_cell_names_, 
									   sc_atac_peaks_ann = sc_atac_peaks_ann_){

	if (!require("rlist")) install.packages("rlist")
	if (!require("data.table")) install.packages("data.table")

	# Convert to matrix for fast transpose
	# Q: Do I need as.data.frame()?
	sc_exp_data = as.data.frame(t(as.matrix(sc_exp_data))) 
	sc_exp_data$gene = rownames(sc_exp_data)

	# Create version of ATAC matrix with only promoter peaks


	gene_regions = sc_atac_peaks_ann[grepl(paste(gene_regions_in_atac, collapse="|"), sc_atac_peaks_ann$X7),]
	gene_regions_vec = gene_regions[order(gene_regions$X13),]

	sc_atac_data = merge(sc_atac_data, gene_regions_vec[c('X13', 'PeakFile_Peak_ID')], by = 'PeakFile_Peak_ID')


	sc_atac_data <- data.table::setDT(sc_atac_data)
	sc_atac_data = sc_atac_data[, lapply(.SD, sum), by = "X13", .SDcols = - "PeakFile_Peak_ID"]
	# sc_atac_data = aggregate(. ~ X13, subset(sc_atac_data, select = - PeakFile_Peak_ID), sum)
	sc_atac_genes = sc_atac_data$X13
	sc_atac_data$X13 = NULL

	# Select only TFs intersected between ATAC and RNA
	genes_common = intersect(sc_atac_genes, sc_exp_data$gene)

	sc_atac_data = sc_atac_data[sc_atac_genes %in% genes_common, ]
	sc_atac_data = sc_atac_data[order(genes_common),]

	sc_exp_data_genes = sc_exp_data[sc_exp_data$gene %in% genes_common, ]
	sc_exp_data_genes = sc_exp_data_genes[order(genes_common),]

	sc_exp_data$gene = NULL

	celltypes_regex = paste(celltypes, collapse = "|")
	selected_atac_cells = sc_atac_cell_names$cell_id[grepl(celltypes_regex, sc_atac_cell_names$cell_types)]
	sc_atac_data = as.data.frame(sc_atac_data)
	rownames(sc_atac_data) = genes_common
	test_atac = sc_atac_data[gene, selected_atac_cells]

	selected_exp_cells = grepl(celltypes_regex, colnames(sc_exp_data))
	test_exp = sc_exp_data[gene, selected_exp_cells]

	par(mfrow = c(1,2))
	hist(as.matrix(test_atac), main = NULL, xlab = sprintf("scATAC-seq counts \n (%s)", paste(gene_regions_in_atac, collapse = ", ")))
	hist(as.matrix(test_exp), main = NULL, xlab = "scRNA-seq counts")
	title(main = sprintf("%s in %s", gene, paste(celltypes, collapse = ", ")), outer = TRUE, line = -1)
}


compare_access_exp_of_gene(sc_exp_data, sc_atac_data, celltypes = "HSC", gene = "CD34", gene_regions_in_atac = c("promoter", "TSS"))