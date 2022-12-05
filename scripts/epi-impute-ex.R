#!/usr/bin/Rscript
# -*- coding: utf-8 -*

if (!require("pacman")) install.packages("pacman")
pacman::p_load("readr", "rlist", "reshape2", "feather")

sc_atac_cell_names <- read_csv("/home/mraevsky/MIPT/data/GSE96769_cell_names_matrix.csv")
sc_atac_data <- read_feather("/home/mraevsky/MIPT/data/GSE96769_scATACseq_matrix.feather")
sc_atac_peaks_ann <- read_csv("/home/mraevsky/MIPT/data/GSE96769_PeakFile.csv")
sc_exp_data = get(load("/home/mraevsky/MIPT/epi-impute/data/GSE117498_scRNAseq_genes.Rdata"))


# Import scATACseq data and metadata

sc_atac_cell_names_ <- read_csv("/home/mraevsky/MIPT/data/GSE96769_cell_names_matrix.csv")
sc_atac_data_ <- read_feather("/home/mraevsky/MIPT/data/GSE96769_scATACseq_matrix.feather")
sc_atac_peaks_ann_ <- read_csv("/home/mraevsky/MIPT/data/GSE96769_PeakFile.csv")


# Load scRNAseq dataset that will be imputed

# sc_exp_data = read_feather('scRNA_ATAC_task/data/GSE117498_scRNAseq_dataset.feather')
# rownames(sc_exp_data) = sc_exp_data$gene
# sc_exp_data$gene = NULL

epi_impute <- function(sc_exp_data, atac_bin_thrld = 100, 
									sc_atac_data = sc_atac_data_, 
									sc_atac_cell_names = sc_atac_cell_names_, 
									sc_atac_peaks_ann = sc_atac_peaks_ann_){
	# sc_exp_data should be a cells x genes matrix 
	if (!require("rlist")) install.packages("rlist")
	if (!require("data.table")) install.packages("data.table")

	# Convert to matrix for fast transpose
	# Q: Do I need as.data.frame()?
	sc_exp_data = as.data.frame(t(as.matrix(sc_exp_data))) 
	sc_exp_data$gene = rownames(sc_exp_data)

	# Create version of ATAC matrix with only promoter peaks

	promoters = sc_atac_peaks_ann[grepl('promoter|TSS', sc_atac_peaks_ann$X7),]
	promoters_genes = promoters[order(promoters$X13),]

	sc_atac_data = merge(sc_atac_data, promoters_genes[c('X13', 'PeakFile_Peak_ID')], by = 'PeakFile_Peak_ID')


	sc_atac_data <- setDT(sc_atac_data)
	sc_atac_data = sc_atac_data[, lapply(.SD, sum), by = "X13", .SDcols = - "PeakFile_Peak_ID"]
	# sc_atac_data = aggregate(. ~ X13, subset(sc_atac_data, select = - PeakFile_Peak_ID), sum)
	sc_atac_genes = sc_atac_data$X13
	sc_atac_data$X13 = NULL

	dump = sc_atac_data

	# Select only TFs intersected between ATAC and RNA
	genes_common = intersect(sc_atac_genes, sc_exp_data$gene)

	sc_atac_data = sc_atac_data[sc_atac_genes %in% genes_common, ]
	sc_atac_data = sc_atac_data[order(genes_common),]

	sc_exp_data_genes = sc_exp_data[sc_exp_data$gene %in% genes_common, ]
	sc_exp_data_genes = sc_exp_data_genes[order(genes_common),]

	# Filter ATAC matrix from unused celltypes and cells

	sc_atac_cell_names = sc_atac_cell_names[grepl('HSC|CMP|GMP|MLP', sc_atac_cell_names$cell_types),]
	sc_atac_data = as.data.frame(sc_atac_data) # next line will generate bool if comment this because of data.frame or use ", with=FALSE"
	sc_atac_data = sc_atac_data[,colnames(sc_atac_data) %in% sc_atac_cell_names$cell_id]
	sc_atac_cell_names = sc_atac_cell_names[sc_atac_cell_names$cell_id %in% colnames(sc_atac_data),]

	# Sum peaks counts relative to celltypes
	 
	celltype_list = split(sc_atac_cell_names, sc_atac_cell_names$cell_types)
	celltype_list = celltype_list[c("HSC", "CMP", "GMP")] # there is no MLP

	bulk_peaks_list = list()
	for(celltype in celltype_list){
		bulk_peaks_df = rowSums(sc_atac_data[, celltype$cell_id])
		# threshold = (dim(celltype_list)[1] * 0.25)
		bulk_peaks_df = as.data.frame(as.matrix((bulk_peaks_df > atac_bin_thrld) + 0)) # binarize it with t
		bulk_peaks_list = list.append(bulk_peaks_list, bulk_peaks_df)
	}
	names(bulk_peaks_list) = c("HSC", "CMP", "GMP")

	# Split scRNAseq by celltype

	print('Binirize the matrix')
	# Binirize the matrix
	sc_exp_data_genes_values = subset(sc_exp_data_genes, select = - gene)
	sc_exp_data_genes_values = as.data.frame(as.matrix((sc_exp_data_genes_values > 0) + 0))
	sc_exp_data_genes = cbind(sc_exp_data_genes$gene, sc_exp_data_genes_values)
	colnames(sc_exp_data_genes)[1] = "gene"

	sc_exp_data_genes = dcast(melt(sc_exp_data_genes, id.vars = "gene"), variable ~ gene)
	rownames(sc_exp_data_genes) = sc_exp_data_genes$variable
	sc_exp_data_genes$cells = gsub('_.*', '', sc_exp_data_genes$variable)
	sc_exp_data_genes = sc_exp_data_genes[order(sc_exp_data_genes$cells),]
	sc_exp_data_genes$variable = NULL

	sc_data_list <- split(sc_exp_data_genes , f = sc_exp_data_genes$cells)
	sc_data_list = lapply(sc_data_list, function(x){x$cells = NULL; as.data.frame(t(x))})

	print('Imputation...')
	runtime = system.time({
		## Imputation
		celltypes_to_impute = c('HSC', 'CMP', 'GMP')
		scRNAseq_imputed_list = list()

		for(i in celltypes_to_impute){
			imputed_i = as.data.frame(sapply(as.data.frame(sc_data_list[[i]]), function(x){x + bulk_peaks_list[[i]]}))
			scRNAseq_imputed_list = list.append(scRNAseq_imputed_list, as.data.frame(imputed_i))
		}

		scRNAseq_imputed_list = as.data.frame(scRNAseq_imputed_list)
		# scRNAseq_imputed_list$gene = genes_common
		rownames(scRNAseq_imputed_list) = genes_common
	})

	return(list(result = as.data.frame(t(scRNAseq_imputed_list)),
				runtime = runtime))
}



# # Export processed scRNAseq and scATACseq matrixes for other methods

# sc_exp_data = as.data.frame(sc_data_list[c("HSC", "CMP", "GMP")])
# sc_exp_data = as.data.frame(t(sc_exp_data))
# rownames(sc_exp_data) = sub('[A-Z]{3}\\.', '', rownames(sc_exp_data))

# save(sc_exp_data, file = "~/MIPT/epi_impute/data/GSE117498_scRNAseq_genes.Rdata")
# # write_feather(sc_exp_data, "~/MIPT/epi_impute/data/GSE117498_scRNAseq_genes.feather")


# bulk_peaks_data = as.data.frame(bulk_peaks_list[c("HSC", "CMP", "GMP")])
# bulk_peaks_data = as.data.frame(t(bulk_peaks_data))
# rownames(bulk_peaks_data) = c("HSC", "CMP", "GMP")
# colnames(bulk_peaks_data) = colnames(sc_exp_data)

# save(bulk_peaks_data, file = "~/MIPT/epi_impute/data/bulk_ATAC_genes.Rdata")
# # write_feather(bulk_peaks_data, "~/MIPT/epi_impute/data/bulk_ATAC_genes.feather")