library(readr)
library(rlist)
library(reshape2)
library(feather)

setwd('~/MIPT/')

sc_atac_cell_names <- read_csv("data/GSE96769_cell_names_matrix.csv")
sc_atac_data <- read_feather("data/GSE96769_scATACseq_matrix.feather")
sc_atac_peaks_ann <- read_csv("data/GSE96769_PeakFile.csv")
knonwn_tfs <- read_csv("scRNA_ATAC_task/INTREGNET_Reproducible/AnimalTFDB_TF.bed", col_names = FALSE)
knonwn_tfs = knonwn_tfs$X1
sc_exp_data = get(load("~/MIPT/epi_impute/data/GSE117498_scRNAseq_TFs.Rdata"))


# Import scATACseq data and metadata

sc_atac_cell_names_ <- read_csv("data/GSE96769_cell_names_matrix.csv")
sc_atac_data_ <- read_feather("data/GSE96769_scATACseq_matrix.feather")
sc_atac_peaks_ann_ <- read_csv("data/GSE96769_PeakFile.csv")
knonwn_tfs_ <- read_csv("scRNA_ATAC_task/INTREGNET_Reproducible/AnimalTFDB_TF.bed", col_names = FALSE)
knonwn_tfs_ = knonwn_tfs_$X1

# Load scRNAseq dataset that will be imputed

# sc_exp_data = read_feather('scRNA_ATAC_task/data/GSE117498_scRNAseq_dataset.feather')
# rownames(sc_exp_data) = sc_exp_data$gene
# sc_exp_data$gene = NULL

epi_impute <- function(sc_exp_data, atac_bin_thrld, 
									sc_atac_data = sc_atac_data_, 
									sc_atac_cell_names = sc_atac_cell_names_, 
									sc_atac_peaks_ann = sc_atac_peaks_ann_, 
									knonwn_tfs = knonwn_tfs_){
	# sc_exp_data should be a cells x genes matrix 

	sc_exp_data = as.data.frame(t(sc_exp_data))
	sc_exp_data$gene = rownames(sc_exp_data)

	# Create version of ATAC matrix with only promoter peaks

	promoters = sc_atac_peaks_ann[grepl('promoter|TSS', sc_atac_peaks_ann$X7),]
	promoters_tfs = promoters[promoters$X13 %in% knonwn_tfs,]
	promoters_tfs = promoters_tfs[order(promoters_tfs$X13),]

	sc_atac_data = merge(sc_atac_data, promoters_tfs[c('X13', 'PeakFile_Peak_ID')], by = 'PeakFile_Peak_ID')

	sc_atac_data = aggregate(. ~ X13, subset(sc_atac_data, select = - PeakFile_Peak_ID), sum)
	sc_atac_tfs = sc_atac_data$X13
	sc_atac_data$X13 = NULL

	# Select only TFs intersected between ATAC and RNA
	tfs_common = intersect(sc_atac_tfs, sc_exp_data$gene)

	sc_atac_data = sc_atac_data[sc_atac_tfs %in% tfs_common, ]
	sc_atac_data = sc_atac_data[order(tfs_common),]

	sc_exp_data_TFs = sc_exp_data[sc_exp_data$gene %in% tfs_common, ]
	sc_exp_data_TFs = sc_exp_data_TFs[order(tfs_common),]

	# Filter ATAC matrix from unused celltypes and cells

	sc_atac_cell_names = sc_atac_cell_names[grepl('HSC|CMP|GMP|MLP', sc_atac_cell_names$cell_types),]
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
	sc_exp_data_TFs_values = subset(sc_exp_data_TFs, select = - gene)
	sc_exp_data_TFs_values = as.data.frame(as.matrix((sc_exp_data_TFs_values > 0) + 0))
	sc_exp_data_TFs = cbind(sc_exp_data_TFs$gene, sc_exp_data_TFs_values)
	colnames(sc_exp_data_TFs)[1] = "gene"

	sc_exp_data_TFs = dcast(melt(sc_exp_data_TFs, id.vars = "gene"), variable ~ gene)
	rownames(sc_exp_data_TFs) = sc_exp_data_TFs$variable
	sc_exp_data_TFs$cells = gsub('_.*', '', sc_exp_data_TFs$variable)
	sc_exp_data_TFs = sc_exp_data_TFs[order(sc_exp_data_TFs$cells),]
	sc_exp_data_TFs$variable = NULL

	sc_data_list <- split(sc_exp_data_TFs , f = sc_exp_data_TFs$cells)
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
		# scRNAseq_imputed_list$gene = tfs_common
		rownames(scRNAseq_imputed_list) = tfs_common
	})

	return(list(result = as.data.frame(t(scRNAseq_imputed_list)),
				runtime = runtime))
}



# # Export processed scRNAseq and scATACseq matrixes for other methods

# sc_exp_data = as.data.frame(sc_data_list[c("HSC", "CMP", "GMP")])
# sc_exp_data = as.data.frame(t(sc_exp_data))
# rownames(sc_exp_data) = sub('[A-Z]{3}\\.', '', rownames(sc_exp_data))

# save(sc_exp_data, file = "~/MIPT/epi_impute/data/GSE117498_scRNAseq_TFs.Rdata")
# # write_feather(sc_exp_data, "~/MIPT/epi_impute/data/GSE117498_scRNAseq_TFs.feather")


# bulk_peaks_data = as.data.frame(bulk_peaks_list[c("HSC", "CMP", "GMP")])
# bulk_peaks_data = as.data.frame(t(bulk_peaks_data))
# rownames(bulk_peaks_data) = c("HSC", "CMP", "GMP")
# colnames(bulk_peaks_data) = colnames(sc_exp_data)

# save(bulk_peaks_data, file = "~/MIPT/epi_impute/data/bulk_ATAC_TFs.Rdata")
# # write_feather(bulk_peaks_data, "~/MIPT/epi_impute/data/bulk_ATAC_TFs.feather")