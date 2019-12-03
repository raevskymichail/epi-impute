library(readr)
library(rlist)
library(reshape2)
library(feather)

setwd('~/MIPT/')

# Import scATACseq data and metadata

GSE96769_cell_names_matrix <- read_csv("data/GSE96769_cell_names_matrix.csv")
GSE96769_scATACseq_matrix <- read_feather("data/GSE96769_scATACseq_matrix.feather")
GSE96769_PeakFile <- read_csv("data/GSE96769_PeakFile.csv")


# Load scRNAseq dataset that will be imputed

GSE117498_scRNAseq_dataset = read_feather('scRNA_ATAC_task/data/GSE117498_scRNAseq_dataset.feather')

# Create version of ATAC matrix with only promoter peaks

GSE96769_PeakFile_promoters = GSE96769_PeakFile[grepl('promoter|TSS', GSE96769_PeakFile$X7),]
GSE96769_PeakFile_promoters_genes = GSE96769_PeakFile_promoters[order(GSE96769_PeakFile_promoters$X13),]

GSE96769_scATACseq_matrix = merge(GSE96769_scATACseq_matrix, GSE96769_PeakFile_promoters_genes[c('X13', 'PeakFile_Peak_ID')], by = 'PeakFile_Peak_ID')

	# aggregate(subset(GSE96769_scATACseq_matrix, select = - PeakFile_Peak_ID), by=list(TF=GSE96769_scATACseq_matrix$X13), FUN=sum)
GSE96769_scATACseq_matrix = aggregate(. ~ X13, subset(GSE96769_scATACseq_matrix, select = - PeakFile_Peak_ID), sum)
GSE96769_genes = GSE96769_scATACseq_matrix$X13
GSE96769_scATACseq_matrix$X13 = NULL

# Select only TFs intersected between ATAC and RNA
genes_common = intersect(GSE96769_genes, GSE117498_scRNAseq_dataset$gene)

GSE96769_scATACseq_matrix = GSE96769_scATACseq_matrix[GSE96769_genes %in% genes_common, ]
GSE96769_scATACseq_matrix = GSE96769_scATACseq_matrix[order(genes_common),]

GSE117498_scRNAseq_dataset_genes = GSE117498_scRNAseq_dataset[GSE117498_scRNAseq_dataset$gene %in% genes_common, ]
GSE117498_scRNAseq_dataset_genes = GSE117498_scRNAseq_dataset_genes[order(genes_common),]

# Filter ATAC matrix from unused celltypes and cells

GSE96769_cell_names_matrix = GSE96769_cell_names_matrix[grepl('HSC|CMP|GMP|MLP', GSE96769_cell_names_matrix$cell_types),]
GSE96769_scATACseq_matrix = GSE96769_scATACseq_matrix[,colnames(GSE96769_scATACseq_matrix) %in% GSE96769_cell_names_matrix$cell_id]
GSE96769_cell_names_matrix = GSE96769_cell_names_matrix[GSE96769_cell_names_matrix$cell_id %in% colnames(GSE96769_scATACseq_matrix),]

# Sum peaks counts relative to celltypes
 
celltype_list = split(GSE96769_cell_names_matrix, GSE96769_cell_names_matrix$cell_types)
celltype_list = celltype_list[c("HSC", "CMP", "GMP")] # there is no MLP

bulk_peaks_list = list()
for(celltype in celltype_list){
	bulk_peaks_df = rowSums(GSE96769_scATACseq_matrix[, celltype$cell_id])
	# threshold = (dim(celltype_list)[1] * 0.25)
	bulk_peaks_df = as.data.frame(as.matrix((bulk_peaks_df > 0) + 0)) # binarize it with t
	bulk_peaks_list = list.append(bulk_peaks_list, bulk_peaks_df)
}
names(bulk_peaks_list) = c("HSC", "CMP", "GMP")


# Split scRNAseq by celltype

# # Binirize the matrix
# GSE117498_scRNAseq_dataset_genes_values = subset(GSE117498_scRNAseq_dataset_genes, select = - gene)
# GSE117498_scRNAseq_dataset_genes_values = as.data.frame(as.matrix((GSE117498_scRNAseq_dataset_genes_values > 0) + 0))
# GSE117498_scRNAseq_dataset_genes = cbind(GSE117498_scRNAseq_dataset_genes$gene, GSE117498_scRNAseq_dataset_genes_values)
# colnames(GSE117498_scRNAseq_dataset_genes)[1] = "gene"

GSE117498_scRNAseq_dataset_genes = dcast(melt(GSE117498_scRNAseq_dataset_genes, id.vars = "gene"), variable ~ gene)
rownames(GSE117498_scRNAseq_dataset_genes) = GSE117498_scRNAseq_dataset_genes$variable
GSE117498_scRNAseq_dataset_genes$cells = gsub('_.*', '', GSE117498_scRNAseq_dataset_genes$variable)
GSE117498_scRNAseq_dataset_genes = GSE117498_scRNAseq_dataset_genes[order(GSE117498_scRNAseq_dataset_genes$cells),]
GSE117498_scRNAseq_dataset_genes$variable = NULL

sc_data_list <- split(GSE117498_scRNAseq_dataset_genes , f = GSE117498_scRNAseq_dataset_genes$cells)
sc_data_list = lapply(sc_data_list, function(x){x$cells = NULL; as.data.frame(t(x))})




# Export processed scRNAseq and scATACseq matrixes for other methods

sc_exp_data = as.data.frame(sc_data_list[c("HSC", "CMP", "GMP")])
sc_exp_data = as.data.frame(t(sc_exp_data))
rownames(sc_exp_data) = sub('[A-Z]{3}\\.', '', rownames(sc_exp_data))

save(sc_exp_data, file = "~/MIPT/epi-impute/data/GSE117498_scRNAseq_genes.Rdata", compress = TRUE)
# write_feather(sc_exp_data, "~/MIPT/epi-impute/data/GSE117498_scRNAseq_genes.feather")


bulk_peaks_data = as.data.frame(bulk_peaks_list[c("HSC", "CMP", "GMP")])
bulk_peaks_data = as.data.frame(t(bulk_peaks_data))
rownames(bulk_peaks_data) = c("HSC", "CMP", "GMP")
colnames(bulk_peaks_data) = colnames(sc_exp_data)

save(bulk_peaks_data, file = "~/MIPT/epi-impute/data/bulk_ATAC_genes.Rdata", compress = TRUE)
# write_feather(bulk_peaks_data, "~/MIPT/epi-impute/data/bulk_ATAC_genes.feather")


# Preprocess bulk data

bulk_exp_data = read_feather('~/MIPT/scRNA_ATAC_task/data/GSE87195_bulk_RNAseq_dataset.feather')
bulk_exp_data = as.data.frame(bulk_exp_data)
# bulk_exp_data = bulk_exp_data[bulk_exp_data$gene %in% homo_sapiens.TFs,]
bulk_exp_data = bulk_exp_data[bulk_exp_data$gene %in% colnames(sc_exp_data),]
rownames(bulk_exp_data) = bulk_exp_data$gene
bulk_genes = bulk_exp_data$gene
bulk_exp_data$gene = NULL
bulk_exp_data = bulk_exp_data[,-1] # drop CLP if needed
bulk_exp_data = as.data.frame(t(bulk_exp_data))

save(bulk_exp_data, file = "~/MIPT/epi-impute/data/GSE87195_bulk_RNAseq_dataset.Rdata", compress = TRUE)
