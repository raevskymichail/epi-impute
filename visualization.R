sc_exp_data = get(load("/home/mraevsky/MIPT/epi-impute/data/GSE117498_scRNAseq_genes.Rdata"))
sc_atac_cell_names_ <- read_csv("/home/mraevsky/MIPT/data/GSE96769_cell_names_matrix.csv")
sc_atac_data <- read_feather("/home/mraevsky/MIPT/data/GSE96769_scATACseq_matrix.feather")
sc_atac_peaks_ann_ <- read_csv("/home/mraevsky/MIPT/data/GSE96769_PeakFile.csv")



compare_access_exp_of_gene <- function(sc_exp_data, sc_atac_data,
									   celltype, gene,
									   sc_atac_cell_names = sc_atac_cell_names_, 
									   sc_atac_peaks_ann = sc_atac_peaks_ann_){

	if (!require("rlist")) install.packages("rlist")
	if (!require("data.table")) install.packages("data.table")
	# if (!require("plotly")) install.packages("plotly")

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

	# Select only TFs intersected between ATAC and RNA
	genes_common = intersect(sc_atac_genes, sc_exp_data$gene)

	sc_atac_data = sc_atac_data[sc_atac_genes %in% genes_common, ]
	sc_atac_data = sc_atac_data[order(genes_common),]

	sc_exp_data_genes = sc_exp_data[sc_exp_data$gene %in% genes_common, ]
	sc_exp_data_genes = sc_exp_data_genes[order(genes_common),]
	sc_exp_data$gene = NULL


	cells = sc_atac_cell_names[grepl(celltype, sc_atac_cell_names$cell_types),]
	sc_atac_data = as.data.frame(sc_atac_data)
	test_atac = sc_atac_data[genes_common %in% gene, colnames(sc_atac_data) %in% cells$cell_id]

	test_exp = sc_exp_data[genes_common %in% gene, colnames(sc_atac_data) %in% cells$cell_id]

	par(mfrow=c(1,2))
	hist(as.matrix(test_atac), main = "scATAC-seq", xlab = "values")
	hist(as.matrix(test_exp), main = "scRNA-seq", xlab = "values")
}


compare_access_exp_of_gene(sc_exp_data, sc_atac_data, celltype = "HSC", gene = "CD34")