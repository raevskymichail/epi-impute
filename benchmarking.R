# devtools::install_github("eddelbuettel/rbenchmark")
# library(benchmark)


# benchmark("MAGIC" = {})

# splatter - Simple Simulation of Single-cell RNA Sequencing Data

library(feather)
library(reticulate)
use_python('/home/mraevsky/miniconda3/bin/python3.7', required = TRUE)

# library(future)
# plan("multiprocess")

## Default - matrixes - cells-by-gene (cell name include celltype name)
## Only TFs
# sc_exp_data = get(load("~/MIPT/epi-impute/data/GSE117498_scRNAseq_TFs.Rdata"))
# bulk_peaks_data = get(load("~/MIPT/epi-impute/data/bulk_ATAC_TFs.Rdata"))
# bulk_exp_data = get(load("~/MIPT/epi-impute/data/GSE87195_bulk_RNAseq_TFs_dataset.Rdata"))
## All genes
sc_exp_data = get(load("~/MIPT/epi-impute/data/GSE117498_scRNAseq_genes.Rdata"))
bulk_peaks_data = get(load("~/MIPT/epi-impute/data/bulk_ATAC_genes.Rdata"))

# Define list of surface markers

# CD59 and D49f/ITGA6 theoretically should be expressed in HSC but zeros in data
# CD38 should be expressed at CMP and GMP but zeros in data
# CD38 should be expressed at CMP but zeros in data
# CD164 - is leaneage marker - but in data it's not expressed

markers_list = list(HSC = list(plus = c("CD34", "CD59", "ITGA6"), minus = c("CD38", "PTPRC", "FUT4", "TFRC", "ITGA2B", "CD19")),
					CMP = list(plus = c("CD34", "CD38", "FLT3"), minus = c("MME", "PTPRC", "FUT4", "TFRC", "ITGA2B", "CD19")),
					GMP = list(plus = c("CD34", "CD38", "PTPRC"), minus = c("MME", "FLT3", "FUT4", "TFRC", "ITGA2B", "CD19")))



magic_pipeline <- function(data){
	if (!require(viridis)) install.packages("viridis")
	if (!require(ggplot2)) install.packages("ggplot2")
	if (!require(phateR)) install.packages("phateR")
	if (!require(Rmagic)) install.packages("Rmagic")

	runtime = system.time({
		data_imputed = magic(data, genes = "all_genes")
	})

	return(list(result = as.data.frame(data_imputed),
				runtime = runtime))
}


# # Run magic pipeline
# magic_results = magic_pipeline(sc_exp_data)
# save(magic_results, file = './magic_results.Rdata')


scrabble_preprocess <- function(sc_data, bulk_data){
	if (!require(stringr)) install.packages("stringr")

	Ens_TF = read.table("/home/mraevsky/MIPT/scRNA_ATAC_task/INTREGNET_Reproducible/TF_list_hs_AnimalTFDB3_Nov18.txt", sep = "\t", stringsAsFactors=F) # V1 is symbol V2 is EnsID
	homo_sapiens.TFs = as.character(Ens_TF[, 1])

	# Subset only TFs
	bulk_data = as.data.frame(bulk_data)
	bulk_data = bulk_data[bulk_data$gene %in% homo_sapiens.TFs,]
	rownames(bulk_data) = bulk_data$gene
	bulk_genes = bulk_data$gene
	bulk_data$gene = NULL

	bulk_data = bulk_data[,grepl('HSC|CMP|GMP', colnames(bulk_data))]
	# metadata = read_feather('~/MIPT/scRNA_ATAC_task/data/GSE87195_bulk_RNAseq_metadata.feather')

	# Split bulk data by celltype
	genes = rownames(bulk_data)
	cells = colnames(bulk_data)

	bulk_data$gene = rownames(bulk_data)
	bulk_data = dcast(melt(bulk_data, id.vars = "gene"), variable ~ gene)
	rownames(bulk_data) = bulk_data$variable
	bulk_data$cells = str_extract(bulk_data$variable, 'HSC|CMP|GMP')
	bulk_data = bulk_data[order(bulk_data$cells),]
	bulk_data$variable = NULL

	scramble_bulk_data_list <- split(bulk_data, f = bulk_data$cells)
	scramble_bulk_data_list = lapply(scramble_bulk_data_list, function(x){x$cells = NULL; t(x)})

	# Split sc data by celltype
	sc_data = as.data.frame(t(sc_data))
	genes = rownames(sc_data)
	cells = colnames(sc_data)

	sc_data$gene = rownames(sc_data)
	sc_data = dcast(melt(sc_data, id.vars = "gene"), variable ~ gene)
	rownames(sc_data) = sc_data$variable
	sc_data$cells = gsub('_.*', '', sc_data$variable)
	sc_data = sc_data[order(sc_data$cells),]
	sc_data$variable = NULL

	scramble_sc_data_list <- split(sc_data , f = sc_data$cells)
	scramble_sc_data_list = lapply(scramble_sc_data_list, function(x){x$cells = NULL; t(x)})

	return(list(scramble_sc_data_list, scramble_bulk_data_list))
}


scrabble_pipeline <- function(sc_data, bulk_data){
	if (!require("SCRABBLE")) install.packages("SCRABBLE")
	if (!require("parallel")) install.packages("parallel")

	data = scrabble_preprocess(sc_data, bulk_data)
	scramble_sc_data_list = data[[1]]
	scramble_bulk_data_list = data[[2]]

	scrabble_data = list()
	for(i in 1:length(scramble_bulk_data_list)){
	    celltype_list = list(scramble_sc_data_list[[i]], scramble_bulk_data_list[[i]])
	    scrabble_data = list.append(scrabble_data, celltype_list)
	}

	parameter <- c(1, 1e-6, 1e-4)

	runtime = system.time({
		result = mclapply(scrabble_data, function(x){scrabble(x, parameter = parameter)}, mc.cores = detectCores() - 2)
	})
	
	return(list(result = as.data.frame(result),
				runtime = runtime))
}

scrabble_results = scrabble_pipeline(sc_exp_data, bulk_exp_data)
save(scrabble_results, file = "~/MIPT/epi-impute/scrabble_results.Rdata")


saver_pipeline <- function(sc_data){
	if (!require('devtools')) install.packages('devtools')
	if (!require("SAVER")) devtools::install_github("mohuangx/SAVER@*release")
	if (!require("parallel")) install.packages("parallel")

	# convert to gene-by-cells matrix
	sc_data = as.data.frame(t(sc_data))
	runtime = system.time({
		result = saver(sc_data, ncores = detectCores() - 2)
	})
	return(list(result = result,
			runtime = runtime))
}

# Run saver pipeline
saver_results = saver_pipeline(sc_exp_data)
save(saver_results, file = './saver_results.Rdata')


scimpute_pipeline <- function(sc_data){
	if (!require('devtools')) install.packages('devtools')
	if (!require("scImpute")) devtools::install_github("Vivianstats/scImpute")
	if (!require("dplyr")) install.packages("dplyr")
	if (!require("parallel")) install.packages("parallel")

	# convert ot gene-by-cells matrix
	sc_data = as.data.frame(t(sc_data))	

	# sc_data['genes'] = rownames(sc_data)
	# sc_data %>% select(genes, everything()) -> sc_data
	write.table(sc_data, './scimpute_raw_count_data.csv', sep=',')
	labels = gsub("_.*", "", colnames(sc_data))
	ncores = detectCores() - 2

	runtime = system.time({
		scimpute(# full path to raw count matrix
		         count_path = system.file("./scimpute_raw_count_data.csv"), 
		         infile = "csv",           # format of input file
		         outfile = "csv",          # format of output file
		         out_dir = "./",           # full path to output directory
		         labeled = TRUE,           # are cell type labels available
		         labels = labels,          # vectors specifing labels for each cell
		         drop_thre = 0.5,          # threshold set on dropout probability
		         Kcluster = 3,             # 3 cell subpopulations
		         ncores = ncores)          # number of cores used in parallel computation

	})
	result = read.table('./scimpute_count.csv', sep=',')
	file.remove('./scimpute_count.csv')

	return(list(result = result,
				runtime = runtime))
}


drimpute_pipeline <- function(sc_data){
	if (!require('devtools')) install.packages('devtools')
	if (!require("DrImpute")) devtools::install_github("gongx030/DrImpute")

	# convert ot gene-by-cells matrix
	sc_data = as.data.frame(t(sc_data))
	sc_data_log <- log(sc_data + 1)

	runtime = system.time({
		result <- DrImpute(sc_data_log)
	})

	return(list(result = result,
				runtime = runtime))
}


dca_pipeline <- function(sc_data){
	if (!require('reticulate')) install.packages('reticulate')
	source_python("run_dca.py")
	# convert to gene-by-cells matrix
	sc_exp = as.data.frame(t(sc_data))
	write.csv(sc_exp, "tmp_sc_data.csv", quote = FALSE)

	results = run_dca("tmp_sc_data.csv")

	file.remove("tmp_sc_data.csv")

	return(list(result = results[[0]],
				runtime = results[[1]]))
}


# coupledNMF_pipeline <- function(){

# }


define_classes_via_markers <- function(sc_data, .markers_list = markers_list){
	# Generate dataframe filled with FALSE
	sc_data_FALSEs = data.frame(matrix(FALSE, ncol = ncol(sc_data), nrow = nrow(sc_data)),
								   row.names = rownames(sc_data))
	colnames(sc_data_FALSEs) = colnames(sc_data)

	is_positive_class = sc_data_FALSEs
	is_negative_class = sc_data_FALSEs
	is_TP_class = sc_data_FALSEs
	is_TN_class = sc_data_FALSEs

	for (celltype in names(.markers_list)){
		sc_i = grepl(celltype, rownames(sc_data))
		positive_markers = .markers_list[[celltype]]$plus
		negative_markers = .markers_list[[celltype]]$minus
		# positive class (all positive markers expression) -> TRUE, others -> FALSE
		is_positive_class[sc_i, positive_markers] = TRUE
		# negative class (all negative markers expression) -> TRUE, others -> FALSE
		is_negative_class[sc_i, negative_markers] = TRUE
		# TP -> TRUE, drop-outs, others -> FALSE
		is_TP_class[sc_i, positive_markers] = sc_data[sc_i, positive_markers] > 2 # more than 2 counts
		# TN -> TRUE, wrong alignment, small expression -> FALSE
		is_TN_class[sc_i, negative_markers] = sc_data[sc_i, negative_markers] == 0
	}
	# NOTE: is_negative_class is equal is_TN_class, because, it could be not a
	# bad alignment but gene itself could be expressed but we don't see protein on the surface because
	# of regulation on other level (mRNA expressed but protein suppresed)
	return(list(is_positive_class = as.matrix(is_TP_class), # NOTE
				is_TP_class = as.matrix(is_TP_class),
				is_negative_class = as.matrix(is_TN_class), # NOTE
				is_TN_class = as.matrix(is_TN_class)))
}


define_classes_via_bulk <- function(sc_data, bulk_exp_data, .celltypes){
	sc_data = sc_data[, order(colnames(sc_data))]
	bulk_exp_data = bulk_exp_data[, order(colnames(bulk_exp_data))]

	# Consider positive class as all non-zero values in sc RNA data
	# positive -> True, negative & drop-outs -> False
	is_positive_class = sc_data != 0

	# Consider negative class as zeros both in sc and bulk RNA data
	# negative -> True, positive & drop-outs -> False
	for (celltype in .celltypes){
		sc_i = grepl(celltype, rownames(sc_data))
		bulk_i = grepl(celltype, rownames(bulk_exp_data))
		sc_data[sc_i, ] = sc_data[sc_i, ] + colSums(bulk_exp_data[bulk_i, ])
	}

	is_negative_class = sc_data == 0
	return(list(is_positive_class = is_positive_class,
				is_negative_class = is_negative_class))
}


# # Generate datasets of different sparcity from given dataset
# eval_metrics_on_sparced_data <- function(sc_data, simulated_sparcity, bulk_data_for_TN, imputation_func, ...){
# 	if (!require("dplyr")) install.packages("dplyr")

# 	set.seed(42)

# 	initial_sparcity = sum(sc_data == 0) / (dim(sc_data)[1] * dim(sc_data)[2])
# 	print(paste0("Initial sparcity of the datasets is ", round(initial_sparcity, 2), '%'))

# 	# Define Positive and Negative classes

# 	# Consider negative class as zeros both in sc and bulk RNA data
# 	# negative -> True, positive & drop-outs -> False
# 	is_negative_class = define_negative_class_via_bulk(sc_data, bulk_data_for_TN)
# 	negative_vals = sc_data[is_negative_class]
# 	negative_class_size = length(negative_vals)

# 	# Consider positive class as all non-zero values in sc RNA data
# 	is_positive_class = sc_data != 0 # positive -> True, negative & drop-outs -> False
# 	positive_vals = sc_data[is_positive_class]

# 	# generate drop-outs in sc_data
# 	vals_to_dropout = sample(positive_vals, as.integer(simulated_sparcity * length(positive_vals)))
# 	sc_data[is_positive_class][vals_to_dropout] <- 0
# 	n_sim_droupouts = length(vals_to_dropout)

# 	results = imputation_func(sc_data, ...)
# 	imputed_data = results[[1]]
# 	runtime = results[[2]]

# 	# Binarize results (to filter out non-confident predictions)
# 	# imputed_data = round(imputed_data)

# 	# recall (or TPR) is TP/(TP+FN)
# 	true_positive = sum(imputed_data[is_positive_class][vals_to_dropout] != 0)
# 	false_negative = n_sim_droupouts - true_positive
# 	recall = true_positive / n_sim_droupouts

# 	# Specifity = TNR = (1 - FPR) = 1 - FP/(FP + TN) = TN/(FP + TN) 
# 	true_negative = sum(imputed_data[is_negative_class] == 0)
# 	false_positive = negative_class_size - true_negative
# 	specifity = true_negative / negative_class_size

# 	# FDR = 1 - precision = FP/(FP + TP)
# 	false_discovery_rate = false_positive / (false_positive + true_positive)

# 	return(list(recall = recall,
# 				specifity = specifity,
# 				FDR = false_discovery_rate,
# 				result = imputed_data, 
# 				runtime = runtime))
# }


# Generate datasets of different sparcity from given dataset
eval_metrics_on_sparced_data <- function(sc_data, simulated_sparcity, added_droputs_ratio, define_classes_via = c("bulk", "markers"), .markers_list,
																			bulk_data_for_TN, imputation_func, binarize_results = FALSE, atac_bin_thrld, ...){
	if (!require("dplyr")) install.packages("dplyr")

	set.seed(42)
	
	initial_sparcity = sum(sc_data == 0) / prod(dim(sc_data))
	print(paste0("Initial sparcity of the datasets is ", round(initial_sparcity, 2), '%'))
	
	if (!missing(simulated_sparcity)){
		if (initial_sparcity > simulated_sparcity){
			stop("Simulated_sparcity cannot be less than initial sparcity of the input sc_data")
		}	
	}

	# Define Positive and Negative classes
	via_type = match.arg(define_classes_via)
	classes_list = switch(via_type,
					      bulk = define_classes_via_bulk(sc_data, bulk_exp_data = bulk_data_for_TN, ...),
					      markers = define_classes_via_markers(sc_data, ...))

	is_positive_class = classes_list[['is_positive_class']]
	is_negative_class = classes_list[['is_negative_class']]

	positive_vals = sc_data[is_positive_class]

	negative_vals = sc_data[is_negative_class]
	negative_class_size = length(negative_vals)

	# Generate drop-outs in sc_data
	if (missing(simulated_sparcity) & missing(added_droputs_ratio)){
		stop("Specify simulated_sparcity or added_droputs_ratio")
	}

	ifelse(missing(added_droputs_ratio), 
		   added_droputs_ratio = simulated_sparcity - initial_sparcity / length(positive_vals),
		   simulated_sparcity = added_droputs_ratio + initial_sparcity / length(positive_vals)
	)


	vals_to_dropout = sample(positive_vals, as.integer(added_droputs_ratio * length(positive_vals)))
	sc_data[is_positive_class][vals_to_dropout] <- 0
	n_sim_droupouts = length(vals_to_dropout)

	print(paste0("Simulated drop-outs ratio is ", round(added_droputs_ratio, 2), "%"))

	results = imputation_func(sc_data, ...)
	imputed_data = results[[1]]
	runtime = results[[2]]

	# Binarize results (to filter out non-confident predictions)
	if (binarize_results){
		imputed_data = round(imputed_data)
	}

	# recall (or TPR) is TP/(TP+FN)
	true_positive = sum(imputed_data[is_positive_class][vals_to_dropout] != 0)
	false_negative = n_sim_droupouts - true_positive
	recall = true_positive / n_sim_droupouts

	print(paste0("n_sim_droupouts: ", n_sim_droupouts))
	print(paste0("negative_class_size: ", negative_class_size))


	# Specifity = TNR = (1 - FPR) = 1 - FP/(FP + TN) = TN/(FP + TN) 
	true_negative = sum(imputed_data[is_negative_class] == 0)
	false_positive = negative_class_size - true_negative
	specifity = true_negative / negative_class_size

	# FDR = 1 - precision = FP/(FP + TP)
	false_discovery_rate = false_positive / negative_class_size

	return(list(recall = recall,
				specifity = specifity,
				FDR = false_discovery_rate,
				result = imputed_data, 
				runtime = runtime))
}


# for markers only^ subset only markers to speed up becnhmarking
sc_exp_data_markers = sc_exp_data[, unique(unlist(markers_list))]



magic_benchmarking_results = eval_metrics_on_sparced_data(sc_exp_data_markers, 
														  # simulated_sparcity = 0.95,
														  added_droputs_ratio = 0.75,
														  define_classes_via = "markers",
														  .markers_list = markers_list,
														  # .celltypes = c('HSC', 'CMP', 'GMP') 
														  # bulk_data_for_TN = bulk_exp_data, 
														  imputation_func = magic_pipeline,
														  binarize_results = TRUE)
View(magic_benchmarking_results)

epi_impute_benchmarking_results = eval_metrics_on_sparced_data(sc_exp_data_markers, 
															   # simulated_sparcity = 0.95,
															   added_droputs_ratio = 0.75,
															   define_classes_via = "markers",
															   .markers_list = markers_list,
															   # .celltypes = c('HSC', 'CMP', 'GMP')
															   # bulk_data_for_TN = bulk_exp_data,
															   imputation_func = epi_impute,
															   atac_bin_thrld = 100)
View(epi_impute_benchmarking_results)


# Find optimal thresholds
thresholds = c(20, 50, 100, 150, 200)
for (i in 1:length(thresholds)) {
	epi_impute_benchmarking_results = eval_metrics_on_sparced_data(sc_exp_data_markers, 
																   added_droputs_ratio = 0.75,
																   define_classes_via = "markers",
																   .markers_list = markers_list,
																   imputation_func = epi_impute,
																   atac_bin_thrld = thresholds[i])
	print('------------')
	print(paste0("used threshold: ", thresholds[i]))
	print(paste0("recall: ", epi_impute_benchmarking_results[["recall"]]))
	print(paste0("specifity: ", epi_impute_benchmarking_results[["specifity"]]))
	print(paste0("FDR: ", epi_impute_benchmarking_results[["FDR"]]))
	print('------------')
}


# # Generate datasets of different sparcity from given dataset
# eval_recall_on_sparced_data <- function(sc_data, simulated_sparcity, imputation_func){
# 	if (!require("dplyr")) install.packages("dplyr")

# 	initial_sparcity = sum(sc_data == 0) / (dim(sc_data)[1] * dim(sc_data)[2])
# 	print(paste0("Initial sparcity of the datasets is ", round(initial_sparcity, 2), '%'))

# 	is_positive_class = sc_data != 0 # postive = True, negative = False
# 	positive_vals = sc_data[is_positive_class]
# 	positive_class_size = length(positive_vals)

# 	# generate drop-outs in sc_data
# 	vals_to_dropout = sample(positive_vals, as.integer(simulated_sparcity * length(positive_vals)))
# 	sc_data[is_positive_class][vals_to_dropout] <- 0

# 	results = imputation_func(sc_data)
# 	imputed_data = results[[1]]
# 	runtime = results[[2]]

# 	# Binarize results (to filter out non-confident predictions)
# 	imputed_data = round(imputed_data)

# 	# recall (or TPR) is TP/(TP+FN)
# 	true_positive = sum(imputed_data[is_positive_class] != 0)
# 	recall = true_positive / positive_class_size


# 	return(list(recall = recall, 
# 				result = imputed_data, 
# 				runtime = runtime))
# }


# magic_benchmarking_results = eval_recall_on_sparced_data(sc_exp_data, 
# 														 simulated_sparcity = 0.25,  
# 														 magic_pipeline)

# epi_impute_benchmarking_results = eval_recall_on_sparced_data(sc_exp_data, 
# 															  simulated_sparcity = 0.25, 
# 															  epi_impute)




# library(gplots)
# heatmap.2(as.matrix(sc_exp_data), xlab = "cells", ylab = "genes",
# 									            labRow = FALSE, labCol = FALSE,
# 							                    margins = c(2, 2),
# 							                    dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none',
# 									            main = "Origianl data")

# heatmap.2(as.matrix(epi_impute_benchmarking_results[[2]]), xlab = "cells", ylab = "genes",
# 									            labRow = FALSE, labCol = FALSE,
# 							                    margins = c(2, 2),
# 							                    dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none',
# 									            main = "epi_imputed data")

