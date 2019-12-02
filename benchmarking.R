# devtools::install_github("eddelbuettel/rbenchmark")
# library(benchmark)


# benchmark("MAGIC" = {})

# splatter - Simple Simulation of Single-cell RNA Sequencing Data

library(feather)
library(reticulate)
use_python('/home/mraevsky/miniconda3/bin/python3.7', required = TRUE)

# Default - matrixes - cells-by-gene (cell name include celltype name)
sc_exp_data = get(load("~/MIPT/epi-impute/data/GSE117498_scRNAseq_TFs.Rdata"))
bulk_peaks_data = get(load("~/MIPT/epi-impute/data/bulk_ATAC_TFs.Rdata"))
bulk_exp_data = get(load("~/MIPT/epi-impute/data/GSE87195_bulk_RNAseq_dataset.Rdata"))


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



define_negative_class <- function(sc_data, bulk_exp_data){
	sc_data = sc_data[, order(colnames(sc_data))]
	bulk_exp_data = bulk_exp_data[, order(colnames(bulk_exp_data))]

	celltypes = c('HSC', 'CMP', 'GMP')

	for (celltype in celltypes){
		sc_i = grepl(celltype, rownames(sc_data))
		bulk_i = grepl(celltype, rownames(bulk_exp_data))
		sc_data[sc_i, ] = sc_data[sc_i, ] + colSums(bulk_exp_data[bulk_i, ])
	}
	is_negative_class = sc_data == 0
	# negative -> True, positive & drop-outs -> False
	return(is_negative_class)
}


# Generate datasets of different sparcity from given dataset
eval_metrics_on_sparced_data <- function(sc_data, simulated_dropouts_ratio, bulk_data_for_TN, imputation_func, ...){
	if (!require("dplyr")) install.packages("dplyr")

	set.seed(42)

	initial_sparcity = sum(sc_data == 0) / (dim(sc_data)[1] * dim(sc_data)[2])
	print(paste0("Initial sparcity of the datasets is ", round(initial_sparcity, 2), '%'))

	# Define Positive and Negative classes

	# Consider negative class as zeros both in sc and bulk RNA data
	# negative -> True, positive & drop-outs -> False
	is_negative_class = define_negative_class(sc_data, bulk_data_for_TN)
	negative_vals = sc_data[is_negative_class]
	negative_class_size = length(negative_vals)

	# Consider positive class as all non-zero values in sc RNA data
	is_positive_class = sc_data != 0 # positive -> True, negative & drop-outs -> False
	positive_vals = sc_data[is_positive_class]

	# generate drop-outs in sc_data
	vals_to_dropout = sample(positive_vals, as.integer(simulated_dropouts_ratio * length(positive_vals)))
	sc_data[is_positive_class][vals_to_dropout] <- 0
	n_sim_droupouts = length(vals_to_dropout)

	# ifelse(missing(atac_bin_thrld), 
	# 	results = imputation_func(sc_data), 
	# 	results = imputation_func(sc_data, atac_bin_thrld)
	# )

	results = imputation_func(sc_data, ...)
	imputed_data = results[[1]]
	runtime = results[[2]]

	# Binarize results (to filter out non-confident predictions)
	# imputed_data = round(imputed_data)

	# recall (or TPR) is TP/(TP+FN)
	true_positive = sum(imputed_data[is_positive_class][vals_to_dropout] != 0)
	false_negative = n_sim_droupouts - true_positive
	recall = true_positive / n_sim_droupouts

	# Specifity = TNR = (1 - FPR) = 1 - FP/(FP + TN) = TN/(FP + TN) 
	true_negative = sum(imputed_data[is_negative_class] == 0)
	false_positive = negative_class_size - true_negative
	specifity = true_negative / negative_class_size

	# FDR = 1 - precision = FP/(FP + TP)
	false_discovery_rate = false_positive / (false_positive + true_positive)

	return(list(recall = recall,
				specifity = specifity,
				FDR = false_discovery_rate,
				result = imputed_data, 
				runtime = runtime))
}


magic_benchmarking_results = eval_metrics_on_sparced_data(sc_exp_data, 
														  simulated_dropouts_ratio = 0.25, 
														  bulk_data_for_TN = bulk_exp_data, 
														  magic_pipeline)
View(magic_benchmarking_results)

epi-impute_benchmarking_results = eval_metrics_on_sparced_data(sc_exp_data, 
															   simulated_dropouts_ratio = 0.25, 
															   bulk_data_for_TN = bulk_exp_data,
															   epi-impute,
															   atac_bin_thrld = 100)
View(epi-impute_benchmarking_results)


# Find optimal thresholds
thresholds = c(20, 50, 100, 150, 200)
for (i in 1:length(thresholds)) {
	epi-impute_benchmarking_results = eval_metrics_on_sparced_data(sc_exp_data, 
															   simulated_dropouts_ratio = 0.25, 
															   bulk_data_for_TN = bulk_exp_data,
															   epi-impute,
															   atac_bin_thrld = thresholds[i])
	print('------------')
	print(paste0("used threshold: ", thresholds[i]))
	print(paste0("recall: ", epi-impute_benchmarking_results[["recall"]]))
	print(paste0("specifity: ", epi-impute_benchmarking_results[["specifity"]]))
	print(paste0("FDR: ", epi-impute_benchmarking_results[["FDR"]]))
	print('------------')
}


# # Generate datasets of different sparcity from given dataset
# eval_recall_on_sparced_data <- function(sc_data, simulated_dropouts_ratio, imputation_func){
# 	if (!require("dplyr")) install.packages("dplyr")

# 	initial_sparcity = sum(sc_data == 0) / (dim(sc_data)[1] * dim(sc_data)[2])
# 	print(paste0("Initial sparcity of the datasets is ", round(initial_sparcity, 2), '%'))

# 	is_positive_class = sc_data != 0 # postive = True, negative = False
# 	positive_vals = sc_data[is_positive_class]
# 	positive_class_size = length(positive_vals)

# 	# generate drop-outs in sc_data
# 	vals_to_dropout = sample(positive_vals, as.integer(simulated_dropouts_ratio * length(positive_vals)))
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
# 														 simulated_dropouts_ratio = 0.25,  
# 														 magic_pipeline)

# epi-impute_benchmarking_results = eval_recall_on_sparced_data(sc_exp_data, 
# 															  simulated_dropouts_ratio = 0.25, 
# 															  epi-impute)




# library(gplots)
# heatmap.2(as.matrix(sc_exp_data), xlab = "cells", ylab = "genes",
# 									            labRow = FALSE, labCol = FALSE,
# 							                    margins = c(2, 2),
# 							                    dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none',
# 									            main = "Origianl data")

# heatmap.2(as.matrix(epi-impute_benchmarking_results[[2]]), xlab = "cells", ylab = "genes",
# 									            labRow = FALSE, labCol = FALSE,
# 							                    margins = c(2, 2),
# 							                    dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none',
# 									            main = "Epi-imputed data")

