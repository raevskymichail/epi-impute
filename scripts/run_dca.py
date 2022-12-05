def run_dca(sc_data_filename):
	from dca.api import dca
	import numpy as np
	import pandas as pd
	import scanpy as sc
	import time
	count_matrix = pd.read_csv(sc_data_filename, index_col=0)
	count_matrxi = count_matrix.astype('int')
	adata_ae = sc.AnnData(count_matrix.values, obs = count_matrix.index, var = count_matrix.columns)
	sc.pp.filter_genes(adata_ae, min_counts=1)
	# Run imputation
	start = time.time()
	dca(adata_ae, threads=1)
	end = time.time()
	runtime = end - start

	imputed_sc_data = pd.DataFrame(adata_ae.data, 
	                               index=adata_ae.obs[0].values, 
	                               columns=adata_ae.var[0].values)
	return [imputed_sc_data, runtime]

