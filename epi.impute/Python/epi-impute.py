
import pandas as pd
import anndata as ad
import os
import numpy as np
import scipy.sparse as sps
import episcanpy.api as epi
import gc
from pyensembl import EnsemblRelease
from tqdm import tqdm
 

def check_atac_anndata(atac_ann):
    '''
    func to check requared metadata in ATAC anndata object to perform epi-impute process
    '''
    celltype_col_name = 'cell_type'
    assert celltype_col_name in atac_an.obs.columns, f'No needed cell type column found in ATAC file, expected {celltype_col_name}, but found {atac_an.obs.columns}'

def check_rna_anndata(rna_ann, c_types_cell_sums_df):
    '''
    func to check requared metadata in RNA anndata object to perform epi-impute process
    '''
    assert rna_an[:, rna_an.var['gene'].isin(c_types_cell_sums_df.index)].shape[1] != 0, f'No common genes found between RNA var and ATAC annonated peaks'
    celltype_col_name = 'cell_type'
    assert celltype_col_name in rna_an.obs.columns, f'No needed cell type column found in RNA file, expected {celltype_col_name}, but found {rna_an.obs.columns}'

def epi_impute(
    sc_exp_anndata, 
    sc_atac_anndata, 
    gtf_file_path,
    atac_raw = True,
    atac_bin_thrld = 100
):
    '''
    Main function to run imutation. Mimics function, written 
    '''

    # **processing ATAC data
    atac_an = ad.read_h5ad(sc_atac_anndata)
    # check needed metadata
    check_atac_anndata(atac_an)
    # annotate ATAC to get gene names from peaks
    epi.tl.find_genes(
        atac_an,
        gtf_file= gtf_file_path,
        key_added='transcript_annotation',
        upstream=2000,
        feature_type='transcript',
        annotation='HAVANA',
        raw=atac_raw
    )

    # drop regions, that have no connection to genes
    drop_items = ['unassigned', 'intergenic']
    atac_an = atac_an[:, ~atac_an.var.transcript_annotation.isin(drop_items)]

    # split ENSG ids and genes from 'transcript_annotation' to separate columns 'ensg' and 'gene'
    atac_an.var['ensg'] = None
    atac_an.var['gene'] = None
    for row in atac_an.var.index:
        gene = atac_an.var.at[row, 'transcript_annotation']
        g_parts = gene.split(';')
        for g_part in g_parts:
            if 'ENSG' in g_part:
                atac_an.var.at[row, 'ensg'] = g_part
            else:
                atac_an.var.at[row, 'gene'] = g_part
    # drop peaks with None in 'gene' column
    atac_an = atac_an[:, atac_an.var['gene'].notna()]

    # sum peaks over cell types
    c_types_cell_sums_df = pd.DataFrame(columns = atac_an.obs.cell_type.unique())
    c_types_cell_sums_df['gene'] = atac_an.var['gene'].tolist()
    c_types_cell_sums_df = c_types_cell_sums_df.set_index('gene')

    for c_type in atac_an.obs.cell_type.unique():
        c_type_sum = atac_an[atac_an.obs.cell_type == c_type].X.sum(axis = 0).reshape([-1, 1])
        c_types_cell_sums_df[c_type] = np.squeeze(np.asarray(c_type_sum))

    # sum over duplicated genes
    c_types_cell_sums_df = c_types_cell_sums_df.groupby('gene').sum()
    # binarize values by atac_bin_thrld. In original R epi-impute it is equal to 100
    c_types_cell_sums_df[c_types_cell_sums_df < atac_bin_thrld] = 0
    c_types_cell_sums_df[c_types_cell_sums_df >= atac_bin_thrld] = 1

    # then we will need only c_types_cell_sums_df, so we can delete atac anndata
    del atac_an
    gc.collect()
    # **processing RNA file
    rna_an = ad.read_h5ad(sc_exp_anndata)
    # TODO: add check to gene columns
    check_rna_anndata(rna_an, c_types_cell_sums_df)
    # # drop duplicated 
    # rna_an = rna_an[:, ~rna_an.var.duplicated('gene')]
    # intersect genes in ATAC and RNA
    rna_an = rna_an[:, rna_an.var.gene.isin(c_types_cell_sums_df.index)]
    c_types_cell_sums_df = c_types_cell_sums_df[c_types_cell_sums_df.index.isin(rna_an.var.gene)]
    print(f'Common genes count: {rna_an.shape[1]}')

    # convert to lil, does this matter ?
    rna_an.X = rna_an.X.tolil()

    # main imputation loop: adding 0 ar 1 to cell\gene by cell type
    for c_type in tqdm(rna_an.obs.cell_type.unique()):
        rna_an[rna_an.obs.cell_type == c_type].X = \
            rna_an[rna_an.obs.cell_type == c_type].X + c_types_cell_sums_df[c_type].to_numpy()

    # save results
    rna_an.write_h5ad(f'{rna_an}_imputed.h5')




    