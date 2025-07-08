# -*- coding: utf-8 -*-
"""
@File    :   build_json.py
@Time    :   2024/09/28 
@Author  :   Dawn
@Version :   1.0
@Desc    :   Build JSON from LRS and SRS data
"""

import os
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy.io import mmread
import pandas as pd
import json
import argparse
import logging


# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


# def _numerical_to_dna(numerical_str, length):
#     """
#     Convert numerical string to DNA sequence.
#     """
#     NT_COMP = {'0': 'A', '1': 'C', '2': 'G', '3': 'T'}
#     numerical = int(numerical_str, 16)  # Convert hex string to integer
#     dna = []

#     for _ in range(length):
#         numerical, remainder = divmod(numerical, 4)
#         dna.append(NT_COMP[str(remainder)])
#     dna = dna[::-1]  # Reverse the list to get the correct order
#     return ''.join(dna)


def _filter_non_zero_indices(adata, batch_size=2000):
    """
    Filters non-zero indices for each column in a DataFrame.
    """
    result = {}

    for j in range(0, adata.shape[1], batch_size):
        data = adata[:, j:j + batch_size].to_df()
        result.update({col: data.index[data[col] != 0].tolist() for col in data.columns})
        del data

    return result


def _invert_dict(data):
    """
    Inverts a dictionary where values are lists.
    """
    transformed_data = {}
    for key, values in data.items():
        for value in values:
            if value not in transformed_data:
                transformed_data[value] = []
            transformed_data[value].append(key)
    return transformed_data


def _generate_adata(path, outdir, sample_type):
    """
    Generates an AnnData object from the input directory path.
    """
    if not os.path.isdir(path):
        raise ValueError(f"Invalid path provided: {path}. Please provide a valid directory path.")

    obs, var, mtx = None, None, None

    for file in os.listdir(path):
        file_path = os.path.join(path, file)
        if not os.path.exists(file_path):
            logging.error(f"File not found: {file_path}")
            continue

        if file == "barcodes.tsv.gz":
            obs = pd.read_csv(file_path, header=None, index_col=0, sep="\t")
            obs.index.name = "barcode"
        elif file == "features.tsv.gz":
            var = pd.read_csv(file_path, header=None, index_col=0, sep="\t")
            var.index.name = "gene"
        elif file == "matrix.mtx.gz":
            mtx = csr_matrix(mmread(file_path).T)

    if obs is None or var is None or mtx is None:
        raise ValueError("Missing one or more required files: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz")

    adata = sc.AnnData(mtx, obs=obs, var=var)
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.filter_cells(adata, min_genes=1)
    adata.write(f'{outdir}/{sample_type}.h5ad')

    return adata


def generate_json(lrs_rawmtx_path, srs_rawmtx_path, outdir, batch_size, min_genes=None):
    """
    Processes input data and generates JSON files.
    """
    logging.info("Generating AnnData objects...")
    l_adata = _generate_adata(lrs_rawmtx_path, outdir, sample_type="lrs")
    s_adata = _generate_adata(srs_rawmtx_path, outdir, sample_type="srs")

    logging.info("Modifying barcode indices...")
    l_adata.obs.index = [i[:10] + i[-10:] for i in l_adata.obs.index]
    # s_adata.obs.index = [_numerical_to_dna(i, 20) for i in s_adata.obs.index]

    logging.info("Filtering barcodes based on min_genes...")
    if min_genes is not None:
        l_adata_trust = l_adata[l_adata.obs['n_genes'] >= min_genes]

        # Save filtered barcodes
        _barcode_list=list(set(l_adata_trust.obs.index))
        best_barcode = {i: i for i in _barcode_list}
        with open(f'{outdir}/l_barcode_trust_dict.json', 'w') as f:
            json.dump(best_barcode, f)
            
        l_adata = l_adata[~l_adata.obs.index.isin(_barcode_list)]

    logging.info("Filtering genes and generating dictionaries...")
    gene_list = list(set(s_adata.var.index) & set(l_adata.var.index))
    l_adata = l_adata[:, l_adata.var.index.isin(gene_list)]
    s_adata = s_adata[:, s_adata.var.index.isin(gene_list)]

    gene_barcode_dict = _filter_non_zero_indices(s_adata, batch_size)
    with open(f'{outdir}/s_gene_barcode_valid.json', 'w') as f:
        json.dump(gene_barcode_dict, f)

    l_gene_barcode_dict = _filter_non_zero_indices(l_adata, batch_size)
    barcode_gene_dict = _invert_dict(l_gene_barcode_dict)
    with open(f'{outdir}/l_barcode_gene_valid.json', 'w') as f:
        json.dump(barcode_gene_dict, f)

    logging.info("JSON files generated successfully.")


def main():
    parser = argparse.ArgumentParser(description='Process data and generate JSON files.')
    parser.add_argument('--lrs_rawmtx_path', type=str, required=True, help='Path to the LRS raw matrix directory')
    parser.add_argument('--srs_rawmtx_path', type=str, required=True, help='Path to the SRS raw matrix directory')
    parser.add_argument('--outdir', type=str, required=True, help='Output directory to store JSON files')
    parser.add_argument('--batch_size', type=int, default=2000, help='Batch size for processing large datasets')
    parser.add_argument('--min_genes', type=int, default=None, help='Minimum number of genes required for a barcode')
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
        logging.info(f"Created output directory: {args.outdir}")

    generate_json(args.lrs_rawmtx_path, args.srs_rawmtx_path, args.outdir, args.batch_size, args.min_genes)


if __name__ == '__main__':
    main()