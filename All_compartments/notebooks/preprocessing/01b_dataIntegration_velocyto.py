import os
import sys
from datetime import datetime
today = datetime.today().strftime('%Y-%m-%d')

import pandas as pd
import numpy as np

import scanpy as sc
import anndata as ad
import hdf5plugin

import argparse

# Add repo path to sys path (allows to access scripts and metadata from repo)
repo_path = '/lustre/scratch126/cellgen/team205/lm25/thymus_projects/thymus_ageing_atlas/General_analysis'
sys.path.insert(1, f'{repo_path}/scripts')

from utils import velocyto_to_anndata

parser = argparse.ArgumentParser()
parser.add_argument('--barcodes_path', type=str, help='Path to the barcodes file.')
parser.add_argument('--velocyto_meta_path', type=str, help='Path to the velocyto meta file.')
parser.add_argument('--out_file_name', type=str, help='Name of the output file.')
parser.add_argument('--n_cpu', type=int, help='Number of CPUs to use.', default=1)

args = parser.parse_args()
barcodes_path = args.barcodes_path
velocyto_meta_path = args.velocyto_meta_path
out_file_name = args.out_file_name
n_cpu = args.n_cpu

# Read barcodes and meta
velocyto_meta = pd.read_csv(velocyto_meta_path, index_col = 0)
barcodes = pd.read_csv(barcodes_path, sep='\t', header=None)[0].tolist()

velocyto_adata = velocyto_to_anndata(meta = velocyto_meta, subset_barcodes=barcodes, n_cpu = n_cpu)

print(f'Velocyto AnnData object shape: {velocyto_adata.shape}')
print(f'Writing AnnData object to {out_file_name}...')

velocyto_adata.write_h5ad(out_file_name,
                compression=hdf5plugin.FILTERS["zstd"],
                compression_opts=hdf5plugin.Zstd(clevel=5).filter_options,
        )