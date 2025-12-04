import os
import sys
import session_info
from datetime import datetime
today = datetime.today().strftime('%Y-%m-%d')

import numpy as np
import pandas as pd
import anndata as ad
import hdf5plugin

# Imoprt custom scvi wrapper
import run_scvi from scvi_wrapper

# Define paths and variables
script_dir = '.' # Path to directory where scvi_wrapper.py is
adata_path = '.h5ad'
object_name = f'DATASETNAME_scvi_{today}'

data_path = '.'
model_path = '.'
plots_path = '.'

# Add script dir to path (otherwise you cannot load the wrapper)
sys.path.append(script_dir)

# Load adata
adata = ad.read_h5ad(adata_path)

# Run scvi
scvi_run = run_scvi(adata, 
                    layer_raw = 'X', 
                    # Excluded genes
                    include_genes=[], exclude_cc_genes=True, exclude_mt_genes=True, 
                    exclude_vdjgenes = True, remove_cite = False,
                    # Highly variable gene selection
                    batch_hv="study", hvg = 3000, span = 0.5,
                    # scVI 
                    batch_scvi="sample",
                    cat_cov_scvi=["donor", "chemistry_simple", "sex"], 
                    cont_cov_scvi=[], 
                    max_epochs=400, batch_size=2000, early_stopping = True, early_stopping_patience = 15, early_stopping_min_delta = 10.0,
                    plan_kwargs = {'lr': 0.001, 'reduce_lr_on_plateau' : True, 'lr_patience' : 10, 'lr_threshold' : 20}, 
                    n_layers = 3, n_latent = 30, dispersion = 'gene-batch',
                    # Leiden clustering
                    leiden_clustering = None,
                    col_cell_type = ['final_anno', 'yoshida_l1', 'crossTissue_l0'], 
                    fig_dir = plots_path, fig_prefix = object_name,
                    )


# Save adata and scvi model
overwrite = True

for c in scvi_run['data'].obs.columns:
    if scvi_run['data'].obs[c].dtype == 'O':
        scvi_run['data'].obs[c] = scvi_run['data'].obs[c].astype('|S')
        
if not os.path.exists(f'{data_path}/{object_name}.zarr') or overwrite:
    scvi_run['data'].write_h5ad(
        f'{data_path}{object_name}.zarr',
        compression=hdf5plugin.FILTERS["zstd"],
        compression_opts=hdf5plugin.Zstd(clevel=5).filter_options,
    )
    scvi_run['vae'].save(f'{model_path}{object_name}', save_anndata=False, overwrite=overwrite)
else:
    print('File already exists')