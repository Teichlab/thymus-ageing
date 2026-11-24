#!/usr/bin/env python

import os
import numpy as np
import anndata as ad
import scanpy as sc
import scvi
import torch

def main():
    scvi.settings.seed = 0
    torch.set_float32_matmul_precision("high")

    adata = ad.read_h5ad(
        "/home/nd590/rfs/rfs-teichlab-nfs-iCNyzSAaucw/nd590/projects/aging/objects/xenium_sc_sn_proseg_aging.h5ad"
    )

    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key = 'cells_spatial',
        categorical_covariate_keys=['donor','chemistry_simple']
    )

    xsc_model = scvi.model.SCVI(
        adata,
        n_latent=30,
        n_layers=2,
        n_hidden=256,
        dispersion="gene-batch",
        gene_likelihood="nb",
    )

    for var in ["SLURM_NTASKS", "SLURM_JOB_ID", "SLURM_NODELIST"]:
        os.environ.pop(var, None)

    xsc_model.train(
        accelerator="gpu",
        devices=1,
        batch_size=2000,
        max_epochs=400,
        early_stopping=True,
    )

    model_dir = (
        "/home/nd590/rfs/rfs-teichlab-nfs-iCNyzSAaucw/nd590/projects/aging/models/scvi_sc_sn_xenium_proseg_aging"
    )
    xsc_model.save(model_dir, overwrite=True)


if __name__ == "__main__":
    main()
