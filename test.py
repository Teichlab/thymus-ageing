import os
import scipy.io
import pandas as pd
import scanpy as sc
import anndata as ad

prefices = list(set([
    "_".join(f.split("_")[:-2]) for f in os.listdir(".") if f.endswith("_matrix.mtx.gz")
]))
adata_all = []
longest_gene_list = None
for p in prefices:
    mtx = scipy.io.mmread(f"./{p}_matrix.mtx.gz").T.tocsc()
    barcodes = pd.read_csv(f"./{p}_barcodes.tsv.gz", header=None, sep="\t")
    barcodes = [p + "_" + bc for bc in barcodes[0].values]
    genes_file = f"./{p}_genes.tsv.gz"
    features_file = f"./{p}_features.tsv.gz"
    if os.path.exists(genes_file):
        genes = pd.read_csv(genes_file, header=None, sep="\t")
    elif os.path.exists(features_file):
        genes = pd.read_csv(features_file, header=None, sep="\t")
    else:
        raise FileNotFoundError(f"Neither {genes_file} nor {features_file} found")
    
    if longest_gene_list is None or len(genes.iloc[0].tolist()) > len(longest_gene_list.iloc[0].tolist()):
        longest_gene_list = genes
        
    adata = sc.AnnData(mtx)
    adata.var["gene_ids"] = genes[0].values
    adata.var["gene_symbols"] = genes[1].values
    adata.obs["barcodes"] = barcodes
    adata_all.append(adata)
big_adata = ad.concatenate(adata_all, join="outer")
big_adata.var["gene_name"] = longest_gene_list[1].values
big_adata.var["gene_id"] = longest_gene_list[0].values
big_adata.obs["barcode"] = big_adata.obs["barcodes"].astype(str)
big_adata.obs['sample_id'] = [bc.split("_")[0] for bc in big_adata.obs['barcode']]
patient_meta = pd.read_csv("./patient_meta.csv", dtype=str)
patient_meta['Age'] = patient_meta['Age'].astype(float)
big_adata.obs = big_adata.obs.merge(patient_meta, left_on="sample_id", how ='left')

cnt_ls = [c for c in os.listdir() if 'count.csv' in c]
ad_cnt = []
for c in cnt_ls:
    test = pd.read_csv(f"./{c}", index_col =0)
    lib,sample = c.split('.')[0].split('_')
    ad_tmp = ad.AnnData(csr_matrix(test.values.astype(int).T))
    ad_tmp.obs_names = test.columns
    ad_tmp.var_names = test.index
    ad_tmp.obs['library'] = lib
    ad_tmp.obs['sample_id'] = sample
    ad_cnt.append(ad_tmp)
    
cnt_ls = [c for c in os.listdir() if 'counts.mtx' in c]   
ad_cnt = []
for c in cnt_ls:  
    ad_tmp = sc.read_mtx(f'{c}')
    ad_tmp = ad_tmp.T
    lib,sample = c.split('.')[0].split('_')
    genes = pd.read_csv(f'{c.replace("counts.mtx.gz","genes.csv.gz")}').x.tolist()
    barcodes = pd.read_csv(f'{c.replace("counts.mtx.gz","barcode.csv.gz")}').iloc[:,0].tolist()
    ad_tmp.X = csr_matrix(ad_tmp.X.astype(int))
    ad_tmp.obs_names = barcodes
    ad_tmp.var_names = genes
    ad_tmp.obs['library'] = lib
    ad_tmp.obs['sample_id'] = sample
    ad_cnt.append(ad_tmp)
    