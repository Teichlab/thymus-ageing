#!/usr/bin/env python
"""
ST2HE integration demo (CPU).

Re-runs the cross-platform integration and niche discovery of
04a_morphoNiches_integration.ipynb / 04b_morphoNiches_discovery.ipynb on the released
H&E tile object (ThyAge_he_2026-09-05.h5ad, https://thymusageingatlas.cellgeni.sanger.ac.uk/):

  X_pca (PCA(50) of L2-normalised grayscale UNI2-h tile embeddings)
    -> CORAL alignment of post-Xenium to standalone H&E, anchored on the paired donors
    -> within-paired-donor platform AUC / iLISI before vs after CORAL
    -> batch-balanced kNN (bbknn) -> Leiden (resolution 2.0) -> UMAP
    -> agreement with the published clusters / morphological niches

Optionally (--demo-h5ad) it also places the tiles embedded from the demo H&E crop
(demo/data/A72_demo.ome.tif) into the atlas: their embeddings are compared with the
atlas embeddings of the same tiles, projected into the atlas PCA + CORAL space, and
given a niche label by kNN label transfer.

The released object contains only the tiles kept after removing donor-/sample-specific
clusters, so CORAL statistics and clusters are close to, but not bit-identical with,
the published ones.
"""
import argparse
import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import adjusted_rand_score
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler, normalize

PAIRED_DONORS = ['A67', 'A70', 'A72', 'U48', 'A45']   # donors with post-Xenium AND standalone H&E
ATLAS_URL = 'https://thymusageingatlas.cellgeni.sanger.ac.uk/'


# ---- CORAL (as in 04a) -------------------------------------------------------------
def sym_pow(C, p, eps=1e-6):
    w, V = np.linalg.eigh(C)
    w = np.clip(w, eps, None)
    return (V * (w ** p)) @ V.T


def fit_coral(Z, is_xen, donors, paired=PAIRED_DONORS, reg=1e-4):
    """Source = paired-donor Xenium tiles, target = paired-donor non-Xenium tiles."""
    xp = np.isin(donors, paired) & is_xen
    npd = np.isin(donors, paired) & ~is_xen
    mu_s, mu_t = Z[xp].mean(0), Z[npd].mean(0)
    Cs = np.cov(Z[xp], rowvar=False) + reg * np.eye(Z.shape[1])
    Ct = np.cov(Z[npd], rowvar=False) + reg * np.eye(Z.shape[1])
    A = sym_pow(Cs, -0.5) @ sym_pow(Ct, 0.5)          # whiten by source, recolour to target
    return mu_s, mu_t, A


def apply_coral(Z, is_xen, mu_s, mu_t, A):
    Z = Z.copy()
    Z[is_xen] = (Z[is_xen] - mu_s) @ A + mu_t
    return Z


# ---- integration metrics (as in 04a) -----------------------------------------------
def ilisi(rep, batch, k=30, subsample=50000, seed=0):
    """Mean inverse-Simpson index of the batch label over each tile's kNN (2 batches: 1..2)."""
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(rep), min(subsample, len(rep)), replace=False)
    R, b = rep[idx], np.asarray(batch)[idx]
    neigh = NearestNeighbors(n_neighbors=k + 1).fit(R).kneighbors(R, return_distance=False)[:, 1:]
    labs = b[neigh]
    p = np.stack([(labs == v).mean(axis=1) for v in np.unique(b)], axis=1)
    return float(np.mean(1.0 / np.sum(p ** 2, axis=1)))


def paired_auc(X, y, seed=0):
    """Within-paired-donor platform separability (0.5 = mixed, 1 = separable)."""
    X = StandardScaler().fit_transform(X)
    return float(cross_val_score(LogisticRegression(max_iter=1000), X, y,
                                 cv=StratifiedKFold(5, shuffle=True, random_state=seed),
                                 scoring='roc_auc').mean())


def paired_ilisi_balanced(X, y, k=30, seed=0):
    rng = np.random.default_rng(seed)
    n = min((y == 0).sum(), (y == 1).sum())
    idx = np.concatenate([rng.choice(np.where(y == 0)[0], n, False),
                          rng.choice(np.where(y == 1)[0], n, False)])
    return ilisi(X[idx], y[idx], k=k, subsample=len(idx))


def log(msg, t0):
    print(f'[{time.time() - t0:7.1f}s] {msg}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--atlas', type=Path, required=True,
                        help=f'ThyAge_he_2026-09-05.h5ad (download from {ATLAS_URL}).')
    parser.add_argument('--out-dir', type=Path, default=Path('demo_output'))
    parser.add_argument('--demo-h5ad', type=Path, default=None,
                        help='Optional: <sample>_THY_tile_all_parts_mpp025_uni2h_morpho.h5ad produced '
                             'from demo/data/A72_demo.ome.tif by the pipeline scripts.')
    parser.add_argument('--resolution', type=float, default=2.0)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = args.out_dir
    t0 = time.time()

    # Only obs/obsm are needed for the integration; X (tiles x 1536) is read lazily.
    atlas = ad.read_h5ad(args.atlas, backed='r')
    adata = ad.AnnData(obs=atlas.obs[['SampleID', 'DonorID', 'age_group', 'has_xenium', 'leiden_r2.0',
                                      'morpho_niche_l1', 'morpho_niche_l2', 'morpho_niche_l3',
                                      'grid_row', 'grid_col']].copy(),
                       obsm={'X_pca': np.asarray(atlas.obsm['X_pca'], dtype=np.float64),
                             'X_pca_coral_published': np.asarray(atlas.obsm['X_pca_coral'])})
    obs = adata.obs
    is_xen = obs['has_xenium'].astype(bool).values
    donors = obs['DonorID'].astype(str).values
    log(f'Loaded {adata.n_obs:,} tiles from {obs.SampleID.nunique()} slides / {obs.DonorID.nunique()} donors '
        f'({is_xen.sum():,} post-Xenium, {(~is_xen).sum():,} standalone H&E)', t0)

    # 1) CORAL
    mu_s, mu_t, A = fit_coral(adata.obsm['X_pca'], is_xen, donors)
    adata.obsm['X_pca_coral'] = apply_coral(adata.obsm['X_pca'], is_xen, mu_s, mu_t, A).astype(np.float32)
    r = np.corrcoef(adata.obsm['X_pca_coral'][is_xen].ravel(),
                    adata.obsm['X_pca_coral_published'][is_xen].ravel())[0, 1]
    log(f'CORAL done (correlation with published X_pca_coral on Xenium tiles: r = {r:.3f})', t0)

    # 2) Batch-mixing metrics within the paired donors
    m = np.isin(donors, PAIRED_DONORS)
    y = is_xen[m].astype(int)
    metrics = {}
    for rep in ['X_pca', 'X_pca_coral']:
        X = np.asarray(adata.obsm[rep])[m]
        metrics[rep] = {'paired_donor_platform_AUC': round(paired_auc(X, y, args.seed), 3),
                        'paired_donor_balanced_iLISI': round(paired_ilisi_balanced(X, y, seed=args.seed), 3)}
        log(f'{rep:12s} AUC = {metrics[rep]["paired_donor_platform_AUC"]:.3f}  '
            f'iLISI = {metrics[rep]["paired_donor_balanced_iLISI"]:.3f}  (mixed: AUC 0.5, iLISI 2.0)', t0)

    # 3) bbknn graph on the CORAL representation -> Leiden -> UMAP
    sc.external.pp.bbknn(adata, batch_key='has_xenium', use_rep='X_pca_coral',
                         neighbors_within_batch=3, n_pcs=50)
    log('bbknn graph built', t0)
    key = f'leiden_r{args.resolution}_demo'
    sc.tl.leiden(adata, resolution=args.resolution, key_added=key, random_state=args.seed,
                 flavor='igraph', n_iterations=2, directed=False)
    log(f'Leiden r={args.resolution}: {obs[key].nunique()} clusters', t0)
    sc.tl.umap(adata, min_dist=0.3, random_state=args.seed)
    log('UMAP done', t0)

    ari_cluster = adjusted_rand_score(obs['leiden_r2.0'], obs[key])
    ari_niche = adjusted_rand_score(obs['morpho_niche_l2'], obs[key])
    # Purity: share of tiles whose demo cluster's majority niche equals their published niche.
    maj = obs.groupby(key, observed=True)['morpho_niche_l2'].agg(lambda s: s.value_counts().index[0])
    purity = float((obs[key].map(maj).astype(str) == obs['morpho_niche_l2'].astype(str)).mean())
    metrics['clustering'] = {'n_clusters': int(obs[key].nunique()),
                             'ARI_vs_published_leiden_r2.0': round(ari_cluster, 3),
                             'ARI_vs_morpho_niche_l2': round(ari_niche, 3),
                             'niche_l2_purity_of_demo_clusters': round(purity, 3)}
    log(f'Agreement with published clustering: ARI = {ari_cluster:.3f}; '
        f'niche_l2 purity of demo clusters = {purity:.3f}', t0)

    pd.crosstab(obs[key], obs['morpho_niche_l2']).to_csv(args.out_dir / 'demo_clusters_vs_morpho_niche_l2.csv')
    sc.pl.umap(adata, color=['has_xenium', 'DonorID', 'morpho_niche_l2', key], ncols=2, wspace=0.4,
               show=False, save='_demo_integration.png')

    # 4) Optional: place the demo crop's tiles into the atlas
    if args.demo_h5ad is not None:
        metrics['demo_crop'] = place_demo_tiles(args, atlas, adata, mu_s, mu_t, A, t0)

    with open(args.out_dir / 'demo_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)
    adata.obs[[key]].to_csv(args.out_dir / 'demo_leiden_labels.csv')
    log(f'Done. Outputs in {args.out_dir}/', t0)
    print(json.dumps(metrics, indent=2))


def place_demo_tiles(args, atlas, adata, mu_s, mu_t, A, t0, k=15):
    """Match demo-crop tiles to the same A72.ome tiles in the atlas, then kNN-transfer niche labels."""
    demo = ad.read_h5ad(args.demo_h5ad)
    Zd = normalize(np.asarray(demo.obsm['mf2_tile_uni2h_gray'], dtype=np.float64), norm='l2')

    # The demo crop starts at atlas grid tile (row 123, col 128) of slide A72.ome.
    meta = json.loads(Path(__file__).with_name('data').joinpath('A72_demo_crop.json').read_text())
    rows = demo.obs['grid_row'].to_numpy() + meta['atlas_grid_row0']
    cols = demo.obs['grid_col'].to_numpy() + meta['atlas_grid_col0']
    a72 = adata.obs[adata.obs['SampleID'] == meta['atlas_sample']]
    pos = pd.Series(np.arange(adata.n_obs), index=adata.obs_names)
    lut = {(r, c): pos[i] for i, r, c in zip(a72.index, a72['grid_row'], a72['grid_col'])}
    hit = np.array([(r, c) in lut for r, c in zip(rows, cols)])
    idx = np.array([lut[(r, c)] for r, c, h in zip(rows, cols, hit) if h])

    Xa = atlas.X[np.sort(idx)]
    Xa = np.asarray(Xa.todense() if hasattr(Xa, 'todense') else Xa, dtype=np.float64)
    Xa = Xa[np.argsort(np.argsort(idx))]                                 # back to demo order
    cos = np.sum(Zd[hit] * Xa, axis=1)
    log(f'Demo crop: {len(Zd)} tiles, {hit.sum()} also in atlas (rest removed by atlas QC); '
        f'median cosine similarity to published embedding = {np.median(cos):.3f}', t0)

    # Project into the atlas PCA space. X_pca = (X - m) @ PCs, so the centring term m @ PCs is
    # recovered exactly as the mean of (X @ PCs - X_pca) over the atlas tiles.
    PCs = np.asarray(atlas.varm['PCs'], dtype=np.float64)
    offset = np.zeros(PCs.shape[1])
    for s in range(0, atlas.n_obs, 50000):
        blk = atlas.X[s:s + 50000]
        offset += np.asarray(blk @ PCs - adata.obsm['X_pca'][s:s + 50000]).sum(axis=0)
    offset /= atlas.n_obs
    Zp = Zd @ PCs - offset
    Zc = (Zp - mu_s) @ A + mu_t                                          # demo crop is post-Xenium H&E

    nn = NearestNeighbors(n_neighbors=k).fit(adata.obsm['X_pca_coral'])
    neigh = nn.kneighbors(Zc, return_distance=False)
    niches = adata.obs['morpho_niche_l2'].astype(str).to_numpy()
    pred = np.array([pd.Series(niches[n]).value_counts().index[0] for n in neigh])
    truth = niches[idx]
    acc = float((pred[hit] == truth).mean())
    out = demo.obs[['grid_row', 'grid_col']].copy()
    out['predicted_morpho_niche_l2'] = pred
    out['published_morpho_niche_l2'] = pd.Series(truth, index=out.index[hit])
    out.to_csv(args.out_dir / 'demo_crop_niche_transfer.csv')
    log(f'kNN niche transfer (k={k}) agrees with published niche for {acc:.1%} of matched tiles', t0)
    return {'n_tiles': int(len(Zd)), 'n_matched_to_atlas': int(hit.sum()),
            'median_cosine_to_published_embedding': round(float(np.median(cos)), 3),
            'niche_transfer_agreement': round(acc, 3)}


if __name__ == '__main__':
    main()
