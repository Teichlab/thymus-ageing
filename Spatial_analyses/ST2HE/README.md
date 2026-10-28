# ST2HE: H&E morphology niches anchored by spatial transcriptomics

ST2HE turns whole-slide H&E images into **morphological niches**. Tissue is tiled into 56 µm patches and embedded with the
[UNI2-h](https://huggingface.co/MahmoodLab/UNI2-h) pathology foundation model. Post-Xenium and standalone H&E are
integrated with CORAL, and patches are clustered into niches. Xenium single-cell data from paired donors then provides
the cell-type composition of each niche. The niche labels can therefore be applied to H&E-only slides, which have no
spatial transcriptomics.

See *Methods → H&E image processing and integration (ST2HE)* in the manuscript for the full description.

```
WSI (.ndpi / .svs / .ome.tif)
  │  scripts/preprocess_morpho.py         resample to 0.25 mpp → background clean-up → 224 px grid → tile QC → component filter
  ▼
<sample>_processed/*.parquet, *.npy
  │  scripts/generate_embeddings_morpho.py  UNI2-h (bf16) on RGB + grayscale tiles   [GPU]
  ▼
<sample>_THY_tile_all_parts_mpp025_uni2h_morpho.h5ad   (obsm: mf2_tile_uni2h, mf2_tile_uni2h_gray)
  │  04a_morphoNiches_integration.ipynb    concat grayscale embeddings → L2 → PCA(50) → CORAL (paired-donor anchors)
  │                                        → batch-mixing metrics (AUC / iLISI) → bbknn
  ▼
thyAging_CAM_morphoFM_mpp025_mf2_he_coral.h5ad   (released as ThyAge_he_2026-09-05.h5ad)
  │  04b_morphoNiches_discovery.ipynb      Leiden r=2.0 → QC violins + Gini filter → QuPath GeoJSON → annotation
  │                                        → Xenium cell ↔ patch mapping → log2 O/E cell-type enrichment
  ▼
morpho niches (morpho_niche_l1 / l2 / l3)  →  04c (abundance with age), 04d (figures), 05a (patch × gene matrix)
```

## Folder contents

| Path | Purpose |
|---|---|
| `scripts/preprocess_morpho.py` | Reads the slide (OpenSlide or tifffile), resamples it to 0.25 mpp and turns red padding and black background white. It then tiles the image (224 × 224 px, stride 224), applies tile QC and removes small or line-like connected components. |
| `scripts/generate_embeddings_morpho.py` | Computes UNI2-h embeddings for each kept tile and writes one `.h5ad` per slide. `--grayscale` adds the colour-free embeddings used downstream. `--stain-norm` (Macenko) is an optional variant that we tested; it did not remove the platform effect and was not used in the paper. |
| `scripts/stain_norm.py` | Self-contained Macenko stain normaliser. Imported by the embedding script, but only used with `--stain-norm` (not part of the final analysis). |
| `scripts/diagnose_tile_filtering.py` | Optional. Renders per-tile QC metric maps, histograms and galleries of rejected tiles. Used to tune the QC thresholds. |
| `scripts/validate_stain_norm.py` | Optional, exploratory. Before/after check for stain normalisation (CPU). |
| `scripts/plot_morpho.py` | Helper functions. Maps tiles back to original WSI pixel coordinates and exports QuPath GeoJSON. |
| `scripts/run_morpho_{xenium,nonXenium,nih}.sh` | LSF array jobs, one task per slide (preprocess → embed). |
| `scripts/{xenium,nonXenium,nih}_slides.txt` | Slide lists used for the paper, one path per line. The line number is the array index. |
| `scripts/requirements.txt` | Pinned environment for the image-processing scripts (`morpho-env`). |
| `environment_analysis.yml` | Conda environment for the analysis notebooks (exported from the environment used for the paper). |
| `metadata/thyAgeing_HnE.csv` | Per-slide metadata (`DonorID`, `tissue_embedding`, `has_xenium`, `source`, native resolution). |
| `metadata/morpho_niches.xlsx` | Manual annotation of the Leiden r2.0 clusters, made from the QuPath overlays. Columns: `leiden_r2.0` → `morpho_niche_l1/l2/l3`. `morpho_niche_keep = False` marks artifact or donor/sample-specific clusters, which are excluded downstream. The per-donor columns hold the notes made while annotating each slide. |
| `demo/` | Demo data (`data/A72_demo.ome.tif`), demo script (`demo_integration.py`), minimal CPU environment (`requirements_demo.txt`) and reference outputs (`expected_output/`). |
| `04a_morphoNiches_integration.ipynb` | Builds the integrated H&E object (see Reproduction, step 3). |
| `04b_morphoNiches_discovery.ipynb` | Clustering, QC and niche annotation, plus Xenium → niche mapping and enrichment (see Reproduction, step 4). |
| `04c_morphoNiches_sccoda.ipynb` | Changes in niche abundance with age (Milo and scCODA). |
| `04d_morphoNiches_plotting.ipynb` | Manuscript figures: niche composition, density and overlays on H&E and Xenium. |
| `05a_patchExpression.ipynb` | Aggregates Xenium transcripts per H&E patch into a patch × gene `AnnData`. |
| `00_{mcTEC,pDC,inflammatory}_niche.ipynb` | Downstream applications. Xenium cell-neighbourhood niches around selected populations, related back to the morphological niches. |

## 1. System requirements

**Operating system.** All steps were developed and run on Linux: Ubuntu 22.04.5 LTS (kernel 5.15), x86_64. The demo
needs no Linux-specific features, but other operating systems have not been tested.

**Software dependencies.** There are two environments, each with pinned versions:

| Environment | Used for | Python | Key packages (tested versions) |
|---|---|---|---|
| `morpho-env`: `scripts/requirements.txt` | Image processing and UNI2-h embedding (`preprocess_morpho.py`, `generate_embeddings_morpho.py`) | 3.10.20 | torch 2.6.0 (CUDA 12.4), torchvision 0.21.0, timm 1.0.9, huggingface-hub 1.21.0, openslide-python 1.4.6, openslide-bin 4.0.1.2, tifffile 2025.5.10, opencv-python 5.0.0.93, numpy 2.0.2, pandas 2.2.3, anndata 0.11.4 |
| `st2he-analysis`: `environment_analysis.yml` | Integration and niche notebooks (04a–05a, 00_*) | 3.12.13 | scanpy 1.12.1, anndata 0.12.16, bbknn 1.6.0, rapids-singlecell 0.15.0 (RAPIDS 26.4, CUDA 12), scikit-learn 1.8.0, statsmodels 0.14.6, pertpy 1.1.1, spatialdata 0.7.3, spatialdata-io 0.7.0, leidenalg 0.12.0, igraph 1.0.0 |
| `demo/requirements_demo.txt` | CPU-only integration demo | 3.12 | The CPU subset of `st2he-analysis`, with the same versions |

The notebooks also import the project-wide helper modules `utils`, `plotting` and `anndata_compat`. These live in
`General_analysis/scripts` in a parent folder of this module and are added via `sys.path` in each notebook's first
cell. The demo does not need them.

**Hardware.**
- *Demo:* a standard desktop or laptop. The whole demo runs on CPU and needs at most 6 GB of RAM and about 10 GB of
  disk space (environments, UNI2-h weights and the atlas object).
- *Full-size slides (non-standard hardware):* an NVIDIA GPU with CUDA 12 for UNI2-h inference. We tested an NVIDIA
  L40S (48 GB) and an A100 80 GB MIG slice (2g.20gb). At least 10 GB of GPU memory is enough at batch size 384 with
  bf16. The whole slide is read into memory, so each slide needs up to about 55 GB of RAM (peak in our runs: 52 GB).
  On an L40S, one slide (8,000–46,000 tiles) took 8–48 min. The notebooks use RAPIDS (rapids-singlecell) on a CUDA 12
  GPU for PCA, bbknn and Leiden on about 360,000 tiles; the CPU equivalents (see the demo) work but are slower.
- *UNI2-h access:* UNI2-h is a gated model. You need a Hugging Face account that has been granted access at
  <https://huggingface.co/MahmoodLab/UNI2-h>.

## 2. Installation guide

```bash
git clone <this repository>
cd thymus-ageing/Spatial_analyses/ST2HE

# (a) Image processing / embedding environment (morpho-env)
conda create -n morpho-env python=3.10 -y
conda activate morpho-env
pip install torch==2.6.0 torchvision==0.21.0 \
    --index-url https://download.pytorch.org/whl/cu124 --extra-index-url https://pypi.org/simple
pip install -r scripts/requirements.txt
huggingface-cli login            # token of an account with UNI2-h access (or: export HF_TOKEN=...)

# (b) Analysis environment for the notebooks (GPU, CUDA 12)
conda env create -f environment_analysis.yml
conda activate st2he-analysis

# (c) Alternatively, CPU-only environment for the integration demo
conda create -n st2he-demo python=3.12 -y
conda activate st2he-demo
pip install -r demo/requirements_demo.txt
```

The CUDA wheels of torch also run on machines without a GPU, falling back to CPU.

**Typical install time.** Measured with an empty pip cache on our compute cluster:
- (a) took 3.6 min (about 5.7 GB, mostly torch with its CUDA libraries). The UNI2-h weights (2.6 GB) are downloaded
  the first time the embedding script runs.
- (c) took under 1 min (about 0.8 GB).
- (b) is a full export of our environment including RAPIDS, so expect about 15–30 min.

On a typical desktop connection, expect (a) to take about 10–20 min and (c) about 2–5 min; most of this is download
time.

## 3. Demo

The demo has two parts:
- **(A)** runs the image-processing pipeline on a small H&E crop that ships with the repository;
- **(B)** reproduces the CORAL integration and niche discovery on the released H&E tile object, then places the
  crop's tiles into the atlas.

### Demo data

- `demo/data/A72_demo.ome.tif` (28 MB). A lossless 4,050 × 4,050 px crop (0.9 × 0.9 mm, native 0.221 µm/px) of the
  post-Xenium H&E slide of donor A72. Its position in the full slide is recorded in `demo/data/A72_demo_crop.json`.
  The crop was chosen to contain medulla, cortex, septal cortex, vessels, adipose tissue and fibrosis.
- `ThyAge_he_2026-09-05.h5ad` (5.1 GB). The integrated H&E tile object of the paper: 358,603 tiles from 32 slides of
  17 donors. It holds the L2-normalised grayscale UNI2-h embeddings (`X`), `X_pca`, `X_pca_coral`, UMAP, Leiden
  clusters and the `morpho_niche_l1/l2/l3` labels. **It can be downloaded from our website:
  <https://thymusageingatlas.cellgeni.sanger.ac.uk/>.**

### (A) Tile and embed the demo crop (`morpho-env`)

```bash
cd thymus-ageing/Spatial_analyses/ST2HE
mkdir -p demo_run
echo "$PWD/demo/data/A72_demo.ome.tif" > demo_run/slides.txt

python scripts/preprocess_morpho.py         --data-dir demo_run --file-list demo_run/slides.txt
python scripts/generate_embeddings_morpho.py --data-dir demo_run --file-list demo_run/slides.txt --grayscale --batch-size 32
```

**Expected output** in `demo_run/A72_demo.ome_processed/`:
- `A72_demo.ome_part0_he_tiles_mpp025_tile224_stride224.parquet`: 256 tiles (a 16 × 16 grid). All tiles pass QC.
- `A72_demo.ome_THY_tile_all_parts_mpp025_uni2h_morpho.h5ad`: 256 × 1536 embeddings in `obsm['mf2_tile_uni2h']`
  (RGB) and `obsm['mf2_tile_uni2h_gray']` (grayscale). A reference copy is in `demo/expected_output/`.
- QC images, for example `tiles_overlay_crop_mpp025_.png`. Compare it with `demo/expected_output/A72_demo_tiles_overlay.png`.

### (B) Integration and niche discovery (`st2he-demo` or `st2he-analysis`)

```bash
python demo/demo_integration.py \
    --atlas /path/to/ThyAge_he_2026-09-05.h5ad \
    --out-dir demo_run/demo_output \
    --demo-h5ad demo_run/A72_demo.ome_processed/A72_demo.ome_THY_tile_all_parts_mpp025_uni2h_morpho.h5ad
```

`--demo-h5ad` is optional. Without UNI2-h access, point it at the reference copy in
`demo/expected_output/A72_demo.ome_THY_tile_all_parts_mpp025_uni2h_morpho.h5ad`.

The script runs the following steps:
1. Fits CORAL on the five paired donors and maps all post-Xenium tiles.
2. Measures platform mixing within the paired donors before and after CORAL.
3. Builds a bbknn graph on `X_pca_coral`, then runs Leiden (resolution 2.0) and UMAP.
4. Compares the result with the published clusters and niches.
5. Embeds the demo crop's tiles in the atlas PCA + CORAL space and transfers niche labels to them by kNN (k = 15).

**Expected output** in `demo_run/demo_output/`. Reference copies of all files are in `demo/expected_output/`.

`demo_metrics.json`:

| Metric | Expected | Meaning |
|---|---|---|
| Platform AUC within paired donors, `X_pca` → `X_pca_coral` | 1.000 → 0.485 | 1 = platforms fully separable; 0.5 = mixed |
| Balanced iLISI within paired donors, `X_pca` → `X_pca_coral` | 1.000 → 1.135 | 1 = no shared neighbourhoods; 2 = perfect mixing |
| CORAL vs published `X_pca_coral` (Xenium tiles) | r = 0.989 | |
| Leiden r2.0 clusters | 40 | |
| ARI vs published `leiden_r2.0` / niche purity of the demo clusters | 0.354 / 0.730 | |
| Demo crop: tiles matched to the atlas / median cosine similarity to the published embedding | 251 of 256 / 0.984 | |
| Demo crop: kNN niche label agrees with the published niche | 77.3 % | |

The other output files are:
- `umap_demo_integration.png`: UMAP coloured by platform, donor, `morpho_niche_l2` and the demo clusters. Platforms
  should be mixed and niches should form contiguous regions.
- `demo_clusters_vs_morpho_niche_l2.csv`: cross-tabulation of the demo clusters against the published niches.
- `demo_crop_niche_transfer.csv`: predicted and published niche for each demo tile.
- `demo_leiden_labels.csv`: the demo cluster of every tile.

Notes on reproducibility of these numbers:
- The released object only contains the tiles kept after removing donor- and sample-specific clusters. CORAL is
  therefore fitted on a slightly different tile set than in the paper, which is why r is 0.989 rather than 1.
- The paper's graph and clusters were computed on GPU (rapids-singlecell, CAGRA nearest neighbours). The demo uses
  the CPU implementations (bbknn with annoy, igraph Leiden). The clusters are therefore similar to the published ones
  but not identical.
- Leiden and UMAP results can also vary slightly between machines.
- The 5 crop tiles that don't match the atlas are tiles whose published cluster was removed.

### Expected run time

Measured on 8 CPU cores with no GPU. We used a cluster node restricted to 8 cores to approximate a desktop.

| Step | Run time | Peak RAM |
|---|---|---|
| (A) `preprocess_morpho.py` | 12 s | 0.5 GB |
| (A) `generate_embeddings_morpho.py` on CPU (256 tiles, RGB + grayscale, UNI2-h weights cached) | 10 min | 6.0 GB |
| (B) `demo_integration.py` (358,603 tiles, including crop placement) | 8 min | 3.5 GB |

On a GPU, step (A) embedding takes under 1 min. The first run of (A) also downloads the UNI2-h weights (2.6 GB),
and (B) needs the atlas object (5.1 GB) to be downloaded first.

## 4. Instructions for use

### Run ST2HE on your own H&E slides

**1. Prepare a slide list.** Put one slide path per line in a text file. Use absolute paths, or paths relative to
`--data-dir`. Formats can be mixed, because each file is handled according to its suffix:
- `.ome.tif` / `.ome.tiff` are read with tifffile. Use this for post-Xenium H&E exported from the Xenium alignment.
- `.ndpi`, `.svs`, `.tif`, `.mrxs`, etc. are read with OpenSlide.

The sample name is the file stem, with any trailing `" - <scanner timestamp>"` removed. The native resolution is read
from the slide metadata. Slides coarser than 0.25 µm/px are skipped, because they would need upscaling.

**2. Tile and embed each slide** (`morpho-env`, GPU recommended).

Single slide or no scheduler:

```bash
cd ST2HE/scripts
python preprocess_morpho.py         --data-dir $DATA_DIR --file-list my_slides.txt --index 1
python generate_embeddings_morpho.py --data-dir $DATA_DIR --file-list my_slides.txt --index 1 --grayscale --batch-size 384
```

If you omit `--index`, every slide in the list is processed one after another. Both scripts skip slides whose outputs
already exist, so it is safe to re-submit. To run without network access from cached UNI2-h weights, set
`HF_HUB_OFFLINE=1`.

On an LSF cluster, with one array task per slide:

```bash
cd ST2HE/scripts
mkdir -p _logs
export DATA_DIR=/path/to/output_dir      # <sample>_processed/ folders are written here
bsub < run_morpho_xenium.sh              # edit the slide list, #BSUB -q queue and array range [1-N] first
```

Tile QC thresholds can be changed on the command line. The defaults are the values used in the paper:

| Flag | Default | Keeps tile if |
|---|---|---|
| `--valid-tissue-fraction-min` | 0.15 | non-white, non-red fraction ≥ value |
| `--red-fraction-max` | 0.05 | red-artifact fraction ≤ value |
| `--saturation-fraction-min` | 0.01 | colour-saturated fraction ≥ value |
| `--gray-std-min` | 4.5 | grayscale standard deviation ≥ value |
| `--dark-fraction-min` | off | (nuclear fraction; disabled because it removes adipose) |

To inspect what the filters remove on a slide that is already processed, run
`python diagnose_tile_filtering.py --sample-name <sample>`.

**Outputs, written to `$DATA_DIR/<sample>_processed/`:**
- `<sample>_part0_image_he_mpp025.npy`: the resampled, background-cleaned image.
- `<sample>_part0_he_tiles_mpp025_tile224_stride224.parquet`: one row per kept tile, with its grid position, pixel box
  (`x0,x1,y0,y1`) in 0.25-mpp space and its QC metrics.
- `<sample>_THY_tile_all_parts_mpp025_uni2h_morpho.h5ad`: tiles × 1536 UNI2-h embeddings in `obsm`
  (`mf2_tile_uni2h` = RGB, `mf2_tile_uni2h_gray` = grayscale and used downstream; `mf2_tile_uni2h_stainMacenko` only
  with `--stain-norm`).
- QC PNGs: the background-cleaned preview, tile-grid overlays and a random tile gallery.

**3. Assign niches to the new tiles.** Use `place_demo_tiles()` in `demo/demo_integration.py` as a template:
1. L2-normalise the grayscale embeddings.
2. Project them into the atlas PCA space using `varm['PCs']` of `ThyAge_he_2026-09-05.h5ad`.
3. Apply the CORAL map if the slides are post-Xenium H&E. Standalone H&E is the CORAL target and is not transformed.
4. Transfer `morpho_niche_l1/l2/l3` by kNN against `X_pca_coral`.

Alternatively, to discover niches de novo in a new cohort, follow steps 3–4 of the reproduction instructions below.

### Reproduction instructions

These notebooks reproduce the analyses in the paper (`st2he-analysis` environment, CUDA GPU). Before running them,
change the absolute paths in each notebook's first cell (`repo_path`, `data_path`, `plots_path`) to your local copies.

**Steps 1–2.** Run the three LSF scripts on the slide lists in `scripts/`: 8 post-Xenium slides, 38 standalone H&E
slides and 6 NIH slides.

**Step 3. Integration (`04a_morphoNiches_integration.ipynb`).**
1. Collect every `*_all_parts_mpp025_uni2h_morpho.h5ad`, take the **grayscale** embedding (`mf2_tile_uni2h_gray`) and
   attach the slide metadata from `metadata/thyAgeing_HnE.csv`.
2. Remove OCT-embedded and NIH slides. Only FFPE slides from the CAM cohort are clustered.
3. L2-normalise each tile, then run PCA with 50 components.
4. **CORAL**: the five donors with both post-Xenium and standalone H&E (A67, A70, A72, U48, A45) serve as anchors.
   Estimate the mean and covariance (+10⁻⁴·I) of the paired Xenium tiles (source) and the paired non-Xenium tiles
   (target). Then map *all* Xenium tiles with `Z = (Z − μ_s)·C_s^{-1/2}·C_t^{1/2} + μ_t`. The result is stored as
   `obsm['X_pca_coral']`.
5. Check the integration with platform AUC and iLISI within paired donors. For transparency, the notebook keeps the
   corrections we tried before CORAL: Harmony, a global mean-shift, location/scale matching and an MNN step on top of
   CORAL. None of them fixed the platform separation within paired donors, or they distorted the embedding, so they
   were not used. Only `X_pca_coral` is used downstream.
6. Build a bbknn graph on `X_pca_coral` (`batch_key='has_xenium'`) and compute a UMAP.

Output: `morphoFM/thyAging_CAM_morphoFM_mpp025_mf2_he_coral.h5ad`.

**Step 4. Niche discovery and annotation (`04b_morphoNiches_discovery.ipynb`).**
1. Leiden clustering at resolution 2.0 (`leiden_r2.0`).
2. Plot violins of the tile QC metrics for each cluster. Compute Gini coefficients of donor and sample contributions per
   cluster. Clusters dominated by a single donor or sample are removed.
3. Export the clusters as QuPath GeoJSON polygons in original WSI pixel coordinates
   (`plot_morpho.export_niches_to_geojson`). Load them in QuPath on top of the slide and annotate each cluster from its
   histology. Record the labels in `metadata/morpho_niches.xlsx`, which the notebook merges back into the object.
4. **Xenium → niche**: for the paired donors, Xenium cell coordinates are transformed into H&E pixel space with the
   Xenium alignment matrix (`<donor>_alignment_files/matrix.csv`). Each cell is assigned to its 224-px patch by integer
   grid indexing (`patch_col`, `patch_row`). Output: `xenium/thyAging_xenium_morpho_niches.csv`.
5. **Enrichment**: compute log₂ observed/expected per donor for each cell type in each niche. Combine donors with a
   one-sample t-test (H₀: log₂ O/E = 0) and correct with BH-FDR.

**Step 5. Downstream analyses.**
- `04c_morphoNiches_sccoda.ipynb`: changes in niche abundance with age, tested with Milo and scCODA.
- `04d_morphoNiches_plotting.ipynb`: manuscript figure panels.
- `05a_patchExpression.ipynb`: patch × gene matrix (`thyAging_CAM_morphoFM_mpp025_mf2_patchExpr.h5ad`) built from the
  Xenium `spatialdata` zarrs and the alignment matrices.
- `00_*_niche.ipynb`: neighbourhood niches of selected Xenium populations, cross-referenced to the morphological niches.

## Data availability

The integrated H&E tile object (`ThyAge_he_2026-09-05.h5ad`) is available for download at
<https://thymusageingatlas.cellgeni.sanger.ac.uk/>. The raw H&E images and Xenium data are available as described in
the manuscript's *Data availability* section. The paths in the slide lists, `metadata/thyAgeing_HnE.csv` and the
notebooks point to our compute environment. They are kept for illustration and must be changed to your local copies.

## Citation

If you use ST2HE, please cite the manuscript *A Human Thymus Ageing Cell Atlas Defines A Spatially Resolved Model of
Involution*, together with UNI (Chen et al., *Nat. Med.* 2024), CORAL (Sun et al., *AAAI* 2016) and bbknn
(Polański et al., *Bioinformatics* 2020).

## Contributors

ST2HE was developed by Lisa Marie Milchsack, with code contributions from Hansheng Xue and supervision by Nadav Yayon.

## License

ST2HE is released under the BSD 3-Clause License; see [`LICENSE`](LICENSE). This covers the code in this folder only.
Third-party models and tools it uses, such as the UNI2 weights, are distributed under their own licenses.
