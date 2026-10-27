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
<sample>_THY_tile_all_parts_mpp025_uni2h_morpho.h5ad   (obsm: mf2_tile_uni2h, mf2_tile_uni2h_gray, …)
  │  04a_morphoNiches_integration.ipynb    concat grayscale embeddings → L2 → PCA(50) → CORAL (paired-donor anchors)
  │                                        → batch-mixing metrics (AUC / iLISI) → bbknn
  ▼
thyAging_CAM_morphoFM_mpp025_mf2_he_coral.h5ad
  │  04b_morphoNiches_discovery.ipynb      Leiden r=2.0 → QC violins + Gini filter → QuPath GeoJSON → annotation
  │                                        → Xenium cell ↔ patch mapping → log2 O/E cell-type enrichment
  ▼
morpho niches (morpho_niche_l2 / l3)  →  04c (abundance with age), 04d (figures), 05a (patch × gene matrix)
```

## Folder contents

| Path                                             | Purpose                                                                                                                                                                                                                                                          |
| ------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `scripts/preprocess_morpho.py`                 | Reads the slide (OpenSlide or tifffile), resamples it to 0.25 mpp and turns red padding and black background white. It then tiles the image (224 × 224 px, stride 224), applies tile QC and removes small or line-like connected components.                    |
| `scripts/generate_embeddings_morpho.py`        | Computes UNI2-h embeddings for each kept tile and writes one `.h5ad` per slide. `--grayscale` adds the colour-free embeddings used downstream. `--stain-norm` (Macenko) is an optional variant that we tested; it did not remove the platform effect and was not used in the paper. |
| `scripts/stain_norm.py`                        | Self-contained Macenko stain normaliser. Imported by the embedding script, but only used with `--stain-norm` (not part of the final analysis).                                                                                                                    |
| `scripts/diagnose_tile_filtering.py`           | Optional. Renders per-tile QC metric maps, histograms and galleries of rejected tiles. Used to tune the QC thresholds.                                                                                                                                           |
| `scripts/validate_stain_norm.py`               | Optional, exploratory. Before/after check for stain normalisation (CPU).                                                                                                                                                                                                   |
| `scripts/plot_morpho.py`                       | Helper functions. Maps tiles back to original WSI pixel coordinates and exports QuPath GeoJSON.                                                                                                                                                                  |
| `scripts/run_morpho_{xenium,nonXenium,nih}.sh` | LSF array jobs, one task per slide (preprocess → embed).                                                                                                                                                                                                        |
| `scripts/{xenium,nonXenium,nih}_slides.txt`    | Slide lists, one path per line. The line number is the array index.                                                                                                                                                                                              |
| `scripts/requirements.txt`                     | Pinned environment for the two pipeline scripts.                                                                                                                                                                                                                 |
| `metadata/thyAgeing_HnE.csv`                   | Per-slide metadata (`DonorID`, `tissue_embedding`, `has_xenium`, `source`, native resolution).                                                                                                                                                           |
| `metadata/morpho_niches.xlsx`                  | Manual annotation of the Leiden r2.0 clusters, made from the QuPath overlays. Columns: `leiden_r2.0` → `morpho_niche_l1/l2/l3`. `morpho_niche_keep = False` marks artifact or donor/sample-specific clusters, which are excluded downstream. The per-donor columns hold the notes made while annotating each slide. |
| `04a_morphoNiches_integration.ipynb`           | Builds the integrated H&E object (see step 3).                                                                                                                                                                                                                   |
| `04b_morphoNiches_discovery.ipynb`             | Clustering, QC and niche annotation, plus Xenium → niche mapping and enrichment (see step 4).                                                                                                                                                                   |
| `04c_morphoNiches_sccoda.ipynb`                | Changes in niche abundance with age (Milo and scCODA).                                                                                                                                                                                                           |
| `04d_morphoNiches_plotting.ipynb`              | Manuscript figures: niche composition, density and overlays on H&E and Xenium.                                                                                                                                                                                   |
| `05a_patchExpression.ipynb`                    | Aggregates Xenium transcripts per H&E patch into a patch × gene `AnnData`.                                                                                                                                                                                    |
| `00_{mcTEC,pDC,inflammatory}_niche.ipynb`      | Downstream applications. Xenium cell-neighbourhood niches around selected populations, related back to the morphological niches.                                                                                                                                 |

## Requirements

**Pipeline scripts (steps 1–2):** a Linux machine with a CUDA GPU (≥10 GB; bf16 autocast) and enough RAM to hold a
whole slide at full resolution (we requested 80 GB per slide).

```bash
conda create -n morpho-env python=3.10 -y
conda activate morpho-env
pip install torch==2.6.0 torchvision==0.21.0 --index-url https://download.pytorch.org/whl/cu124
pip install -r scripts/requirements.txt
```

UNI2-h is a gated model. Request access at [https://huggingface.co/MahmoodLab/UNI2-h](https://huggingface.co/MahmoodLab/UNI2-h) and then authenticate with
`huggingface-cli login` or by setting `HF_TOKEN`. Running `generate_embeddings_morpho.py` calls `huggingface_hub.login()`,
which uses the cached token.

**Notebooks (steps 3–5):** these were run on a GPU with Python 3.12 using the following versions: scanpy 1.12.1,
anndata 0.12.16, rapids-singlecell 0.15.0 (GPU PCA, bbknn and Leiden), bbknn 1.6.0, scikit-learn 1.8.0,
statsmodels 0.14.6, pertpy 1.1.1, spatialdata 0.7.3 and spatialdata-io 0.7.0. On a machine without a GPU, replace
`rsc.pp.pca`, `rsc.pp.bbknn`, `rsc.tl.leiden` and `rsc.tl.umap` with the `scanpy` / `bbknn` equivalents.

The notebooks also import the project-wide helper modules `utils`, `plotting` and `anndata_compat`, which live in
`General_analysis/scripts` in a parent folder of this module and are added via `sys.path` in each notebook's first
cell. They also use absolute paths defined in that first cell (`repo_path`, `data_path`,
`plots_path`). Change these paths to your local layout before running.

## How to run

### 1. Prepare the slide lists

Put one slide path per line in a `*_slides.txt` file. Use absolute paths, or paths relative to `--data-dir`. Formats
can be mixed, because each file is handled according to its suffix:

- `.ome.tif` / `.ome.tiff` are read with tifffile. Use this for post-Xenium H&E exported from the Xenium alignment.
- `.ndpi`, `.svs`, `.tif`, `.mrxs`, etc. are read with OpenSlide.

The sample name is the file stem, with any trailing `" - <scanner timestamp>"` removed.
The lists used for the paper are included: 8 post-Xenium slides, 38 standalone H&E slides and 6 NIH slides.

### 2. Tile and embed each slide

**On LSF (one array task per slide):**

```bash
cd ST2HE/scripts
mkdir -p _logs
export DATA_DIR=/path/to/output_dir      # <sample>_processed/ folders are written here
bsub < run_morpho_xenium.sh              # array [1-8]
bsub < run_morpho_nonXenium.sh           # array [1-38]
bsub < run_morpho_nih.sh                 # array [1-6]
```

Edit the `#BSUB -q` queue name and the array range (`[1-N]`, where N is the number of lines in the list) to match
your cluster and slide list.

**Single slide / no scheduler:**

```bash
cd ST2HE/scripts
python preprocess_morpho.py         --data-dir $DATA_DIR --file-list xenium_slides.txt --index 1
python generate_embeddings_morpho.py --data-dir $DATA_DIR --file-list xenium_slides.txt --index 1 --grayscale --batch-size 384
```

If you omit `--index`, every slide in the list is processed one after another. Both scripts skip slides whose outputs
already exist, so it is safe to re-submit. Slides whose native resolution is coarser than 0.25 mpp are skipped, because
they would need upscaling.

Tile QC thresholds can be changed on the command line. The defaults are the values used in the paper:

| Flag                            | Default | Keeps tile if                                           |
| ------------------------------- | ------- | ------------------------------------------------------- |
| `--valid-tissue-fraction-min` | 0.15    | non-white, non-red fraction ≥ value                    |
| `--red-fraction-max`          | 0.05    | red-artifact fraction ≤ value                          |
| `--saturation-fraction-min`   | 0.01    | colour-saturated fraction ≥ value                      |
| `--gray-std-min`              | 4.5     | grayscale standard deviation ≥ value                   |
| `--dark-fraction-min`         | off     | (nuclear fraction; disabled because it removes adipose) |

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

### 3. Integrate post-Xenium and standalone H&E (`04a_morphoNiches_integration.ipynb`)

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

### 4. Discover and annotate niches (`04b_morphoNiches_discovery.ipynb`)

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

### 5. Downstream analyses

- `04c_morphoNiches_sccoda.ipynb`: changes in niche abundance with age, tested with Milo and scCODA.
- `04d_morphoNiches_plotting.ipynb`: manuscript figure panels.
- `05a_patchExpression.ipynb`: patch × gene matrix (`thyAging_CAM_morphoFM_mpp025_mf2_patchExpr.h5ad`) built from the
  Xenium `spatialdata` zarrs and the alignment matrices.
- `00_*_niche.ipynb`: neighbourhood niches of selected Xenium populations, cross-referenced to the morphological niches.

## Applying the niches to new H&E slides

1. Run steps 1–2 on the new slides with the same parameters (0.25 mpp, 224 px tiles, `--grayscale`).
2. L2-normalise the grayscale embeddings and project them into the PCA space fitted in step 3.
3. Transfer the niche labels to the new tiles, for example by kNN label transfer in `X_pca_coral` space.

Standalone H&E slides are the CORAL *target*, so new non-Xenium slides do not need correcting. New post-Xenium slides should go through the same CORAL map (`A`, `μ_s`, `μ_t`).

## Data availability

The raw H&E images, Xenium data and processed objects are available as described in the manuscript's *Data
availability* section. The paths in the slide lists and notebooks point to our compute environment and must be changed
to your local copies.

## Citation

If you use ST2HE, please cite the manuscript *A Human Thymus Ageing Cell Atlas Defines A Spatially Resolved Model of Involution*, together with UNI (Chen et al., *Nat. Med.* 2024), CORAL (Sun et al., *AAAI* 2016) and bbknn
(Polański et al., *Bioinformatics* 2020).
