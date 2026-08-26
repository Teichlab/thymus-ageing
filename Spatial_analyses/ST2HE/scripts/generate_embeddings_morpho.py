from pathlib import Path
import warnings
import sys
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import torch
from huggingface_hub import login
from PIL import Image, ImageDraw, ImageFont
import timm
from timm.data import resolve_data_config
from timm.data.transforms_factory import create_transform
from timm.layers import SwiGLUPacked
import anndata as ad

from stain_norm import MacenkoNormalizer

warnings.filterwarnings("ignore")


def extract_patches_generator_pil(image_array, coordinates, crop_size=128, resize_to=224, pad_value=255,):
    """ Yield (oid, PIL.Image patch) from image_array using local image pixel coordinates. """
    half = crop_size // 2
    H, W = image_array.shape[:2]

    for oid, (x, y) in coordinates.items():
        x = int(round(x))
        y = int(round(y))

        left = x - half
        top = y - half
        right = left + crop_size
        bottom = top + crop_size

        src_left = max(0, left)
        src_top = max(0, top)
        src_right = min(W, right)
        src_bottom = min(H, bottom)

        dst_left = src_left - left
        dst_top = src_top - top
        dst_right = dst_left + (src_right - src_left)
        dst_bottom = dst_top + (src_bottom - src_top)

        patch = np.full((crop_size, crop_size, 3), pad_value, dtype=np.uint8)

        if src_right > src_left and src_bottom > src_top:
            patch[dst_top:dst_bottom, dst_left:dst_right] = image_array[src_top:src_bottom, src_left:src_right]

        patch_img = Image.fromarray(patch).resize((resize_to, resize_to), Image.Resampling.LANCZOS,)
        yield oid, patch_img
        
# def extract_tile_patches_generator_pil(image_array, tiles_df, resize_to=224,):
#     """ Yield (tile_id, PIL.Image patch) from image_array using tile coordinates. tiles_df requires: tile_id, x0, y0, x1, y1 """
#     for _, row in tiles_df.iterrows():
#         tile_id = row["tile_id"]

#         x0 = int(row["x0"])
#         y0 = int(row["y0"])
#         x1 = int(row["x1"])
#         y1 = int(row["y1"])

#         patch = image_array[y0:y1, x0:x1]
#         patch_img = Image.fromarray(patch).resize((resize_to, resize_to), Image.Resampling.LANCZOS,)

#         yield tile_id, patch_img
def extract_tile_patches_generator_pil(image_array, tiles_df, resize_to=224,):
    """
    Yield (tile_id, PIL.Image patch) from image_array using tile coordinates.
    Works with normal numpy array or np.load(..., mmap_mode="r").

    tiles_df requires: tile_id, x0, y0, x1, y1
    """
    for _, row in tiles_df.iterrows():
        tile_id = row["tile_id"]

        x0 = int(row["x0"])
        y0 = int(row["y0"])
        x1 = int(row["x1"])
        y1 = int(row["y1"])

        # important for mmap: convert only this patch to a normal ndarray
        patch = np.asarray(image_array[y0:y1, x0:x1])

        if patch.shape[0] == 0 or patch.shape[1] == 0:
            continue

        if patch.dtype != np.uint8:
            patch = np.clip(patch, 0, 255).astype(np.uint8)

        if patch.ndim != 3 or patch.shape[-1] != 3:
            raise ValueError(f"Expected RGB patch, got shape={patch.shape}")

        patch_img = Image.fromarray(patch)

        if resize_to is not None and patch_img.size != (resize_to, resize_to):
            patch_img = patch_img.resize((resize_to, resize_to), Image.Resampling.LANCZOS,)

        yield tile_id, patch_img


def to_grayscale_pil(patch_img):
    """
    Convert an H&E patch to luminance grayscale, kept as a 3-channel RGB image
    (UNI2-h expects 3 channels). Removes all colour so the embedding depends on
    morphology/texture alone, not on the stain/scanner colour axis.
    """
    return patch_img.convert("L").convert("RGB")


def build_uni2h(device):
    """
    Load UNI2-h (MahmoodLab/UNI2-h) via timm + hf-hub.
    Note: Don't hardcode your HF token in code. Use login() or env var HUGGINGFACE_HUB_TOKEN.
    """
    print("Preparing UNI2-h...")
    print(f"torch.cuda.is_available(): {torch.cuda.is_available()}")
    print(f"Using device: {device}")
    if device.type == "cuda":
        print("GPU:", torch.cuda.get_device_name(0))

    login()

    timm_kwargs = {
        'img_size': 224,
        'patch_size': 14,
        'depth': 24,
        'num_heads': 24,
        'init_values': 1e-5,
        'embed_dim': 1536,
        'mlp_ratio': 2.66667 * 2,
        'num_classes': 0,
        'no_embed_class': True,
        'mlp_layer': timm.layers.SwiGLUPacked,
        'act_layer': torch.nn.SiLU,
        'reg_tokens': 8,
        'dynamic_img_size': True}

    model = timm.create_model("hf-hub:MahmoodLab/UNI2-h", pretrained=True, **timm_kwargs)
    transform = create_transform(**resolve_data_config(model.pretrained_cfg, model=model))
    model = model.to(device).eval()
    print("UNI2-h loaded.")
    return model, transform

@torch.inference_mode()
def uni_embed_from_patch_generator(patch_generator, model, transform, device, batch_size=32, total=None,
                                   patch_transform=None, desc="UNI embedding"):
    """
    Stream patches -> (optional per-patch transform) -> timm transform -> batch
    -> UNI -> embeddings. Returns ids: list; embs: np.ndarray, shape (N,1536).

    patch_transform, if given, is a callable PIL.Image -> PIL.Image applied to
    every patch before the encoder (e.g. Macenko stain normalization, or
    grayscale conversion). It only changes pixels, not tile identity, so ids are
    identical across transforms for the same tiles_df.
    """
    ids_buf = []
    x_buf = []
    out_ids = []
    out_embs = []

    iterator = patch_generator
    if total is not None:
        iterator = tqdm(patch_generator, total=total, desc=desc)

    # bf16 autocast on CUDA: ~2x faster forward and ~half the activation memory,
    # with no meaningful change to the embeddings (inference only). No-op on CPU.
    use_amp = (device.type == "cuda")

    def _forward(x_list):
        x = torch.stack(x_list, dim=0).to(device, non_blocking=True)  # [B,3,224,224]
        with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=use_amp):
            y = model(x)  # [B,1536]
        return y.detach().float().cpu().numpy().astype(np.float32)

    for oid, patch_img in iterator:
        if patch_transform is not None:
            patch_img = patch_transform(patch_img)
        x_buf.append(transform(patch_img))  # tensor [3,224,224]
        ids_buf.append(oid)

        if len(x_buf) >= batch_size:
            out_embs.append(_forward(x_buf))
            out_ids.extend(ids_buf)
            x_buf.clear()
            ids_buf.clear()

    # last batch
    if len(x_buf) > 0:
        out_embs.append(_forward(x_buf))
        out_ids.extend(ids_buf)

    embs = np.concatenate(out_embs, axis=0) if len(out_embs) > 0 else np.zeros((0, 1536), dtype=np.float32)
    return out_ids, embs

def build_tile_adata(tiles_df, ids, emb_dict):
    """
    Build one AnnData for a set of tiles where each transform's embeddings are a
    separate obsm key. All transforms share identical tiles/obs, so we store the
    obs once and attach every embedding matrix as its own obsm entry (e.g.
    'mf2_tile_uni2h', 'mf2_tile_uni2h_stainMacenko', 'mf2_tile_uni2h_gray').

    ids must be shared across all emb_dict entries (same tile order).
    """
    ids = np.asarray(ids).astype(str)

    obs = tiles_df.copy()
    obs["tile_id"] = obs["tile_id"].astype(str)
    obs = obs.set_index("tile_id").loc[ids].copy()
    obs.index.name = "tile_id"

    adata = ad.AnnData(obs=obs)
    for emb_key, embs in emb_dict.items():
        adata.obsm[emb_key] = np.asarray(embs, dtype=np.float32)

    adata.uns["tile_info"] = {
        "emb_keys": list(emb_dict.keys()),
        "n_tiles": int(adata.n_obs),
        "embedding_dims": {k: int(np.asarray(v).shape[1]) for k, v in emb_dict.items()},
    }
    return adata


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Generate UNI2-h tile embeddings.")
    parser.add_argument(
        "--index",
        type=int,
        default=None,
        help="1-based index of a single sample to embed (e.g. LSF $LSB_JOBINDEX). "
             "Keyed to the same sorted input list as preprocess_morpho.py, so index N "
             "is the same slide in both steps. If omitted, all processed dirs are embedded.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("/lustre/scratch126/cellgen/teichmann/lm25/data/thymus_he/thyAgeing_HnE"),
        help="Directory holding the raw slides and the <sample>_processed output dirs.",
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="*.ndpi",
        help="Glob used ONLY in --index mode to map an index to a sample name. Must "
             "match the --pattern used in preprocess_morpho.py. Ignored with --file-list.",
    )
    parser.add_argument(
        "--file-list",
        type=Path,
        default=None,
        help="Text file of input filenames/paths used ONLY in --index mode to map an "
             "index to a sample name. Must match the --file-list used in preprocessing.",
    )
    parser.add_argument(
        "--stain-norm",
        action="store_true",
        help="Also compute Macenko H&E stain-normalized embeddings (obsm key "
             "'mf2_tile_uni2h_stainMacenko'), stored alongside the raw embeddings in "
             "the same h5ad. Removes the cross-platform (Xenium vs non-Xenium) "
             "stain/scanner colour axis.",
    )
    parser.add_argument(
        "--grayscale",
        action="store_true",
        help="Also compute grayscale (colour-removed) embeddings (obsm key "
             "'mf2_tile_uni2h_gray'), stored alongside the raw embeddings in the same "
             "h5ad. Tests whether morphology alone (no stain colour) integrates.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=128,
        help="UNI2-h forward batch size. Raise on a larger GPU (bf16 autocast keeps "
             "memory low), e.g. 384-512 on a full A100, for higher throughput.",
    )
    args = parser.parse_args()

    base_dir = args.data_dir
    # base_dir = Path("/rds/user/hx305/hpc-work/MorphoFM_data/ThymusAging/Goettingen")

    encoder = "uni2h"
    mpp_tag = "mpp025"
    tile_size = 224
    stride = 224
    resize_to = 224
    batch_size = args.batch_size

    # One h5ad per sample holds every requested transform as a separate obsm key
    # (the tiles/obs are identical across transforms). Raw is always computed;
    # --stain-norm and --grayscale each add another obsm key. Requested keys are
    # resolved up front so the idempotent skip can check they all already exist.
    RAW_KEY = "mf2_tile_uni2h"
    MACENKO_KEY = "mf2_tile_uni2h_stainMacenko"
    GRAY_KEY = "mf2_tile_uni2h_gray"
    requested_keys = [RAW_KEY]
    if args.stain_norm:
        requested_keys.append(MACENKO_KEY)
    if args.grayscale:
        requested_keys.append(GRAY_KEY)
    print(f"Embedding variants (obsm keys): {requested_keys}")

    if args.index is not None:
        # Array mode: resolve the sample name from the same source of truth as
        # preprocessing (same --pattern / --file-list selection), then target its
        # _processed dir. This keeps index N pointing at the same slide across steps
        # regardless of input format (.ndpi, .ome.tif, or a mixed file list).
        from preprocess_morpho import get_tiff_samples, get_samples_from_list

        if args.file_list is not None:
            samples = get_samples_from_list(args.file_list, base_dir=base_dir)
        else:
            samples = get_tiff_samples(base_dir, pattern=args.pattern)
        if not (1 <= args.index <= len(samples)):
            raise IndexError(f"--index {args.index} out of range 1..{len(samples)}")
        sample_name = samples[args.index - 1]["sample_name"]
        sample_dir = base_dir / f"{sample_name}_processed"
        print(f"Array mode: index {args.index} -> sample {sample_name}")

        if not sample_dir.is_dir():
            # Slide was skipped in preprocessing (e.g. coarse mpp) -> nothing to embed.
            print(f"[SKIP] No processed dir for {sample_name} at {sample_dir}; "
                  f"nothing to embed. Exiting cleanly.")
            return
        sample_dirs = [sample_dir]
    else:
        # Find processed sample folders, e.g. A43-THY-0-FFPE-1_processed
        sample_dirs = sorted([p for p in base_dir.glob("*_processed") if p.is_dir()])
    print(f"Found {len(sample_dirs)} sample dirs")

    # Idempotent skip: drop samples whose combined h5ad already holds ALL requested
    # obsm keys, before loading the model, so a fully-done task pays no model-load /
    # HF-login cost. If the file exists but is missing a newly-requested variant
    # (e.g. adding --grayscale to an earlier raw+Macenko run), it is recomputed and
    # overwritten with the full set.
    def existing_obsm_keys(path: Path) -> set:
        import h5py
        try:
            with h5py.File(path, "r") as f:
                return set(f["obsm"].keys()) if "obsm" in f else set()
        except Exception:
            return set()

    pending_dirs = []
    for d in sample_dirs:
        sname = d.name.removesuffix("_processed")
        combined_path = d / f"{sname}_THY_tile_all_parts_{mpp_tag}_uni2h_morpho.h5ad"
        if combined_path.exists() and set(requested_keys).issubset(existing_obsm_keys(combined_path)):
            print(f"[SKIP] Combined embeddings with all requested variants already exist for {sname}; skipping.")
        else:
            pending_dirs.append(d)
    sample_dirs = pending_dirs

    if len(sample_dirs) == 0:
        print("Nothing to embed; all outputs already exist.")
        return

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    if device.type == "cuda":
        print("GPU:", torch.cuda.get_device_name(0))

    model, transform = build_uni2h(device)

    # Map each requested obsm key to its per-patch transform (None = raw pixels).
    # Instantiated once and reused for every sample/part.
    variant_transforms = {RAW_KEY: None}
    if args.stain_norm:
        macenko = MacenkoNormalizer()
        variant_transforms[MACENKO_KEY] = macenko.normalize_pil
        print("Stain normalization: Macenko -> canonical H&E target (ON)")
    if args.grayscale:
        variant_transforms[GRAY_KEY] = to_grayscale_pil
        print("Grayscale (colour-removed) embeddings (ON)")

    for data_path in sample_dirs:
        # Recover the sample name used in the tile/image filenames (dir is "<sample>_processed")
        sample_name = data_path.name.removesuffix("_processed")

        print("\n" + "=" * 100)
        print(f"Processing sample: {sample_name}")
        print(f"Data path: {data_path}")

        tile_files = sorted(data_path.glob(f"{sample_name}_part*_he_tiles_{mpp_tag}_tile{tile_size}_stride{stride}.parquet"))
        print(f"Found {len(tile_files)} part tile files")

        all_part_adatas = []
        for tiles_path in tile_files:
            print("\n" + "-" * 80)
            print(f"Tiles file: {tiles_path.name}")

            tiles_df = pd.read_parquet(tiles_path)
            print("Number of tiles:", len(tiles_df))

            if len(tiles_df) == 0:
                print("[WARN] Empty tiles_df. Skip.")
                continue

            # Find matching image file.
            # Prefer source_image column from preprocessing.
            if "source_image" in tiles_df.columns:
                image_name = str(tiles_df["source_image"].iloc[0])
            else:
                image_name = tiles_path.name.replace(f"_he_tiles_{mpp_tag}_tile{tile_size}_stride{stride}.parquet", f"_image_he_{mpp_tag}.npy",)

            img_path = data_path / image_name

            if not img_path.exists():
                raise FileNotFoundError(f"Image file not found: {img_path}")

            if "part_name" in tiles_df.columns:
                part_name = str(tiles_df["part_name"].iloc[0])
            else:
                part_name = image_name.replace(f"_image_he_{mpp_tag}.npy", "")

            print(f"Part name: {part_name}")
            print(f"Image file: {img_path.name}")

            # Important: use mmap to avoid loading huge .npy fully into RAM
            img = np.load(img_path, mmap_mode="r")
            print("image shape:", img.shape)
            print("image dtype:", img.dtype)
            print(f"image array size: {img.nbytes / 1024**3:.2f} GiB")
            print(tiles_df.head())

            n_tiles = len(tiles_df)

            # Encode each requested variant over the SAME tiles. The generator is
            # deterministic, so ids match across variants; we assert it and store
            # every variant's embeddings as a separate obsm key on one adata.
            emb_dict = {}
            shared_ids = None
            for emb_key, ptf in variant_transforms.items():
                patch_gen = extract_tile_patches_generator_pil(image_array=img, tiles_df=tiles_df, resize_to=resize_to,)
                ids, embs = uni_embed_from_patch_generator(
                    patch_generator=patch_gen, model=model, transform=transform,
                    device=device, batch_size=batch_size, total=n_tiles,
                    patch_transform=ptf, desc=f"UNI embedding [{emb_key}]",
                )
                print(f"  {emb_key}: {embs.shape}")
                if shared_ids is None:
                    shared_ids = ids
                elif ids != shared_ids:
                    raise RuntimeError(
                        f"Tile ids differ between variants for {part_name} "
                        f"({emb_key}); cannot share one obs."
                    )
                emb_dict[emb_key] = embs

            adata = build_tile_adata(tiles_df=tiles_df, ids=shared_ids, emb_dict=emb_dict)
            print(adata)
            all_part_adatas.append(adata)

            # Release memory
            del img, tiles_df, shared_ids, emb_dict, adata
            if device.type == "cuda":
                torch.cuda.empty_cache()

        # Combine all parts for this sample into one h5ad holding every variant.
        if len(all_part_adatas) > 0:
            combined = ad.concat(all_part_adatas, axis=0, join="outer", merge="same", label="part_batch", index_unique=None,)
            combined_path = data_path / f"{sample_name}_THY_tile_all_parts_{mpp_tag}_uni2h_morpho.h5ad"
            combined.write_h5ad(combined_path, compression="gzip", compression_opts=4,)
            print("Saved combined sample h5ad:", combined_path)
            print("  obsm keys:", list(combined.obsm.keys()))

            del combined, all_part_adatas
            if device.type == "cuda":
                torch.cuda.empty_cache()

    print("All done.")



if __name__ == "__main__":
    main()
    
    
    