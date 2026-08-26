"""
CPU-only sanity check for Macenko stain normalization, before committing a full
GPU re-embed.

Samples tiles from Xenium and non-Xenium *_processed dirs and reports the
per-platform mean RGB (and OD) before vs after normalization. If the dominant
UNI2-h batch axis is stain/scanner colour, the between-platform RGB distance
should shrink markedly after normalization. Also writes a before/after montage.

Xenium vs non-Xenium is inferred by matching each processed dir name against the
sample prefixes in xenium_slides.txt.

Usage
-----
    python validate_stain_norm.py \
        --data-dir /lustre/scratch126/cellgen/teichmann/lm25/data/thymus_he/thyAgeing_HnE \
        --xenium-list ./xenium_slides.txt \
        --n-per-slide 200 --n-slides 3 --out ./_logs/stain_norm_check.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt

from stain_norm import MacenkoNormalizer
from generate_embeddings_morpho import extract_tile_patches_generator_pil

MPP_TAG = "mpp025"
TILE = 224
STRIDE = 224


def xenium_prefixes(list_path: Path) -> set[str]:
    """Sample-name stems from the Xenium slide list (e.g. A67, Z16)."""
    prefixes = set()
    for line in Path(list_path).read_text().splitlines():
        line = line.strip()
        if line:
            prefixes.add(Path(line).name.split(".")[0].split("-")[0])
    return prefixes


def is_xenium(dir_name: str, prefixes: set[str]) -> bool:
    stem = dir_name.removesuffix("_processed")
    return stem.split("-")[0] in prefixes


def load_first_part(sample_dir: Path):
    """Return (tiles_df, mmap image) for the first tile-part of a processed dir."""
    sample = sample_dir.name.removesuffix("_processed")
    parts = sorted(sample_dir.glob(
        f"{sample}_part*_he_tiles_{MPP_TAG}_tile{TILE}_stride{STRIDE}.parquet"))
    if not parts:
        return None, None
    tiles_df = pd.read_parquet(parts[0])
    if len(tiles_df) == 0:
        return None, None
    if "source_image" in tiles_df.columns:
        image_name = str(tiles_df["source_image"].iloc[0])
    else:
        image_name = parts[0].name.replace(
            f"_he_tiles_{MPP_TAG}_tile{TILE}_stride{STRIDE}.parquet",
            f"_image_he_{MPP_TAG}.npy")
    img_path = sample_dir / image_name
    if not img_path.exists():
        return None, None
    return tiles_df, np.load(img_path, mmap_mode="r")


def sample_patches(sample_dir: Path, n: int, rng) -> list[np.ndarray]:
    tiles_df, img = load_first_part(sample_dir)
    if tiles_df is None:
        return []
    if len(tiles_df) > n:
        tiles_df = tiles_df.iloc[rng.choice(len(tiles_df), n, replace=False)]
    return [np.asarray(p.convert("RGB"))
            for _, p in extract_tile_patches_generator_pil(img, tiles_df, resize_to=TILE)]


def group_stats(patches, norm):
    """Mean RGB across a list of tiles, before and after normalization."""
    orig = np.stack([p.reshape(-1, 3).mean(0) for p in patches])
    normed = np.stack([norm.normalize(p).reshape(-1, 3).mean(0) for p in patches])
    return orig.mean(0), normed.mean(0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=Path, required=True)
    ap.add_argument("--xenium-list", type=Path,
                    default=Path(__file__).with_name("xenium_slides.txt"))
    ap.add_argument("--n-per-slide", type=int, default=200)
    ap.add_argument("--n-slides", type=int, default=3,
                    help="Max processed dirs to sample per platform.")
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).parent / "_logs" / "stain_norm_check.png")
    args = ap.parse_args()

    rng = np.random.default_rng(0)
    norm = MacenkoNormalizer()
    prefixes = xenium_prefixes(args.xenium_list)

    dirs = sorted(p for p in args.data_dir.glob("*_processed") if p.is_dir())
    xen = [d for d in dirs if is_xenium(d.name, prefixes)][: args.n_slides]
    oth = [d for d in dirs if not is_xenium(d.name, prefixes)][: args.n_slides]
    print(f"Xenium dirs: {[d.name for d in xen]}")
    print(f"Non-Xenium dirs: {[d.name for d in oth]}")

    groups = {"Xenium": xen, "non-Xenium": oth}
    means = {}
    examples = {}
    for label, gdirs in groups.items():
        patches = []
        for d in gdirs:
            patches += sample_patches(d, args.n_per_slide, rng)
        if not patches:
            print(f"[WARN] No patches found for {label}; check paths.")
            continue
        means[label] = group_stats(patches, norm)
        examples[label] = patches[:4]
        print(f"\n[{label}] n={len(patches)} tiles")
        print(f"  mean RGB  before: {means[label][0].round(1)}")
        print(f"  mean RGB  after : {means[label][1].round(1)}")

    if len(means) == 2:
        (xb, xa), (ob, oa) = means["Xenium"], means["non-Xenium"]
        d_before = float(np.linalg.norm(xb - ob))
        d_after = float(np.linalg.norm(xa - oa))
        print("\n=== Between-platform mean-RGB distance ===")
        print(f"  before norm: {d_before:7.2f}")
        print(f"  after  norm: {d_after:7.2f}   "
              f"({100 * (1 - d_after / max(d_before, 1e-9)):.0f}% reduction)")

    # Before/after montage.
    rows = sum(len(v) for v in examples.values())
    if rows:
        fig, axes = plt.subplots(rows, 2, figsize=(4, 2 * rows))
        axes = np.atleast_2d(axes)
        r = 0
        for label, pats in examples.items():
            for p in pats:
                axes[r, 0].imshow(p)
                axes[r, 0].set_ylabel(label, fontsize=8)
                axes[r, 1].imshow(norm.normalize(p))
                r += 1
        axes[0, 0].set_title("original")
        axes[0, 1].set_title("Macenko-normalized")
        for a in axes.ravel():
            a.set_xticks([]); a.set_yticks([])
        args.out.parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        fig.savefig(args.out, dpi=120)
        print(f"\nSaved montage: {args.out}")


if __name__ == "__main__":
    main()
