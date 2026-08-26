#!/usr/bin/env python
"""
Diagnostic for tile-level tissue filtering.

Purpose
-------
The production filter in preprocess_morpho.generate_he_tiles is deliberately
strict, which drops tissue edges, interstitial space and adipose (all pale, few
nuclei, low saturation). Before loosening those thresholds we want to *see* what
each one costs on real patches.

This script does NOT change preprocessing. It:
  1. loads the already-saved full-res image (<sample>_part*_image_he_mpp025.npy)
     so there is no expensive WSI re-read,
  2. tiles it with NO filtering and computes the full metric dict per tile,
  3. writes a per-tile metrics parquet,
  4. renders grid metric maps, per-metric histograms with the current threshold,
     a pass/fail map for the current filter stack, a filter drop-count summary,
     and patch galleries of the tiles each filter removes ("what am I losing")
     plus the borderline tiles straddling each threshold.

You inspect the output, choose thresholds, and then we bake the chosen values
into preprocess_morpho.py as argparse parameters.

Examples
--------
  # by explicit sample name
  python diagnose_tile_filtering.py --sample-name A43-THY-0-FFPE-1

  # by array index, keyed to the same list as preprocess/embeddings
  python diagnose_tile_filtering.py --index 3 --file-list nonXenium_slides.txt
"""
import argparse
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image

# Reuse the exact metric definitions the production filter uses, so the
# diagnostic reflects real behaviour rather than a re-implementation.
from preprocess_morpho import (
    compute_tissue_metrics,
    get_tiff_samples,
    get_samples_from_list,
    save_random_tile_gallery,
    save_image_preview,
)

warnings.filterwarnings("ignore")


# --- current production thresholds (from generate_he_tiles), all overridable ---
# A filter's `op` is the condition under which a tile is REMOVED.
DEFAULT_FILTERS = {
    "valid_tissue_fraction": {"thr": 0.15, "op": "lt", "label": "valid_tissue_fraction < thr"},
    "red_fraction":          {"thr": 0.05, "op": "gt", "label": "red_fraction > thr"},
    "saturation_fraction":   {"thr": 0.01, "op": "lt", "label": "saturation_fraction < thr"},
    "gray_std":              {"thr": 4.5,  "op": "lt", "label": "gray_std < thr"},
    # dark_fraction disabled as a filter (kept only as an informational metric map).
}
# Order in which the production code applies them (matters for cumulative counts).
FILTER_ORDER = ["valid_tissue_fraction", "red_fraction", "saturation_fraction", "gray_std"]


def resolve_sample_name(args) -> str:
    """ Resolve a single sample name from --sample-name, or --index against the same list source as the other scripts. """
    if args.sample_name is not None:
        return args.sample_name

    if args.index is None:
        raise SystemExit("Provide either --sample-name or --index.")

    if args.file_list is not None:
        samples = get_samples_from_list(args.file_list, base_dir=args.data_dir)
    else:
        samples = get_tiff_samples(args.data_dir, pattern=args.pattern)
    if not (1 <= args.index <= len(samples)):
        raise IndexError(f"--index {args.index} out of range 1..{len(samples)}")
    return samples[args.index - 1]["sample_name"]


def tile_all_with_metrics(img, tile_size, stride, sample, white_thr=230):
    """
    Tile `img` on a regular grid with NO filtering and return one row per tile
    with grid indices, box coords and the full compute_tissue_metrics dict.
    """
    H, W = img.shape[:2]
    n_rows = (H - tile_size) // stride + 1
    n_cols = (W - tile_size) // stride + 1

    records = []
    tile_idx = 0
    for gr, y0 in enumerate(range(0, H - tile_size + 1, stride)):
        for gc, x0 in enumerate(range(0, W - tile_size + 1, stride)):
            x1, y1 = x0 + tile_size, y0 + tile_size
            # mmap-friendly: materialise only this patch
            patch = np.asarray(img[y0:y1, x0:x1])
            m = compute_tissue_metrics(patch, white_thr=white_thr)
            records.append({
                "tile_id": f"{sample}_tile_{tile_idx:08d}",
                "sample": sample,
                "grid_row": gr,
                "grid_col": gc,
                "x0": int(x0), "y0": int(y0), "x1": int(x1), "y1": int(y1),
                "cx_img": float((x0 + x1) / 2), "cy_img": float((y0 + y1) / 2),
                **m,
            })
            tile_idx += 1

    df = pd.DataFrame(records)
    return df, (n_rows, n_cols)


def fail_mask(df, metric, spec):
    """ Boolean Series: True where the tile is REMOVED by this single filter. """
    if spec["op"] == "lt":
        return df[metric] < spec["thr"]
    return df[metric] > spec["thr"]


def save_metric_grid_maps(df, grid_shape, metrics, out_path):
    """ Dense heatmap per metric on the tile grid (every grid cell is a tile). """
    n_rows, n_cols = grid_shape
    n = len(metrics)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.2 * nrows), dpi=130)
    axes = np.array(axes).reshape(-1)
    for ax in axes:
        ax.axis("off")
    for ax, metric in zip(axes, metrics):
        grid = np.full((n_rows, n_cols), np.nan, dtype=np.float32)
        grid[df["grid_row"].values, df["grid_col"].values] = df[metric].values.astype(np.float32)
        im = ax.imshow(grid, cmap="viridis", interpolation="nearest")
        ax.set_title(metric, fontsize=10)
        ax.axis("on")
        ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved metric grid maps: {out_path}")


def save_metric_histograms(df, filters, out_path):
    """ Histogram per filtered metric with a vertical line at the current threshold. """
    metrics = list(filters.keys())
    n = len(metrics)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.6 * nrows), dpi=130)
    axes = np.array(axes).reshape(-1)
    for ax in axes:
        ax.axis("off")
    for ax, metric in zip(axes, metrics):
        spec = filters[metric]
        vals = df[metric].values.astype(float)
        ax.axis("on")
        ax.hist(vals, bins=60, color="#4C78A8")
        ax.axvline(spec["thr"], color="red", linestyle="--", linewidth=1.5,
                   label=f"thr={spec['thr']:g}")
        removed = int(fail_mask(df, metric, spec).sum())
        ax.set_title(f"{spec['label']}\nremoves {removed}/{len(df)} "
                     f"({100*removed/max(len(df),1):.1f}%)", fontsize=9)
        ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved metric histograms: {out_path}")


def save_passfail_map(df, grid_shape, filters, out_path):
    """
    Grid map coloured by outcome under the CURRENT filter stack:
      0 = removed, 1 = kept.
    """
    n_rows, n_cols = grid_shape
    removed = np.zeros(len(df), dtype=bool)
    for metric, spec in filters.items():
        removed |= fail_mask(df, metric, spec).values

    grid = np.full((n_rows, n_cols), np.nan, dtype=np.float32)
    grid[df["grid_row"].values, df["grid_col"].values] = (~removed).astype(np.float32)

    fig, ax = plt.subplots(figsize=(8, 8), dpi=140)
    ax.imshow(grid, cmap="RdYlGn", vmin=0, vmax=1, interpolation="nearest")
    ax.set_title(f"Current filter stack: kept (green) {int((~removed).sum())} / "
                 f"removed (red) {int(removed.sum())} of {len(df)}", fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    plt.tight_layout()
    plt.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved pass/fail map: {out_path}")
    return removed


def write_summary(df, filters, out_path):
    """ Per-filter and cumulative drop counts (cumulative follows FILTER_ORDER). """
    n = len(df)
    lines = [f"Total candidate tiles (no filtering): {n}", ""]

    lines.append("Per-filter removal (each filter in isolation):")
    fail_masks = {}
    for metric in FILTER_ORDER:
        spec = filters[metric]
        fm = fail_mask(df, metric, spec).values
        fail_masks[metric] = fm
        lines.append(f"  {spec['label']:<34} removes {fm.sum():>7} ({100*fm.sum()/max(n,1):5.1f}%)")

    lines.append("")
    lines.append("Tiles removed ONLY by this filter (would be rescued if it were dropped):")
    any_other = {m: np.zeros(n, dtype=bool) for m in FILTER_ORDER}
    for m in FILTER_ORDER:
        for m2 in FILTER_ORDER:
            if m2 != m:
                any_other[m] |= fail_masks[m2]
    for metric in FILTER_ORDER:
        only = fail_masks[metric] & (~any_other[metric])
        lines.append(f"  {filters[metric]['label']:<34} unique {only.sum():>7} ({100*only.sum()/max(n,1):5.1f}%)")

    lines.append("")
    lines.append("Cumulative kept as filters are applied in production order:")
    kept = np.ones(n, dtype=bool)
    lines.append(f"  {'(start)':<34} kept {kept.sum():>7}")
    for metric in FILTER_ORDER:
        kept &= ~fail_masks[metric]
        lines.append(f"  after {filters[metric]['label']:<28} kept {kept.sum():>7} ({100*kept.sum()/max(n,1):5.1f}%)")

    text = "\n".join(lines)
    Path(out_path).write_text(text + "\n")
    print("\n" + text + "\n")
    print(f"Saved summary: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Diagnose tile tissue-filtering thresholds (read-only; does not touch preprocessing).")
    parser.add_argument("--data-dir", type=Path,
                        default=Path("/lustre/scratch126/cellgen/teichmann/lm25/data/thymus_he/thyAgeing_HnE"))
    parser.add_argument("--sample-name", type=str, default=None,
                        help="Explicit sample name (its <sample>_processed dir must exist).")
    parser.add_argument("--index", type=int, default=None,
                        help="1-based index into --file-list/--pattern (same source as the other scripts).")
    parser.add_argument("--pattern", type=str, default="*.ndpi")
    parser.add_argument("--file-list", type=Path, default=None)
    parser.add_argument("--tile-size", type=int, default=224)
    parser.add_argument("--stride", type=int, default=224)
    parser.add_argument("--mpp-tag", type=str, default="mpp025")
    parser.add_argument("--white-thr", type=int, default=230)
    parser.add_argument("--n-gallery", type=int, default=25, help="Patches per gallery.")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out-subdir", type=str, default="filter_diagnostics")
    # Threshold overrides (default to current production values).
    for metric, spec in DEFAULT_FILTERS.items():
        parser.add_argument(f"--thr-{metric.replace('_', '-')}", type=float, default=spec["thr"],
                            help=f"Threshold for {spec['label']} (default {spec['thr']}).")
    args = parser.parse_args()

    filters = {m: {**s, "thr": getattr(args, f"thr_{m}")} for m, s in DEFAULT_FILTERS.items()}

    sample_name = resolve_sample_name(args)
    sample_dir = args.data_dir / f"{sample_name}_processed"
    if not sample_dir.is_dir():
        raise FileNotFoundError(
            f"No processed dir for {sample_name} at {sample_dir}. "
            f"Run preprocess_morpho.py first (it saves the image .npy regardless of filtering)."
        )

    out_dir = sample_dir / args.out_subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    image_files = sorted(sample_dir.glob(f"{sample_name}_part*_image_he_{args.mpp_tag}.npy"))
    if not image_files:
        raise FileNotFoundError(
            f"No image .npy matching {sample_name}_part*_image_he_{args.mpp_tag}.npy in {sample_dir}"
        )

    print(f"Sample: {sample_name}")
    print(f"Output dir: {out_dir}")
    print(f"Found {len(image_files)} part image(s): {[p.name for p in image_files]}")

    all_metrics = []
    for img_path in image_files:
        part_name = img_path.name.replace(f"_image_he_{args.mpp_tag}.npy", "")
        print("\n" + "=" * 80)
        print(f"Part: {part_name}  ({img_path.name})")

        # mmap: avoid loading the multi-GB slide fully into RAM.
        img = np.load(img_path, mmap_mode="r")
        print(f"  image shape: {img.shape}, dtype: {img.dtype}, "
              f"size: {img.nbytes / 1024**3:.2f} GiB")

        df, grid_shape = tile_all_with_metrics(
            img, tile_size=args.tile_size, stride=args.stride,
            sample=sample_name, white_thr=args.white_thr,
        )
        df["part_name"] = part_name
        print(f"  candidate tiles (no filtering): {len(df)}  grid {grid_shape}")

        prefix = f"{part_name}_{args.mpp_tag}"

        # --- reference thumbnail so grid maps can be read against real tissue ---
        save_image_preview(np.asarray(img), out_dir / f"{prefix}_thumbnail.png", max_dim=3000)

        # --- tabular + plots ---
        df.to_parquet(out_dir / f"{prefix}_tile_metrics_all.parquet", index=False)
        write_summary(df, filters, out_dir / f"{prefix}_filter_summary.txt")
        save_metric_histograms(df, filters, out_dir / f"{prefix}_metric_histograms.png")
        save_metric_grid_maps(
            df, grid_shape,
            metrics=["valid_tissue_fraction", "tissue_fraction", "saturation_fraction",
                     "dark_fraction", "gray_std", "red_fraction"],
            out_path=out_dir / f"{prefix}_metric_grid_maps.png",
        )
        removed = save_passfail_map(df, grid_shape, filters, out_dir / f"{prefix}_passfail_map.png")

        # --- galleries: what each filter is throwing away ---
        img_full = np.asarray(img)  # galleries need real pixels; materialise once

        kept_df = df.loc[~removed]
        if len(kept_df):
            save_random_tile_gallery(
                img_full, kept_df, out_path=out_dir / f"{prefix}_gallery_KEPT.png",
                n_show=args.n_gallery, seed=args.seed, extra_col="valid_tissue_fraction", show=False,
            )

        # per-filter "removed only by this filter" = tiles rescued if we dropped it
        fail_masks = {m: fail_mask(df, m, filters[m]).values for m in FILTER_ORDER}
        for metric in FILTER_ORDER:
            others = np.zeros(len(df), dtype=bool)
            for m2 in FILTER_ORDER:
                if m2 != metric:
                    others |= fail_masks[m2]
            only = fail_masks[metric] & (~others)
            sub = df.loc[only]
            if len(sub):
                save_random_tile_gallery(
                    img_full, sub, out_path=out_dir / f"{prefix}_gallery_removed_only_by_{metric}.png",
                    n_show=args.n_gallery, seed=args.seed, extra_col=metric, show=False,
                )
            else:
                print(f"  (no tiles removed uniquely by {metric})")

            # borderline tiles straddling the threshold (±15% of thr band)
            thr = filters[metric]["thr"]
            band = max(abs(thr) * 0.15, 1e-3)
            near = df.loc[(df[metric] >= thr - band) & (df[metric] <= thr + band)]
            if len(near):
                save_random_tile_gallery(
                    img_full, near, out_path=out_dir / f"{prefix}_gallery_borderline_{metric}.png",
                    n_show=args.n_gallery, seed=args.seed, extra_col=metric, show=False,
                )

        all_metrics.append(df)
        del img, img_full

    if len(all_metrics) > 1:
        combined = pd.concat(all_metrics, ignore_index=True)
        combined.to_parquet(out_dir / f"{sample_name}_{args.mpp_tag}_tile_metrics_all_parts.parquet", index=False)
        print(f"\nSaved combined metrics for {len(all_metrics)} parts.")

    print("\nDiagnostics done. Inspect:", out_dir)
    print("Key files: *_filter_summary.txt, *_metric_histograms.png, *_passfail_map.png, "
          "*_gallery_removed_only_by_*.png")


if __name__ == "__main__":
    main()
