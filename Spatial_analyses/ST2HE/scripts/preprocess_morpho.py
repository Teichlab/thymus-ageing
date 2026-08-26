#!/usr/bin/env python
import re
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import pickle
from matplotlib.patches import Rectangle
from PIL import Image
import tifffile
import cv2
import openslide
import warnings
warnings.filterwarnings("ignore")


def extract_sample_name(file_name: str) -> str:
    """
    Derive a stable sample name from the WSI filename: drop the extension and,
    if present, the trailing ' - <timestamp>' scanner suffix. If there is no
    ' - ' timestamp, use the whole filename stem.

    Examples:
      'A43-THY-0-FFPE-1 - 2024-10-15 17.24.26.ndpi' -> 'A43-THY-0-FFPE-1'
      'S000000011_A16_Thymus.ndpi'                  -> 'S000000011_A16_Thymus'

    Matches test_morpho.py, which normalizes dir names the same way (split on ' - ').
    """
    stem = Path(file_name).stem

    if " - " in stem:
        sample_name = stem.split(" - ")[0].strip()
    else:
        sample_name = stem.strip()

    # Fall back to the full stem rather than failing (e.g. a leading ' - ').
    if not sample_name:
        sample_name = stem.strip()

    if not sample_name:
        raise ValueError(f"Could not derive sample name from: {file_name}")

    return sample_name


def get_tiff_samples(data_dir: str | Path, pattern: str = "*.tiff") -> list[dict]:
    """ Return all files matching `pattern` as a list of dicts. Each item: {"sample_name": "RE_Tym_000001", "file_name": "...", "file_path": Path(...) } """
    data_dir = Path(data_dir)

    if not data_dir.exists():
        raise FileNotFoundError(f"Directory does not exist: {data_dir}")

    files = sorted(data_dir.glob(pattern))
    samples = []
    for path in files:
        sample_name = extract_sample_name(path.name)
        samples.append({"sample_name": sample_name, "file_name": path.name, "file_path": path,})

    return samples


def get_samples_from_list(list_path: str | Path, base_dir: str | Path | None = None) -> list[dict]:
    """
    Build the sample list from a text file of filenames/paths, one per line.

    Blank lines and lines starting with '#' are ignored. Relative entries are
    resolved against `base_dir` when given. The input format of each sample is
    determined per-file from its suffix (see detect_format), so a single list may
    mix formats (e.g. some .ndpi and some .ome.tif). Each item is a dict:
    {"sample_name": ..., "file_name": ..., "file_path": Path(...)}.
    """
    list_path = Path(list_path)
    if not list_path.exists():
        raise FileNotFoundError(f"File list does not exist: {list_path}")

    samples = []
    for line in list_path.read_text().splitlines():
        entry = line.strip()
        if not entry or entry.startswith("#"):
            continue
        path = Path(entry)
        if base_dir is not None and not path.is_absolute():
            path = Path(base_dir) / path
        # Validate the suffix up front so bad entries fail fast with a clear message.
        detect_format(path)
        samples.append({
            "sample_name": extract_sample_name(path.name),
            "file_name": path.name,
            "file_path": path,
        })

    return samples


def read_wsi_mpp(wsi_path: str | Path) -> float | None:
    """
    Read only the native resolution (microns/pixel) from a WSI's metadata,
    without loading any pixels. Returns None if the slide has no mpp property.

    Used to decide whether to skip a slide before the expensive full read.
    """
    slide = openslide.OpenSlide(str(wsi_path))
    try:
        mpp_x = slide.properties.get(openslide.PROPERTY_NAME_MPP_X)
        return float(mpp_x) if mpp_x is not None else None
    finally:
        slide.close()


def read_wsi_as_rgb(wsi_path: str | Path, level: int = 0) -> tuple[np.ndarray, float | None]:
    """
    Read a whole-slide image (NDPI/SVS/TIFF, anything OpenSlide supports) as an
    RGB uint8 numpy array, and return its native resolution (microns/pixel).

    The full requested pyramid level is read into memory (level 0 = full res).
    NDPI is pyramidal; 0.227 -> 0.25 mpp is only a slight downscale, so level 0
    is the right level to start from.

    Returns (img_rgb, native_mpp). native_mpp is None if the slide has no mpp
    metadata, in which case the caller must supply orig_mpp explicitly.
    """
    wsi_path = Path(wsi_path)
    print(f"Reading WSI via OpenSlide: {wsi_path}")

    slide = openslide.OpenSlide(str(wsi_path))
    try:
        mpp_x = slide.properties.get(openslide.PROPERTY_NAME_MPP_X)
        mpp_y = slide.properties.get(openslide.PROPERTY_NAME_MPP_Y)
        native_mpp = float(mpp_x) if mpp_x is not None else None
        if mpp_x is not None and mpp_y is not None and abs(float(mpp_x) - float(mpp_y)) > 1e-3:
            print(f"[WARN] Non-square pixels: mpp_x={mpp_x}, mpp_y={mpp_y}. Using mpp_x.")

        W, H = slide.level_dimensions[level]
        print(f"  levels: {slide.level_count}, level {level} size (WxH): {W}x{H}")
        print(f"  native mpp: {native_mpp}")

        region = slide.read_region((0, 0), level, (W, H))  # PIL RGBA
        img = np.asarray(region.convert("RGB"))
    finally:
        slide.close()

    print("RGB image shape:", img.shape, "dtype:", img.dtype)
    return img, native_mpp


def to_rgb_uint8(img: np.ndarray) -> np.ndarray:
    """
    Normalize a raw TIFF/OME-TIFF array into an RGB uint8 image (H, W, 3).

    Handles the common layouts we see in H&E data:
      - HWC with 3/4 channels (drops alpha)
      - CHW with 3/4 channels (transposes)
      - 2D grayscale (broadcast to 3 channels)
      - 4D multi-page/pyramid arrays (takes the first image, then recurses)
    """
    img = np.asarray(img)

    # HWC: y, x, c
    if img.ndim == 3 and img.shape[-1] in (3, 4):
        img = img[..., :3]

    # CHW: c, y, x
    elif img.ndim == 3 and img.shape[0] in (3, 4):
        img = np.transpose(img[:3], (1, 2, 0))

    # 2D grayscale (e.g. a single-channel OME-TIFF plane)
    elif img.ndim == 2:
        img = np.stack([img] * 3, axis=-1)

    # Some TIFFs may be multi-page / pyramid-like when read directly
    elif img.ndim == 4:
        return to_rgb_uint8(img[0])

    else:
        raise ValueError(f"Unexpected image shape: {img.shape}")

    if img.dtype != np.uint8:
        img = np.clip(img, 0, 255).astype(np.uint8)

    return img


def read_he_tiff_as_rgb(he_path: str | Path) -> np.ndarray:
    """ Read H&E TIFF as RGB uint8 numpy array. Assumes the image can fit into memory. """
    he_path = Path(he_path)
    print(f"Reading TIFF: {he_path}")

    img = np.asarray(tifffile.imread(str(he_path)))

    print("Raw image shape:", img.shape)
    print("Raw image dtype:", img.dtype)

    img = to_rgb_uint8(img)

    print("RGB image shape:", img.shape)
    print("RGB image dtype:", img.dtype)

    return img


# microns per one unit of the named physical size (OME PhysicalSize units).
_MICRON_PER_UNIT = {
    "µm": 1.0, "um": 1.0, "micron": 1.0, "microns": 1.0,
    "micrometer": 1.0, "micrometre": 1.0,
    "nm": 1e-3, "mm": 1e3, "cm": 1e4, "m": 1e6,
}


def _mpp_from_resolution_tags(tif: "tifffile.TiffFile") -> float | None:
    """ Best-effort mpp from baseline TIFF XResolution/ResolutionUnit tags. """
    try:
        tags = tif.pages[0].tags
        xres = tags.get("XResolution")
        if xres is None:
            return None
        val = xres.value
        res = (val[0] / val[1]) if isinstance(val, tuple) else float(val)  # pixels per unit
        if not res:
            return None
        resunit = tags.get("ResolutionUnit")
        unit_code = int(resunit.value) if resunit is not None else 2
        if unit_code == 3:      # centimeter
            microns_per_unit = 1e4
        elif unit_code == 2:    # inch
            microns_per_unit = 25400.0
        else:                   # no absolute unit -> cannot derive mpp
            return None
        return microns_per_unit / res
    except Exception:
        return None


def extract_ome_mpp(tif: "tifffile.TiffFile") -> float | None:
    """
    Read microns/pixel from an OME-TIFF, preferring the OME-XML PhysicalSizeX
    (converted to microns via its declared unit), falling back to the baseline
    TIFF resolution tags. Returns None if no resolution metadata is present.
    """
    ome = getattr(tif, "ome_metadata", None)
    if ome:
        size = re.search(r'PhysicalSizeX\s*=\s*"([0-9.eE+-]+)"', ome)
        if size:
            mpp = float(size.group(1))
            unit_m = re.search(r'PhysicalSizeXUnit\s*=\s*"([^"]+)"', ome)
            unit = unit_m.group(1).strip() if unit_m else "µm"  # OME default is µm
            scale = _MICRON_PER_UNIT.get(unit)
            if scale is None:
                print(f"[WARN] Unknown OME PhysicalSizeXUnit '{unit}'; assuming microns.")
                scale = 1.0
            return mpp * scale
    return _mpp_from_resolution_tags(tif)


def read_ome_tiff_mpp(path: str | Path) -> float | None:
    """ Read only the native mpp from an OME-TIFF's metadata, without loading pixels. """
    with tifffile.TiffFile(str(path)) as tif:
        return extract_ome_mpp(tif)


def read_ome_tiff_as_rgb(path: str | Path) -> tuple[np.ndarray, float | None]:
    """
    Read an OME-TIFF as an RGB uint8 numpy array and return its native mpp.

    The full-resolution base level of the first image series is read into memory
    (pyramidal OME-TIFFs expose coarser levels, but we start from level 0 to match
    read_wsi_as_rgb). Returns (img_rgb, native_mpp); native_mpp is None if the file
    carries no resolution metadata, in which case the caller supplies orig_mpp.
    """
    path = Path(path)
    print(f"Reading OME-TIFF via tifffile: {path}")

    with tifffile.TiffFile(str(path)) as tif:
        native_mpp = extract_ome_mpp(tif)
        # Base (full-res) level of the first series.
        arr = tif.series[0].asarray()

    print("Raw image shape:", arr.shape, "dtype:", arr.dtype)
    print(f"  native mpp: {native_mpp}")

    img = to_rgb_uint8(arr)
    print("RGB image shape:", img.shape, "dtype:", img.dtype)
    return img, native_mpp


# Suffix -> reader dispatch. OME-TIFF uses a compound suffix and must be matched
# before the plain .tif/.tiff OpenSlide suffixes.
OME_TIFF_SUFFIXES = (".ome.tif", ".ome.tiff")
OPENSLIDE_SUFFIXES = (".ndpi", ".svs", ".tif", ".tiff", ".scn", ".mrxs", ".vms", ".vmu", ".bif")


def detect_format(path: str | Path) -> str:
    """ Return an input-format key ('ome-tiff' or 'openslide') from the filename suffix. """
    name = Path(path).name.lower()
    if name.endswith(OME_TIFF_SUFFIXES):
        return "ome-tiff"
    if name.endswith(OPENSLIDE_SUFFIXES):
        return "openslide"
    raise ValueError(
        f"Unsupported input format for file: {Path(path).name}. "
        f"Supported: {OME_TIFF_SUFFIXES + OPENSLIDE_SUFFIXES}"
    )


def read_native_mpp(path: str | Path) -> float | None:
    """ Format-aware, metadata-only native mpp read (no pixels). """
    if detect_format(path) == "ome-tiff":
        return read_ome_tiff_mpp(path)
    return read_wsi_mpp(path)


def read_image_as_rgb(path: str | Path, level: int = 0) -> tuple[np.ndarray, float | None]:
    """ Format-aware full read: returns (img_rgb_uint8, native_mpp). """
    if detect_format(path) == "ome-tiff":
        return read_ome_tiff_as_rgb(path)
    return read_wsi_as_rgb(path, level=level)


def resample_to_target_mpp(img: np.ndarray, orig_mpp: float, target_mpp: float, tol: float = 1e-4,) -> np.ndarray:
    """
    Resample an RGB image from its native resolution to a target resolution.

    mpp = microns per pixel, so a *larger* mpp means *lower* resolution.
    Scaling factor = orig_mpp / target_mpp:
      - orig_mpp < target_mpp  -> factor < 1 -> downscale (e.g. 0.227 -> 0.25 => 0.908)
      - orig_mpp > target_mpp  -> factor > 1 -> upscale

    Must be applied BEFORE tiling so that tile coordinates and physical sizes
    (tile_size * target_mpp) are correct.

    Uses cv2.INTER_AREA for downscaling (best when shrinking) and INTER_LINEAR
    for upscaling. Falls back to PIL for images whose target dimensions exceed
    OpenCV's 32767 px per-axis limit (common for whole-slide images).
    """
    if abs(orig_mpp - target_mpp) <= tol:
        print(f"orig_mpp ({orig_mpp}) ~= target_mpp ({target_mpp}); no resampling needed.")
        return img

    scale = orig_mpp / target_mpp
    H, W = img.shape[:2]
    new_W = int(round(W * scale))
    new_H = int(round(H * scale))

    print("Resampling to target mpp:")
    print(f"  orig_mpp={orig_mpp}, target_mpp={target_mpp}, scale={scale:.6f}")
    print(f"  {W}x{H} -> {new_W}x{new_H} (WxH)")

    if max(new_W, new_H, W, H) > 32767:
        # OpenCV resize is unreliable above 32767 px on an axis; use PIL instead.
        resample = Image.Resampling.LANCZOS if scale < 1.0 else Image.Resampling.BILINEAR
        img_resized = np.asarray(Image.fromarray(img).resize((new_W, new_H), resample))
    else:
        interp = cv2.INTER_AREA if scale < 1.0 else cv2.INTER_LINEAR
        img_resized = cv2.resize(img, (new_W, new_H), interpolation=interp)

    return np.ascontiguousarray(img_resized, dtype=np.uint8)


def compute_tissue_fraction(patch, white_thr=230):
    """ Simple tissue/background filter. patch: RGB uint8 image, shape (H, W, 3) """
    gray = patch.mean(axis=2)
    tissue_fraction = np.mean(gray < white_thr)
    return float(tissue_fraction)


# def compute_tissue_metrics(patch, white_thr=230, red_fraction_thr=0.5):
#     """ Tissue/background filter for H&E patches. Excludes white background and obvious red artifact/background regions. patch: RGB uint8 image, shape (H, W, 3) """
#     patch = patch.astype(np.uint8)

#     r = patch[..., 0].astype(float)
#     g = patch[..., 1].astype(float)
#     b = patch[..., 2].astype(float)

#     gray = patch.mean(axis=2)
#     # Non-white pixels
#     non_white = gray < white_thr
#     # Detect strong red pixels
#     red_mask = ((r > 150) & (r > g * 1.3) & (r > b * 1.3))
#     # Valid tissue = non-white but not red artifact
#     valid_tissue = non_white & (~red_mask)

#     tissue_fraction = float(np.mean(non_white))
#     red_fraction = float(np.mean(red_mask))
#     valid_tissue_fraction = float(np.mean(valid_tissue))

#     return tissue_fraction, red_fraction, valid_tissue_fraction

def compute_tissue_metrics(patch, white_thr=230, red_r_min=220, red_g_max=80, red_b_max=80, sat_thr=0.08, dark_thr=180,):
    """ More robust tissue/background filter for H&E patches. Filters out pale smooth background that may have tissue_fraction=1.0. """
    if patch.dtype != np.uint8:
        patch = patch.astype(np.uint8, copy=False)

    r = patch[..., 0].astype(np.uint16)
    g = patch[..., 1].astype(np.uint16)
    b = patch[..., 2].astype(np.uint16)

    gray_sum = r + g + b
    gray = gray_sum / 3.0

    # Non-white / non-blank pixels
    non_white = gray_sum < (3 * white_thr)

    # Strict pure-red artifact mask
    red_mask = ((r >= red_r_min) & (g <= red_g_max) & (b <= red_b_max))

    # Saturation-like metric
    maxc = np.maximum(np.maximum(r, g), b).astype(np.float32)
    minc = np.minimum(np.minimum(r, g), b).astype(np.float32)
    saturation = (maxc - minc) / (maxc + 1e-6)

    # Stained pixels: colored enough and not red artifact
    stained = (saturation > sat_thr) & (~red_mask)

    # Dark pixels: likely nuclei / dense tissue
    dark = (gray < dark_thr) & (~red_mask)

    valid_tissue = non_white & (~red_mask)

    tissue_fraction = float(non_white.mean())
    red_fraction = float(red_mask.mean())
    valid_tissue_fraction = float(valid_tissue.mean())

    saturation_fraction = float(stained.mean())
    dark_fraction = float(dark.mean())
    gray_std = float(gray.std())

    return {
        "tissue_fraction": tissue_fraction,
        "red_fraction": red_fraction,
        "valid_tissue_fraction": valid_tissue_fraction,
        "saturation_fraction": saturation_fraction,
        "dark_fraction": dark_fraction,
        "gray_std": gray_std,
    }
    

def generate_he_tiles(
    img,
    tile_size=224,
    stride=224,
    valid_tissue_fraction_min=0.15,
    red_fraction_max=0.05,
    saturation_fraction_min=0.01,
    gray_std_min=4.5,
    dark_fraction_min=None,
    white_thr=230,
    sample=None,
):
    """
    Generate regular grid tiles from aligned H&E image.

    Tile-keep thresholds (tuned via diagnose_tile_filtering.py to keep tissue
    edges / interstitial space / adipose while still rejecting empty background):
      - valid_tissue_fraction_min: min non-white, non-red fraction to keep.
      - red_fraction_max:          drop tiles with more red artifact than this.
      - saturation_fraction_min:   min fraction of colour-saturated (stained) pixels.
      - gray_std_min:              min grayscale std (drops flat/smooth background).
      - dark_fraction_min:         if not None, min fraction of dark (nuclei) pixels.
                                   Left as None (disabled) by default — the dark
                                   filter drops lymphocyte-sparse adipose/edges.

    Returns tiles_df with: tile_id, sample, x0, y0, x1, y1, cx_img, cy_img,
    tile_size, stride, and the per-tile tissue metrics.
    """
    H, W = img.shape[:2]
    records = []
    tile_idx = 0
    for y0 in range(0, H - tile_size + 1, stride):
        for x0 in range(0, W - tile_size + 1, stride):
            x1 = x0 + tile_size
            y1 = y0 + tile_size

            patch = img[y0:y1, x0:x1]
            m = compute_tissue_metrics(patch, white_thr=white_thr,)

            if m["valid_tissue_fraction"] < valid_tissue_fraction_min:
                continue
            if m["red_fraction"] > red_fraction_max:
                continue
            # key filters to remove pale smooth background
            if m["saturation_fraction"] < saturation_fraction_min:
                continue
            if m["gray_std"] < gray_std_min:
                continue
            # optional: useful for lymphocyte-rich thymus (disabled by default)
            if dark_fraction_min is not None and m["dark_fraction"] < dark_fraction_min:
                continue

            tile_id = f"{sample}_tile_{tile_idx:08d}" if sample is not None else f"tile_{tile_idx:08d}"
            records.append({
                "tile_id": tile_id,
                "sample": sample,
                "x0": int(x0),
                "y0": int(y0),
                "x1": int(x1),
                "y1": int(y1),
                "cx_img": float((x0 + x1) / 2),
                "cy_img": float((y0 + y1) / 2),
                "tile_size": int(tile_size),
                "stride": int(stride),
                "tissue_fraction": m["tissue_fraction"],
                "red_fraction": m["red_fraction"],
                "valid_tissue_fraction": m["valid_tissue_fraction"],
                "saturation_fraction": m["saturation_fraction"],
                "dark_fraction": m["dark_fraction"],
                "gray_std": m["gray_std"],
            })
                            
            tile_idx += 1

    tiles_df = pd.DataFrame(records)
    return tiles_df
    

def get_auto_crop_box_from_tiles(tiles_df, crop_size=3000, seed=0):
    if len(tiles_df) == 0:
        return None

    row = tiles_df.sample(1, random_state=seed).iloc[0]
    cx = float(row["cx_img"])
    cy = float(row["cy_img"])

    xmin = max(0, int(cx - crop_size / 2))
    xmax = int(cx + crop_size / 2)
    ymin = max(0, int(cy - crop_size / 2))
    ymax = int(cy + crop_size / 2)

    return (xmin, xmax, ymin, ymax)


def save_debug_overlay(img, centroids_df, gdf_boundaries, out_path, n_show=None, seed=0, outline_color="red", centroid_color="cyan",):
    if len(centroids_df) == 0:
        raise ValueError("No centroids to plot.")

    H, W = img.shape[:2]
    # --- sample centroids ---
    if n_show is None:
        show_df = centroids_df
    elif len(centroids_df) > n_show:
        show_df = centroids_df.sample(n_show, random_state=seed)
    else:
        show_df = centroids_df
    # --- figure ---
    fig = plt.figure(figsize=(W/1000, H/1000), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    # --- image ---
    ax.imshow(img, interpolation="nearest")
    # --- boundaries ---
    gdf_boundaries.boundary.plot(ax=ax, linewidth=0.3, alpha=0.6, color=outline_color,)
    # --- centroids ---
    ax.scatter(show_df["cx_img"], show_df["cy_img"], s=0.5, alpha=0.6, c=centroid_color, rasterized=False,)
    # --- axis ---
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    fig.savefig(out_path, dpi=600)
    plt.close(fig)
    print(f"✅ Saved overlay (centroids + boundaries): {out_path}")
    

def save_tile_grid_overlay(img, tiles_df, out_path=None, crop_box=None, sample_n=None, seed=0, line_width=0.4, alpha=0.9, tile_color="red", color_by=None, cmap_name="viridis", show_tile_id=False, max_labels=30, show=False,):
    """ Visualize tile boundaries on aligned H&E image. tiles_df expected columns: x0, y0, x1, y1, cx_img, cy_img, tile_id; optional: tissue_fraction """
    plot_df = tiles_df.copy()
    # keep tiles overlapping crop
    if crop_box is None:
        xmin, xmax, ymin, ymax = 0, img.shape[1], 0, img.shape[0]
    else:
        xmin, xmax, ymin, ymax = crop_box
        keep = (
            (plot_df["x1"] > xmin) &
            (plot_df["x0"] < xmax) &
            (plot_df["y1"] > ymin) &
            (plot_df["y0"] < ymax))
        plot_df = plot_df.loc[keep].copy()

    if sample_n is not None and len(plot_df) > sample_n:
        plot_df = plot_df.sample(sample_n, random_state=seed).copy()

    img_crop = img[int(ymin):int(ymax), int(xmin):int(xmax)]

    H, W = img_crop.shape[:2]
    fig_w = min(16, max(6, W / 1800))
    fig_h = min(16, max(6, H / 1800))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=200)
    ax.imshow(img_crop, interpolation="nearest")

    # prepare colors
    if color_by is not None:
        vals = plot_df[color_by].astype(float).values
        vmin = np.nanmin(vals)
        vmax = np.nanmax(vals)
        if vmax == vmin:
            vmax = vmin + 1e-8
        cmap = plt.get_cmap(cmap_name)
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
    else:
        cmap = None
        norm = None

    # draw rectangles
    for i, row in enumerate(plot_df.itertuples(index=False)):
        x0 = row.x0 - xmin
        y0 = row.y0 - ymin
        w = row.x1 - row.x0
        h = row.y1 - row.y0

        if color_by is not None:
            edgecolor = cmap(norm(getattr(row, color_by)))
        else:
            edgecolor = tile_color

        rect = Rectangle((x0, y0), w, h, fill=False, edgecolor=edgecolor, linewidth=line_width, alpha=alpha,)
        ax.add_patch(rect)

    # optionally label a few tile IDs
    if show_tile_id:
        n_label = min(max_labels, len(plot_df))
        label_df = plot_df.head(n_label)
        for _, row in label_df.iterrows():
            ax.text(row["x0"]-xmin+3, row["y0"]-ymin+10, str(row["tile_id"]), fontsize=5, color="yellow", bbox=dict(facecolor="black", alpha=0.4, pad=1),)

    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")

    if out_path is not None:
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"✅ Saved tile overlay: {out_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)


def save_random_tile_gallery(img, tiles_df, out_path=None, n_show=16, seed=0, title_col="tile_id", extra_col="tissue_fraction", show=False,):
    """ Show random tile patches in a gallery. """
    if len(tiles_df) == 0:
        raise ValueError("tiles_df is empty.")

    show_df = tiles_df.sample(min(n_show, len(tiles_df)), random_state=seed).copy()
    n = len(show_df)
    ncols = int(np.ceil(np.sqrt(n)))
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows), dpi=150)
    axes = np.array(axes).reshape(-1)

    for ax in axes:
        ax.axis("off")

    for ax, (_, row) in zip(axes, show_df.iterrows()):
        x0, y0, x1, y1 = int(row["x0"]), int(row["y0"]), int(row["x1"]), int(row["y1"])
        patch = img[y0:y1, x0:x1]
        ax.imshow(patch)

        title = str(row[title_col]) if title_col in row else ""
        if extra_col is not None and extra_col in row:
            title += f"\n{extra_col}={row[extra_col]:.2f}"
        ax.set_title(title, fontsize=8)
        ax.axis("off")

    plt.tight_layout()

    if out_path is not None:
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"✅ Saved tile gallery: {out_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)


def save_crop_debug_preview(img, crop_info, out_path, max_dim=2000):
    from PIL import Image, ImageDraw
    import numpy as np

    H, W = img.shape[:2]

    scale = min(max_dim / W, max_dim / H, 1.0)
    thumb_w = int(W * scale)
    thumb_h = int(H * scale)

    thumb = Image.fromarray(img).resize((thumb_w, thumb_h), Image.Resampling.BILINEAR)
    draw = ImageDraw.Draw(thumb)

    x0 = crop_info["crop_x0"] * scale
    y0 = crop_info["crop_y0"] * scale
    x1 = crop_info["crop_x1"] * scale
    y1 = crop_info["crop_y1"] * scale

    draw.rectangle([x0, y0, x1, y1], outline="yellow", width=5)
    thumb.save(out_path)

    print(f"Saved crop debug preview: {out_path}")


def crop_by_red_profile_thumbnail(img, downsample=32, red_r_min=220, red_g_max=80, red_b_max=80, red_frac_max=0.5, margin=200,):
    """ Automatically crop red padding using a downsampled thumbnail. Returns cropped full-res image and crop_info. """
    from PIL import Image
    import numpy as np

    H, W = img.shape[:2]

    # Downsample for robust and fast detection
    thumb_w = max(1, W // downsample)
    thumb_h = max(1, H // downsample)
    thumb = np.asarray(Image.fromarray(img).resize((thumb_w, thumb_h), Image.Resampling.BILINEAR))

    r, g, b = thumb[..., 0], thumb[..., 1], thumb[..., 2]
    red_mask = ((r >= red_r_min) & (g <= red_g_max) & (b <= red_b_max))

    col_red_frac = red_mask.mean(axis=0)
    row_red_frac = red_mask.mean(axis=1)

    good_cols = col_red_frac < red_frac_max
    good_rows = row_red_frac < red_frac_max

    def longest_true_run(mask):
        mask = np.asarray(mask, dtype=bool)
        if mask.sum() == 0:
            return None
        padded = np.r_[False, mask, False]
        changes = np.diff(padded.astype(np.int8))
        starts = np.where(changes == 1)[0]
        ends = np.where(changes == -1)[0]
        lengths = ends - starts
        i = int(np.argmax(lengths))
        return int(starts[i]), int(ends[i])

    x_run = longest_true_run(good_cols)
    y_run = longest_true_run(good_rows)

    if x_run is None or y_run is None:
        print("[WARN] Could not detect crop. Return original image.")
        crop_info = {"crop_x0": 0, "crop_y0": 0, "crop_x1": W, "crop_y1": H, "original_width": W, "original_height": H,}
        return img, crop_info

    x0_t, x1_t = x_run
    y0_t, y1_t = y_run

    # Map thumbnail coordinates back to full resolution
    x0 = max(0, x0_t * downsample - margin)
    x1 = min(W, x1_t * downsample + margin)
    y0 = max(0, y0_t * downsample - margin)
    y1 = min(H, y1_t * downsample + margin)

    img_crop = img[y0:y1, x0:x1].copy()

    crop_info = {
        "crop_x0": int(x0),
        "crop_y0": int(y0),
        "crop_x1": int(x1),
        "crop_y1": int(y1),
        "crop_width": int(x1 - x0),
        "crop_height": int(y1 - y0),
        "original_width": int(W),
        "original_height": int(H),
        "downsample": downsample,
        "red_r_min": red_r_min,
        "red_g_max": red_g_max,
        "red_b_max": red_b_max,
        "red_frac_max": red_frac_max,
        "margin": margin,
    }

    print("Cropping by red profile on thumbnail:")
    print(f"  original: H={H}, W={W}")
    print(f"  thumbnail: H={thumb_h}, W={thumb_w}")
    print(f"  crop box: x0={x0}, y0={y0}, x1={x1}, y1={y1}")
    print(f"  cropped: H={img_crop.shape[0]}, W={img_crop.shape[1]}")

    return img_crop, crop_info
    

def save_image_preview(img, out_path, max_dim=3000):
    """
    Save a downsampled preview of a large RGB image.
    """
    H, W = img.shape[:2]

    scale = min(max_dim / W, max_dim / H, 1.0)
    new_W = int(W * scale)
    new_H = int(H * scale)

    img_pil = Image.fromarray(img.astype("uint8"))
    img_preview = img_pil.resize((new_W, new_H), Image.Resampling.BILINEAR)

    img_preview.save(out_path)
    print(f"Saved preview: {out_path}")
    print(f"Original shape: {img.shape}")
    print(f"Preview size: {new_H} x {new_W}")

def crop_fullres_by_crop_info(img_raw, crop_info):
    x0 = int(crop_info["crop_x0"])
    y0 = int(crop_info["crop_y0"])
    x1 = int(crop_info["crop_x1"])
    y1 = int(crop_info["crop_y1"])

    img_crop = img_raw[y0:y1, x0:x1].copy()

    print("Full-res recrop:")
    print(f"  crop box: x0={x0}, y0={y0}, x1={x1}, y1={y1}")
    print(f"  cropped shape: {img_crop.shape}")

    return img_crop

def replace_red_background_with_white(img: np.ndarray, red_r_min: int = 220, red_g_max: int = 80, red_b_max: int = 80,):
    """ Replace red padding/background with white. """
    img_clean = img.copy()

    r = img_clean[..., 0]
    g = img_clean[..., 1]
    b = img_clean[..., 2]

    red_mask = ((r >= red_r_min) & (g <= red_g_max) & (b <= red_b_max))

    img_clean[red_mask] = 255

    print("Background replacement:")
    print(f"  red pixels replaced: {red_mask.mean():.4f}")

    return img_clean, {"red_fraction_full_image": float(red_mask.mean()),}


def replace_red_black_background_with_white(img: np.ndarray, red_r_min: int = 220, red_g_max: int = 80, red_b_max: int = 80, black_thr: int = 30,):
    """ Replace red padding/background and black/dark background with white. This keeps image shape, mpp, and coordinates unchanged. """
    img_clean = img.copy()

    r = img_clean[..., 0]
    g = img_clean[..., 1]
    b = img_clean[..., 2]

    # Pure/strong red padding
    red_mask = ((r >= red_r_min) & (g <= red_g_max) & (b <= red_b_max))

    # Very dark / black background
    # Use uint16 to avoid overflow when summing uint8 channels
    gray_sum = (r.astype(np.uint16) + g.astype(np.uint16) + b.astype(np.uint16))
    black_mask = gray_sum <= (3 * black_thr)
    bg_mask = red_mask | black_mask
    img_clean[bg_mask] = 255

    print("Background replacement:")
    print(f"  red pixels replaced:   {red_mask.mean():.4f}")
    print(f"  black pixels replaced: {black_mask.mean():.4f}")
    print(f"  total bg replaced:     {bg_mask.mean():.4f}")

    return img_clean, {
        "red_fraction_full_image": float(red_mask.mean()),
        "black_fraction_full_image": float(black_mask.mean()),
        "replaced_bg_fraction_full_image": float(bg_mask.mean()),}


def filter_tile_components(
    tiles_df: pd.DataFrame,
    image_shape: tuple,
    tile_size: int = 224,
    stride: int = 224,
    min_component_tiles: int = 80,
    min_component_height_tiles: int = 4,
    min_component_width_tiles: int = 4,
    thin_height_tiles: int = 3,
    thin_width_tiles: int = 30,
    remove_thin_horizontal: bool = True,
    remove_thin_vertical: bool = True,
):
    """
    Remove small isolated tile components and thin line-like components.

    This is useful for H&E-only tissue tile filtering:
    - removes small preview boxes
    - removes isolated debris/artifacts
    - removes top/bottom thin horizontal scanner lines
    - removes thin vertical artifacts if needed

    Parameters
    ----------
    tiles_df:
        Tile dataframe after initial tissue filtering.
        Must contain x0, y0.
    image_shape:
        Shape of the image used for tiling, e.g. img_mpp.shape.
    tile_size:
        Tile size in pixels.
    stride:
        Tile stride in pixels.
    min_component_tiles:
        Minimum number of connected tiles to keep.
    min_component_height_tiles:
        Minimum component height in tile rows.
    min_component_width_tiles:
        Minimum component width in tile columns.
    thin_height_tiles:
        Components with height <= this and very large width are removed.
    thin_width_tiles:
        Width threshold for thin horizontal components.

    Returns
    -------
    filtered_df:
        Tile dataframe after component filtering.
    component_info:
        DataFrame summarizing detected components.
    """
    if len(tiles_df) == 0:
        print("[WARN] tiles_df is empty before component filtering.")
        return tiles_df.copy(), pd.DataFrame()

    H, W = image_shape[:2]

    n_rows = (H - tile_size) // stride + 1
    n_cols = (W - tile_size) // stride + 1

    df = tiles_df.copy()
    df["grid_row"] = (df["y0"] // stride).astype(int)
    df["grid_col"] = (df["x0"] // stride).astype(int)

    valid = (
        (df["grid_row"] >= 0) &
        (df["grid_row"] < n_rows) &
        (df["grid_col"] >= 0) &
        (df["grid_col"] < n_cols)
    )
    df = df.loc[valid].copy()

    mask = np.zeros((n_rows, n_cols), dtype=np.uint8)
    mask[df["grid_row"].values, df["grid_col"].values] = 1

    n_labels, labels, stats, _ = cv2.connectedComponentsWithStats(
        mask,
        connectivity=8,
    )

    component_records = []
    keep_labels = []

    for label in range(1, n_labels):
        x = int(stats[label, cv2.CC_STAT_LEFT])
        y = int(stats[label, cv2.CC_STAT_TOP])
        w = int(stats[label, cv2.CC_STAT_WIDTH])
        h = int(stats[label, cv2.CC_STAT_HEIGHT])
        area = int(stats[label, cv2.CC_STAT_AREA])

        is_small = area < min_component_tiles
        is_too_short = h < min_component_height_tiles
        is_too_narrow = w < min_component_width_tiles

        is_thin_horizontal = (
            remove_thin_horizontal and
            (h <= thin_height_tiles) and
            (w >= thin_width_tiles)
        )

        is_thin_vertical = (
            remove_thin_vertical and
            (w <= thin_height_tiles) and
            (h >= thin_width_tiles)
        )

        keep = not (
            is_small or
            is_too_short or
            is_too_narrow or
            is_thin_horizontal or
            is_thin_vertical
        )

        component_records.append({
            "component_id": label,
            "grid_x0": x,
            "grid_y0": y,
            "grid_width": w,
            "grid_height": h,
            "n_tiles": area,
            "is_small": is_small,
            "is_too_short": is_too_short,
            "is_too_narrow": is_too_narrow,
            "is_thin_horizontal": is_thin_horizontal,
            "is_thin_vertical": is_thin_vertical,
            "keep": keep,
        })

        if keep:
            keep_labels.append(label)

    component_info = pd.DataFrame(component_records)

    df["component_id"] = labels[df["grid_row"].values, df["grid_col"].values]

    if len(component_info) > 0:
        size_map = component_info.set_index("component_id")["n_tiles"].to_dict()
        height_map = component_info.set_index("component_id")["grid_height"].to_dict()
        width_map = component_info.set_index("component_id")["grid_width"].to_dict()

        df["component_size"] = df["component_id"].map(size_map)
        df["component_height_tiles"] = df["component_id"].map(height_map)
        df["component_width_tiles"] = df["component_id"].map(width_map)

    filtered_df = df[df["component_id"].isin(keep_labels)].copy()

    print("Spatial component filtering:")
    print(f"  components found: {n_labels - 1}")
    print(f"  components kept:  {len(keep_labels)}")
    print(f"  tiles before:      {len(tiles_df)}")
    print(f"  tiles after:       {len(filtered_df)}")

    if len(component_info) > 0:
        print(component_info.sort_values("n_tiles", ascending=False).head(10))

    return filtered_df, component_info
    


import re
from pathlib import Path
import pandas as pd


def extract_valid_data_envelopes_from_fic(fic_path: str | Path):
    """
    Extract valid data envelopes from Philips .fic metadata.

    Each envelope should contain 4 vertices.
    """
    fic_path = Path(fic_path)
    text = fic_path.read_text(errors="ignore")

    blocks = re.findall(
        r'<DataObject ObjectType="UFSImageValidDataEnvelope">'
        r'(.*?<Attribute Name="UFS_IMAGE_VALID_DATA_ENVELOPE_DIMENSIONS".*?</Attribute>\s*</DataObject>)',
        text,
        flags=re.S,
    )

    bboxes = []

    for i, block in enumerate(blocks, start=1):
        pts = re.findall(
            r'<Attribute Name="UFS_IMAGE_OPP_EXTREME_VERTEX"[^>]*>\s*(\d+)\s+(\d+)\s*</Attribute>',
            block,
            flags=re.S,
        )

        pts = [(int(x), int(y)) for x, y in pts]

        if len(pts) == 0:
            continue

        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]

        x0, x1 = min(xs), max(xs)
        y0, y1 = min(ys), max(ys)

        bboxes.append({
            "envelope_id": i,
            "x0": x0,
            "y0": y0,
            "x1": x1,
            "y1": y1,
            "width": x1 - x0,
            "height": y1 - y0,
            "area": (x1 - x0) * (y1 - y0),
            "n_vertices": len(pts),
            "vertices": pts,
        })

    return bboxes


from PIL import Image, ImageDraw, ImageFont
from pathlib import Path
import numpy as np


def save_fic_bboxes_overlay(
    img: np.ndarray,
    bboxes: list[dict],
    out_path: str | Path,
    max_dim: int = 3000,
    line_width: int = 6,
):
    """
    Save a downsampled TIFF preview with .fic valid data envelopes overlaid.

    img:
        RGB image, H x W x 3.
    bboxes:
        Output from extract_valid_data_envelopes_from_fic().
        Each bbox should contain x0, y0, x1, y1, envelope_id.
    """
    out_path = Path(out_path)

    H, W = img.shape[:2]

    # Downsample preview
    scale = min(max_dim / W, max_dim / H, 1.0)
    preview_w = int(W * scale)
    preview_h = int(H * scale)

    preview = Image.fromarray(img.astype(np.uint8)).resize(
        (preview_w, preview_h),
        Image.Resampling.BILINEAR,
    )

    draw = ImageDraw.Draw(preview)

    # Use different colors for different envelopes
    colors = ["yellow", "lime", "cyan", "magenta", "orange", "blue"]

    for i, bbox in enumerate(bboxes):
        color = colors[i % len(colors)]

        x0 = bbox["x0"] * scale
        y0 = bbox["y0"] * scale
        x1 = bbox["x1"] * scale
        y1 = bbox["y1"] * scale

        draw.rectangle(
            [x0, y0, x1, y1],
            outline=color,
            width=line_width,
        )

        label = (
            f"env {bbox.get('envelope_id', i + 1)} "
            f"W={bbox['width']} H={bbox['height']}"
        )

        draw.text(
            (x0 + 10, y0 + 10),
            label,
            fill=color,
        )

    preview.save(out_path)

    print(f"Saved bbox overlay: {out_path}")
    print(f"Original image shape: H={H}, W={W}")
    print(f"Preview size: H={preview_h}, W={preview_w}")



def read_slide_level(tiff_path: str | Path, level: int) -> tuple[np.ndarray, float, float]:
    """
    Read a single pyramid level of a WSI as RGB uint8, plus the exact x/y scale
    factors mapping full-resolution (level-0) pixels into that level.

    Reads ONLY the requested level, so it stays memory-friendly for whole-slide
    images. Format-aware: NDPI (and other OpenSlide-only formats) are read via
    OpenSlide; OME-TIFF and pyramidal TIFF via tifffile's series.levels. Scale
    factors are derived from the actual level shape rather than assuming 2**level,
    since pyramid levels are not always exact powers of two.
    """
    path = Path(tiff_path)
    name = path.name.lower()

    if name.endswith((".ndpi", ".svs", ".scn", ".mrxs", ".vms", ".vmu", ".bif")):
        slide = openslide.OpenSlide(str(path))
        try:
            full_w, full_h = slide.level_dimensions[0]
            lvl_w, lvl_h = slide.level_dimensions[level]
            img = np.asarray(slide.read_region((0, 0), level, (lvl_w, lvl_h)).convert("RGB"))
        finally:
            slide.close()
    else:  # ome-tiff / pyramidal tiff
        with tifffile.TiffFile(str(path)) as tf:
            series = tf.series[0]
            full_h, full_w = series.shape[:2]
            img = to_rgb_uint8(series.levels[level].asarray())
        lvl_h, lvl_w = img.shape[:2]

    sx = lvl_w / full_w
    sy = lvl_h / full_h
    print(f"level {level} shape: {img.shape}, scale x,y: {sx:.6f}, {sy:.6f}")
    return img, sx, sy


def overlay_clusters_on_slide(
    he_adata,
    sample: str,
    tiff_path: str | Path,
    orig_mpp: float,
    clusters=None,
    label_col: str = "leiden_m2_r0.5",
    level: int = 5,
    target_mpp: float = 0.25,
    alpha: float = 0.4,
    cmap_name: str = "tab20",
    sample_col: str = "SampleID",
    coord_cols: tuple[str, str, str, str] = ("x0", "x1", "y0", "y1"),
    palette: dict | str | None = None,
    uns_key: str | None = None,
    out_path: str | Path | None = None,
) -> tuple[np.ndarray, dict]:
    """
    Paint morpho-cluster labels for one sample onto a downsampled level of its WSI.

    Patch boxes are given in the *resampled* mpp-image space (the space tiling was
    done in, at `target_mpp`). They are mapped back to original level-0 WSI pixels
    (via inv_scale = target_mpp / orig_mpp, matching tiles_to_original_coords /
    resample_to_target_mpp) and then scaled into the chosen pyramid `level`.

    Parameters
    ----------
    he_adata:
        AnnData whose .obs holds one row per H&E patch, with `sample_col`,
        `label_col`, and the patch bounds `coord_cols` (in resampled mpp space).
    sample:
        SampleID to plot (subsets `he_adata.obs[sample_col]`).
    tiff_path:
        Path to this sample's WSI. May be OME-TIFF/pyramidal TIFF or NDPI; the
        reader dispatches on the suffix (see read_slide_level).
    orig_mpp:
        Native resolution (microns/pixel) of this sample's slide at level 0. With
        `target_mpp`, sets the inverse resample factor for the coordinate mapping.
    clusters:
        Cluster label(s) to highlight (a single value or a list of values of
        `label_col`). Patches whose label is not in this list are left unpainted
        so the raw H&E shows through. If None, all clusters are painted.
    label_col:
        Categorical column in .obs holding the cluster label per patch.
    level:
        Pyramid level to render (higher = coarser / smaller, faster overview).
    target_mpp:
        Resolution the tiles were generated at (default 0.25, from preprocessing).
    alpha:
        Blend strength of the colour overlay over the H&E (0..1).
    cmap_name / palette / uns_key:
        Colours are assigned from the FULL set of `label_col` categories so each
        cluster keeps a stable colour regardless of which subset is requested.
        `palette` selects the source:
          - None (default): generate from `cmap_name`.
          - "uns": use the scanpy colour list `he_adata.uns[uns_key]` (defaults to
            `he_adata.uns[f"{label_col}_colors"]`), zipped to `.cat.categories`.
          - dict (label -> RGB 0..255 or any matplotlib colour): used as given.
    coord_cols:
        (x0, x1, y0, y1) column names in .obs, in resampled mpp-image space.
    out_path:
        If given, the blended image is saved here (PNG).

    Returns
    -------
    (out_img, palette):
        The blended RGB uint8 image and the label -> RGB colour map used.
    """
    # --- stable colour per cluster, from the full category set ---
    from matplotlib.colors import to_rgb

    col = he_adata.obs[label_col]
    if isinstance(col.dtype, pd.CategoricalDtype):
        cats = list(col.cat.categories)
    else:
        cats = sorted(pd.unique(col.dropna()))

    if palette == "uns":
        key = uns_key or f"{label_col}_colors"
        colors = list(he_adata.uns[key])
        if len(colors) < len(cats):
            raise ValueError(f"he_adata.uns['{key}'] has {len(colors)} colours for {len(cats)} categories")
        palette = {c: tuple(int(v * 255) for v in to_rgb(colors[i])) for i, c in enumerate(cats)}
    elif palette is None:
        cmap = plt.get_cmap(cmap_name, len(cats))
        palette = {c: tuple(int(v * 255) for v in cmap(i)[:3]) for i, c in enumerate(cats)}
    else:
        # dict: normalise any matplotlib colour spec (hex/name/0..1 or 0..255 tuple) to RGB 0..255
        def _to_rgb255(c):
            if isinstance(c, str):
                return tuple(int(v * 255) for v in to_rgb(c))
            if max(c) > 1:  # already 0..255
                return tuple(int(v) for v in c[:3])
            return tuple(int(v * 255) for v in c[:3])
        palette = {k: _to_rgb255(v) for k, v in palette.items()}

    # --- which clusters to paint ---
    if clusters is None:
        show = set(map(str, cats))
    elif np.isscalar(clusters):
        show = {str(clusters)}
    else:
        show = set(map(str, clusters))
    missing = show - set(map(str, cats))
    if missing:
        print(f"[WARN] clusters not found in '{label_col}': {sorted(missing)}")

    # --- read the requested slide level ---
    img, sx, sy = read_slide_level(tiff_path, level)
    lvl_h, lvl_w = img.shape[:2]

    obs = he_adata.obs[he_adata.obs[sample_col] == sample]
    if len(obs) == 0:
        raise ValueError(f"No patches for sample '{sample}' in he_adata.obs['{sample_col}']")

    # resampled mpp-space -> original level-0 pixels (inv_scale), then -> this level.
    # inv_scale matches tiles_to_original_coords / resample_to_target_mpp.
    inv_scale = 1.0 if abs(orig_mpp - target_mpp) <= 1e-4 else target_mpp / orig_mpp
    x0c, x1c, y0c, y1c = coord_cols
    x0 = np.round(obs[x0c].to_numpy() * inv_scale * sx).astype(int).clip(0, lvl_w)
    x1 = np.round(obs[x1c].to_numpy() * inv_scale * sx).astype(int).clip(0, lvl_w)
    y0 = np.round(obs[y0c].to_numpy() * inv_scale * sy).astype(int).clip(0, lvl_h)
    y1 = np.round(obs[y1c].to_numpy() * inv_scale * sy).astype(int).clip(0, lvl_h)
    labels = obs[label_col].astype("object").to_numpy()

    # --- paint only the requested clusters into an overlay raster ---
    overlay = np.zeros((lvl_h, lvl_w, 3), np.float32)
    mask = np.zeros((lvl_h, lvl_w), bool)
    n_painted = 0
    for i in range(len(obs)):
        lab = labels[i]
        if lab is None or (isinstance(lab, float) and np.isnan(lab)) or str(lab) not in show:
            continue
        overlay[y0[i]:y1[i], x0[i]:x1[i]] = palette[lab]
        mask[y0[i]:y1[i], x0[i]:x1[i]] = True
        n_painted += 1

    # --- single alpha blend where painted patches exist ---
    out = img.copy()
    out[mask] = ((1 - alpha) * img[mask] + alpha * overlay[mask]).astype(np.uint8)

    print(f"Sample {sample}: painted {n_painted}/{len(obs)} patches for clusters {sorted(show)}")

    if out_path is not None:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        Image.fromarray(out).save(out_path)
        print(f"✅ Saved cluster overlay: {out_path}")

    return out, palette


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Preprocess H&E WSIs into tiles.")
    parser.add_argument(
        "--index",
        type=int,
        default=None,
        help="1-based index of a single sample to process (e.g. LSF $LSB_JOBINDEX). "
             "If omitted, all samples are processed sequentially.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("/lustre/scratch126/cellgen/teichmann/lm25/data/thymus_he/thyAgeing_HnE"),
        help="Directory containing the raw whole-slide images (and where outputs are written).",
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="*.ndpi",
        help="Glob pattern (relative to --data-dir) selecting input files, e.g. "
             "'*.ndpi' or '*.ome.tif'. The per-file format is inferred from its "
             "suffix. Ignored when --file-list is given.",
    )
    parser.add_argument(
        "--file-list",
        type=Path,
        default=None,
        help="Text file with one input filename/path per line (may mix formats). "
             "Relative entries are resolved against --data-dir. Overrides --pattern.",
    )
    # Tile-keep thresholds. Defaults tuned via diagnose_tile_filtering.py to retain
    # tissue edges / interstitial space / adipose while rejecting empty background.
    parser.add_argument("--valid-tissue-fraction-min", type=float, default=0.15,
                        help="Min non-white, non-red fraction to keep a tile.")
    parser.add_argument("--red-fraction-max", type=float, default=0.05,
                        help="Drop tiles whose red-artifact fraction exceeds this.")
    parser.add_argument("--saturation-fraction-min", type=float, default=0.01,
                        help="Min fraction of colour-saturated (stained) pixels.")
    parser.add_argument("--gray-std-min", type=float, default=4.5,
                        help="Min grayscale std (drops flat/smooth background).")
    parser.add_argument("--dark-fraction-min", type=float, default=None,
                        help="If set, min fraction of dark (nuclei) pixels. Disabled "
                             "by default because it drops adipose/edge tissue.")
    args = parser.parse_args()

    # Directory containing the raw whole-slide images.
    data_dir = args.data_dir

    # Build the sample list either from an explicit file list or a glob pattern.
    # The input format of each sample is determined per-file from its suffix, so
    # .ndpi (OpenSlide) and .ome.tif (tifffile) can be processed by the same run.
    if args.file_list is not None:
        samples = get_samples_from_list(args.file_list, base_dir=data_dir)
        print(f"Found {len(samples)} files in list {args.file_list}")
    else:
        samples = get_tiff_samples(data_dir, pattern=args.pattern)
        print(f"Found {len(samples)} files matching '{args.pattern}' in {data_dir}")

    # In array-job mode, select a single sample by 1-based index. The sorted glob /
    # ordered file list gives a stable ordering, so index N always maps to the
    # same file across array tasks.
    n_total = len(samples)
    if args.index is not None:
        if not (1 <= args.index <= n_total):
            raise IndexError(
                f"--index {args.index} out of range 1..{n_total}"
            )
        samples = [samples[args.index - 1]]
        print(f"Array mode: processing sample {args.index}/{n_total}: "
              f"{samples[0]['sample_name']}")

    for item in samples:
        sample_name,he_path = item["sample_name"],item["file_path"]

        print("\n" + "=" * 80)
        print(f"Processing sample: {sample_name}")
        print(f"H&E path: {he_path}\n")

        mpp = 0.25                 # target resolution to resample to
        fallback_orig_mpp = 0.227  # used only if the slide has no mpp metadata
        tile_size = 224
        stride = 224
        mpp_tag = "mpp025"
        part_name = f"{sample_name}_part0"  # single-part slide; matches test_morpho.py's {sample}_part* glob

        # Idempotent skip: if this sample's outputs already exist, don't redo the
        # expensive WSI read/tiling. These paths must match the ones written at the
        # end of the loop.
        out_dir = data_dir / f"{sample_name}_processed"
        image_path = out_dir / f"{part_name}_image_he_{mpp_tag}.npy"
        tiles_path = out_dir / f"{part_name}_he_tiles_{mpp_tag}_tile{tile_size}_stride{stride}.parquet"
        coord_path = out_dir / f"{part_name}_tile_coordinates_he_{mpp_tag}.pkl"
        if image_path.exists() and tiles_path.exists() and coord_path.exists():
            print(f"[SKIP] Outputs already exist for {sample_name}; skipping preprocessing.")
            continue

        # Check native resolution first (metadata only, no pixel read) so we can
        # skip coarse slides before paying the multi-GB read / creating an output dir.
        # read_native_mpp dispatches on the file suffix (OpenSlide vs OME-TIFF).
        native_mpp = read_native_mpp(he_path)
        orig_mpp = native_mpp if native_mpp is not None else fallback_orig_mpp
        if native_mpp is None:
            print(f"[WARN] No mpp in slide metadata; falling back to orig_mpp={orig_mpp}")

        # Skip slides coarser than the target: native mpp larger than target would
        # require upscaling (inventing detail), which we don't want.
        if orig_mpp > mpp + 1e-4:
            print(f"[SKIP] {sample_name}: native mpp {orig_mpp} is coarser than target {mpp}; "
                  f"would require upscaling. Skipping.")
            continue

        out_dir = data_dir / f"{sample_name}_processed"
        out_dir.mkdir(parents=True, exist_ok=True)

        print(f"Using orig_mpp={orig_mpp}, target_mpp={mpp}")
        print(f"Reading whole-slide image ({detect_format(he_path)})...")
        img_raw, _ = read_image_as_rgb(he_path, level=0)
        print("Raw image shape:", img_raw.shape)


        # Optional Philips .fic valid-data-envelope QC overlay. Hamamatsu NDPI
        # slides won't have a .fic sidecar, so this is skipped when absent.
        fic_path = data_dir / f"{sample_name}-1_1_1.fic"
        if fic_path.exists():
            bboxes = extract_valid_data_envelopes_from_fic(fic_path)
            for b in bboxes:
                print(b)
            df = pd.DataFrame(bboxes)
            print(df[["envelope_id", "x0", "y0", "x1", "y1", "width", "height", "area", "n_vertices"]])
            save_fic_bboxes_overlay(img=img_raw, bboxes=bboxes, out_path=out_dir / "fic_valid_envelopes_overlay.png", max_dim=3000,)
        else:
            print(f"No .fic file at {fic_path}; skipping envelope overlay.")

        # Resample from native orig_mpp to the target mpp BEFORE tiling so that
        # tile coordinates and tile_size * mpp physical sizes are correct.
        # (Kept after the .fic overlay above, whose bbox coords are in native pixels.)
        img_raw = resample_to_target_mpp(img_raw, orig_mpp=orig_mpp, target_mpp=mpp)
        print("Resampled image shape:", img_raw.shape)

        print("Replacing red/black background with white...")
        # img_mpp, bg_info = replace_red_background_with_white(img_raw, red_r_min=220, red_g_max=80, red_b_max=80,)
        img_mpp, bg_info = replace_red_black_background_with_white(img_raw, red_r_min=220, red_g_max=80, red_b_max=80, black_thr=30,)
        save_image_preview(img_mpp, out_dir / "img_white_background_preview.png", max_dim=3000,)
        
        # sys.exit()
        # print("Reading H&E TIFF...")
        # img_raw = read_he_tiff_as_rgb(he_path)
        # print("Raw image shape:", img_raw.shape)
        
        # print("Cropping red padding/background...")
        # img_mpp, crop_info = crop_by_red_profile_thumbnail(img_raw, downsample=32, red_r_min=220, red_g_max=80, red_b_max=80, red_frac_max=0.5, margin=200,)
        # print("Cropped image shape:", img_mpp.shape)
        # print(f"Tile physical size: {tile_size * mpp} µm")
        
        # save_crop_debug_preview(img_raw, crop_info, out_dir / "crop_debug_preview.png",)
        # save_image_preview(img_mpp, out_dir / "img_mpp_preview.png", max_dim=3000,)
                
        # # sys.exit()
        # # img_mpp = crop_fullres_by_crop_info(img_raw, crop_info)
        
        # print("Reading H&E TIFF...")
        # img_mpp = read_he_tiff_as_rgb(he_path)
        
        print("Generating H&E grid tiles...")
        print(f"Tile-keep thresholds: valid_tissue_fraction_min={args.valid_tissue_fraction_min}, "
              f"red_fraction_max={args.red_fraction_max}, "
              f"saturation_fraction_min={args.saturation_fraction_min}, "
              f"gray_std_min={args.gray_std_min}, dark_fraction_min={args.dark_fraction_min}")
        tiles_df = generate_he_tiles(
            img=img_mpp,
            tile_size=tile_size,
            stride=stride,
            valid_tissue_fraction_min=args.valid_tissue_fraction_min,
            red_fraction_max=args.red_fraction_max,
            saturation_fraction_min=args.saturation_fraction_min,
            gray_std_min=args.gray_std_min,
            dark_fraction_min=args.dark_fraction_min,
            white_thr=230,
            sample=sample_name,
        )
        print("Number of valid tiles:", len(tiles_df))

        print("\nNumber of valid tiles before component filtering:", len(tiles_df))
        tiles_df, component_info = filter_tile_components(
            tiles_df=tiles_df,
            image_shape=img_mpp.shape,
            tile_size=tile_size,
            stride=stride,
            min_component_tiles=80,
            min_component_height_tiles=4,
            min_component_width_tiles=4,
            thin_height_tiles=3,
            thin_width_tiles=30,
            remove_thin_horizontal=True,
            remove_thin_vertical=True,
        )
        print("Number of valid tiles after component filtering:", len(tiles_df))
                
        tiles_df["mpp"] = mpp
        tiles_df["tile_physical_size_um"] = tile_size * mpp
        tiles_df["stride_physical_um"] = stride * mpp
        tiles_df["part_name"] = part_name
        tiles_df["source_image"] = f"{part_name}_image_he_{mpp_tag}.npy"
        tiles_df["original_file"] = str(he_path)
        tiles_df["data_type"] = "HE_only"
        tiles_df["has_xenium"] = False
        
        if len(tiles_df) > 0:
            print(tiles_df.head())
            print(tiles_df["tissue_fraction"].describe())
        
        save_tile_grid_overlay(img=img_mpp, tiles_df=tiles_df, out_path=out_dir / f"tiles_overlay_full_{mpp_tag}_.png", crop_box=None, sample_n=None, seed=0, line_width=0.3, alpha=0.8, tile_color="red", color_by=None, show=False,)
        crop_box = get_auto_crop_box_from_tiles(tiles_df, crop_size=3000, seed=0,)
        save_tile_grid_overlay(img=img_mpp, tiles_df=tiles_df, out_path=out_dir / f"tiles_overlay_crop_{mpp_tag}_.png", crop_box=crop_box, sample_n=None, line_width=0.8, alpha=1.0, tile_color="lime", show=False,)        
        save_random_tile_gallery(img=img_mpp, tiles_df=tiles_df, out_path=out_dir / f"tiles_random_gallery_{mpp_tag}_.png", n_show=25, seed=0, show=False,)
        
        image_path = out_dir / f"{part_name}_image_he_{mpp_tag}.npy"
        np.save(image_path, img_mpp)

        tiles_path = out_dir / f"{part_name}_he_tiles_{mpp_tag}_tile{tile_size}_stride{stride}.parquet"
        tiles_df.to_parquet(tiles_path, index=False)
        tile_coordinates = {row["tile_id"]: (row["cx_img"], row["cy_img"]) for _, row in tiles_df.iterrows()}

        coord_path = out_dir / f"{part_name}_tile_coordinates_he_{mpp_tag}.pkl"
        with open(coord_path, "wb") as f:
            pickle.dump(tile_coordinates, f)

        print(f"Saved image: {image_path}")
        print(f"Saved tiles: {tiles_path}")
        print(f"Saved tile coordinates: {coord_path}")

    print("All samples done.")


if __name__ == "__main__":
    main()
    