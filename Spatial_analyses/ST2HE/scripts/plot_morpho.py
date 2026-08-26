import tifffile
import numpy as np
import json
import numpy as np
import matplotlib.cm as cm

def tiles_to_original_coords(tiles_df, orig_mpp, target_mpp=0.25):
    """
    Map tile boxes (in the resampled mpp-image space) back to original
    level-0 WSI pixel coordinates.

    resampled = original * scale, with scale = orig_mpp / target_mpp
      => original = resampled / scale = resampled * (target_mpp / orig_mpp)
    """
    # Matches resample_to_target_mpp: if orig ~= target, no resampling was done.
    if abs(orig_mpp - target_mpp) <= 1e-4:
        inv_scale = 1.0
    else:
        scale = orig_mpp / target_mpp
        inv_scale = 1.0 / scale          # = target_mpp / orig_mpp

    out = tiles_df.copy()
    for col in ["x0", "x1", "y0", "y1", "cx_img", "cy_img"]:
        if col in out.columns:
            out[col + "_orig"] = out[col] * inv_scale

    return out, orig_mpp, inv_scale

from matplotlib.colors import to_rgb

def color_to_rgb255(c):
    """Any matplotlib color (hex, name, rgb tuple) -> [r, g, b] ints 0-255 for QuPath."""
    r, g, b = to_rgb(c)                       # floats 0..1
    return [int(round(r * 255)), int(round(g * 255)), int(round(b * 255))]


def export_niches_to_geojson(he_adata, sample, niche_col, orig_mpp, palette, out_dir = f"./morpho_niches"):
    """
    Export the patches of a given sample to a GeoJSON file with their niche classification.
    
    Parameters:
    - he_adata: AnnData object containing the data.
    - sample: Sample ID to filter the data.
    - niche_col: Column name in he_adata.obs that contains the niche classification.
    - orig_mpp: Original microns per pixel value.
    - palette: Dictionary mapping niche labels to colors.
    """

    obs = he_adata.obs[he_adata.obs['SampleID'] == sample].copy()
    obs = tiles_to_original_coords(obs, orig_mpp)[0]
    
    x0 = np.round(obs["x0_orig"].to_numpy()).astype(int)
    x1 = np.round(obs["x1_orig"].to_numpy()).astype(int)
    y0 = np.round(obs["y0_orig"].to_numpy()).astype(int)
    y1 = np.round(obs["y1_orig"].to_numpy()).astype(int)

    label_col = niche_col
    cats = list(obs[label_col].cat.categories)
    codes = obs[label_col].cat.codes.to_numpy()       # int per patch, -1 if NaN
    labels = obs[label_col].astype(str).to_numpy()

    features = []
    for i in range(len(obs)):
        a, b, c, d = float(x0[i]), float(x1[i]), float(y0[i]), float(y1[i])
        ring = [[a, c], [b, c], [b, d], [a, d], [a, c]]   # closed, CW
        features.append({
            "type": "Feature",
            "geometry": {"type": "Polygon", "coordinates": [ring]},
            "properties": {
                "objectType": "annotation",
                "classification": {"name": labels[i], "color": color_to_rgb255(palette[labels[i]])},
            },
        })

    fc = {"type": "FeatureCollection", "features": features}
    with open(f"{out_dir}/{sample}_{niche_col}_qupath.geojson", "w") as f:
        json.dump(fc, f)
    print(len(features), "patches written")