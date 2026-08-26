"""
Self-contained Macenko H&E stain normalization (no external stain-norm deps).

Rationale
---------
UNI2-h embeds stain/scanner appearance very strongly. In this project the Xenium
H&E (post-Xenium .ome.tif) and the standalone H&E (Hamamatsu .ndpi) come from
different acquisition paths, and PC1 of the L2-normed embeddings almost perfectly
separates the two platforms (t ~ 1e3). Mapping every tile onto a single canonical
H&E stain target before the encoder removes the dominant colour axis so Harmony
only has to close the residual gap.

Every source tile is normalized to a *fixed* canonical target (not to a reference
tile), so both platforms are mapped to the same appearance regardless of which
slide a tile came from. Degenerate tiles (mostly white / padded borders / too
little tissue) are returned unchanged rather than crashing or producing garbage.

Reference: Macenko et al., "A method for normalizing histology slides for
quantitative analysis", ISBI 2009. Target matrix/concentrations are the widely
used staintools defaults.
"""

from __future__ import annotations

import numpy as np
from PIL import Image


class MacenkoNormalizer:
    """
    Macenko normalizer to a fixed canonical H&E target.

    Parameters
    ----------
    Io : int
        Transmitted light intensity (background). 240 is standard for H&E.
    alpha : float
        Percentile for robust min/max of the stain angle (1 -> use 1st/99th).
    beta : float
        Optical-density threshold below which a pixel is treated as background
        and excluded from stain-vector estimation.
    min_tissue_px : int
        If fewer than this many pixels survive the OD/beta mask, the tile is
        considered background-dominated and returned unchanged.
    """

    # Canonical target stain matrix (columns = haematoxylin, eosin) and the
    # 99th-percentile target concentrations. staintools/Macenko defaults.
    HERef = np.array([[0.5626, 0.2159],
                      [0.7201, 0.8012],
                      [0.4062, 0.5581]], dtype=np.float64)
    maxCRef = np.array([1.9705, 1.0308], dtype=np.float64)

    def __init__(self, Io: int = 240, alpha: float = 1.0, beta: float = 0.15,
                 min_tissue_px: int = 100):
        self.Io = float(Io)
        self.alpha = float(alpha)
        self.beta = float(beta)
        self.min_tissue_px = int(min_tissue_px)

    def normalize(self, I: np.ndarray) -> np.ndarray:
        """
        Normalize an RGB uint8 image (H, W, 3) to the canonical H&E target.

        Returns a uint8 (H, W, 3) array. On any degeneracy (too little tissue,
        singular decomposition) the input is returned unchanged so the pipeline
        never breaks on padded/white edge tiles.
        """
        if I.ndim != 3 or I.shape[-1] != 3:
            raise ValueError(f"Expected RGB (H,W,3), got {I.shape}")

        h, w, _ = I.shape
        try:
            return self._normalize_impl(I, h, w)
        except (np.linalg.LinAlgError, ValueError, FloatingPointError):
            # Degenerate tile -> leave it untouched.
            return I if I.dtype == np.uint8 else np.clip(I, 0, 255).astype(np.uint8)

    def _normalize_impl(self, I: np.ndarray, h: int, w: int) -> np.ndarray:
        Io, alpha, beta = self.Io, self.alpha, self.beta

        Iflat = I.reshape((-1, 3)).astype(np.float64)

        # Optical density.
        OD = -np.log((Iflat + 1.0) / Io)

        # Keep only pixels with meaningful stain (drop near-transparent/background).
        mask = ~np.any(OD < beta, axis=1)
        ODhat = OD[mask]
        if ODhat.shape[0] < self.min_tissue_px:
            return I if I.dtype == np.uint8 else np.clip(I, 0, 255).astype(np.uint8)

        # Plane of the two dominant stain directions (eigh -> ascending eigenvalues).
        _, eigvecs = np.linalg.eigh(np.cov(ODhat.T))
        proj = ODhat.dot(eigvecs[:, 1:3])

        # Robust angular extremes -> the two stain vectors.
        phi = np.arctan2(proj[:, 1], proj[:, 0])
        minPhi = np.percentile(phi, alpha)
        maxPhi = np.percentile(phi, 100.0 - alpha)

        vMin = eigvecs[:, 1:3].dot(np.array([np.cos(minPhi), np.sin(minPhi)]))
        vMax = eigvecs[:, 1:3].dot(np.array([np.cos(maxPhi), np.sin(maxPhi)]))

        # Order columns as (haematoxylin, eosin): H has the larger red-channel OD.
        if vMin[0] > vMax[0]:
            HE = np.column_stack((vMin, vMax))
        else:
            HE = np.column_stack((vMax, vMin))

        # Per-pixel stain concentrations for the whole tile (2 x N).
        C = np.linalg.lstsq(HE, OD.T, rcond=None)[0]

        # Robust per-stain max, then rescale to the canonical target.
        maxC = np.array([np.percentile(C[0, :], 99),
                         np.percentile(C[1, :], 99)])
        maxC = np.where(maxC < 1e-6, 1.0, maxC)  # guard empty/degenerate stain
        C = C / maxC[:, None] * self.maxCRef[:, None]

        # Reconstruct in the canonical target stain space.
        Inorm = Io * np.exp(-self.HERef.dot(C))
        Inorm = np.clip(Inorm, 0, 255).T.reshape(h, w, 3).astype(np.uint8)
        return Inorm

    def normalize_pil(self, img: "Image.Image") -> "Image.Image":
        """Normalize a PIL RGB image, returning a PIL RGB image."""
        arr = np.asarray(img.convert("RGB"))
        return Image.fromarray(self.normalize(arr))
