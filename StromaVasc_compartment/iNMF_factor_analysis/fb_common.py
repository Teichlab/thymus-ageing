"""Shared constants and plotting helpers for the final fibroblast-factor notebooks (01-05).

The analysis itself lives in the notebooks. This module only holds what has to be identical
across them: paths, the factor vocabulary, the gene-signature loader, the colour scheme and a
handful of figure helpers (Nature-width figures, heatmaps, UMAP point clouds, the usage
heatmap and the niche dotplot). Nothing here reads a large object.
"""
import csv
import os
import warnings
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, to_hex, to_rgb
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# 1. paths
# ---------------------------------------------------------------------------
REPO = '/nfs/team205/vk8/projects/thymus_ageing_atlas'
FB = f'{REPO}/FB_Vascular_compartment'
FINAL = f'{FB}/notebooks/revisions/final'

# the joint scRNA-seq + multiome fibroblast object (HDF5 despite the .zarr extension)
OBJ = f'{FB}/data/objects/thyAgeing_Fb_scell_multiome_nuclei_anno_2026-08-20.zarr'
# the current harmonised annotation of that object (one label vocabulary in both arms)
HARMONIZED_OBS = f'{FB}/data/objects/thyAgeing_Fb_scell_multiome_obs_harmonized_2026-09-06.csv'
HARMONIZED_LABEL = 'taa_l3_harmonized'
# the editable gene signatures (long format: role, panel, gene)
SIGNATURES_CSV = f'{FB}/data/gene_sets/fb_factor_signatures.csv'

NMF_DIR = f'{FB}/data/nmf_joint'
INMF_DIR = f'{NMF_DIR}/inmf'                         # the final model, k=16 lambda=5
SWEEP_K_DIR = f'{NMF_DIR}/_archive/inmf_other_k'     # k = 6..18 at lambda=5
SWEEP_LAM_DIR = f'{NMF_DIR}/_archive/inmf_lambda_sweep'  # lambda = 0.1..2 at k=16
FINAL_DATA = f'{NMF_DIR}/final'                      # written by notebook 01 / 04
RESULTS = f'{FB}/results/final'
PLOTS = f'{FB}/plots/final'

# Xenium inputs (built by Spatial_analyses/.../final/01-02)
SPATIAL = f'{REPO}/Spatial_analyses'
PB_DIR = f'{SPATIAL}/data/xenium/pseudobulk_final'
XENIUM_COMPOSITION_GLOB = f'{SPATIAL}/results/xenium/final/02_pseudobulks/composition_by_niche_*.csv'

K, LAMBDA = 16, 5
TAG = f'k{K}_lam{LAMBDA:g}'
FACTORS = [f'S{i}' for i in range(1, K + 1)]
ARMS = ['cells', 'nuclei']
ARM_NAMES = {'cells': 'scRNA-seq (cells)', 'nuclei': 'snRNA-seq (multiome nuclei)'}


def out_dirs(section: str) -> Tuple[str, str]:
    """results/final/<section> and plots/final/<section>, created."""
    res, plt_ = f'{RESULTS}/{section}', f'{PLOTS}/{section}'
    for d in (res, plt_, FINAL_DATA):
        os.makedirs(d, exist_ok=True)
    return res, plt_


def light_path() -> str:
    """The light (obs + UMAP, no X) object: notebook 01's copy, else the original."""
    p = f'{FINAL_DATA}/fb_light.h5ad'
    if os.path.exists(p):
        return p
    warnings.warn('final/fb_light.h5ad not found (run notebook 01); using data/nmf_joint/fb_light.h5ad')
    return f'{NMF_DIR}/fb_light.h5ad'


# ---------------------------------------------------------------------------
# 2. factor vocabulary
# ---------------------------------------------------------------------------
# Curated display names (from the factor atlas, notebook 02). The data-driven panel call is kept
# next to them in the annotation table; these are the names used on every figure.
FACTOR_NAMES: Dict[str, str] = {
    'S1': 'TGFb / hypoxia / mechano', 'S2': 'lymphoid (contamination)', 'S3': 'interlobular (CCN5)',
    'S4': 'proliferation', 'S5': 'perilobular', 'S6': 'IEG / dissociation',
    'S7': 'inflammatory', 'S8': 'perilobular (nuclear)', 'S9': 'interlobular (PI16)',
    'S10': 'ECM synthesis', 'S11': 'hormone / vascular tone', 'S12': 'inflammatory stress',
    'S13': 'myofibroblast / contractile', 'S14': 'medullary', 'S15': 'mature adipose',
    'S16': 'modality (nuclear baseline)',
}
FACTOR_LABELS: Dict[str, str] = {f: f'{f} {n}' for f, n in FACTOR_NAMES.items()}

# The four activity programmes selected for the multiome accessibility analysis (notebook 04)
PROGRAMMES: Dict[str, str] = {'TGFb': 'S1', 'Inflammatory': 'S7',
                              'Inflammatory_stress': 'S12', 'Adipose_mature': 'S15'}
PROG_COLORS: Dict[str, str] = {'TGFb': '#785EF0', 'Inflammatory': '#DC267F',
                               'Inflammatory_stress': '#FE6100', 'Adipose_mature': '#FFB000',
                               'unassigned': '#a8a8a8'}

# Technical gene families excluded from the factorisation basis
TECH_PREFIX = ('MT-', 'MT1', 'MT2', 'RPS', 'RPL', 'MRPS', 'MRPL')
TECH_EXACT = {'MALAT1', 'NEAT1', 'MTRNR2L12', 'MTRNR2L8', 'ACTB', 'ACTG1',
              'GAPDH', 'TMSB4X', 'TMSB10', 'EEF1A1', 'FTL', 'FTH1', 'B2M'}


def is_technical(gene: str) -> bool:
    return gene.startswith(TECH_PREFIX) or gene in TECH_EXACT


def natural_factor_order(factors: Sequence[str]) -> List[str]:
    """S1, S2, ... S16 -- numeric, not lexicographic."""
    def key(f):
        d = ''.join(ch for ch in str(f) if ch.isdigit())
        return (0, int(d)) if d else (1, 0)
    return sorted(factors, key=key)


# ---------------------------------------------------------------------------
# 3. gene signatures (target / nuisance / split / spatial panels)
# ---------------------------------------------------------------------------
def load_signatures(path: str = SIGNATURES_CSV, verbose: bool = True) -> Dict[str, Dict[str, List[str]]]:
    """{role: {panel: [genes]}} from the editable long-format CSV.

    role = target (a programme the decomposition must recover), nuisance (absorbs a factor but
    is not a target), split (arm-specific adipose sub-panels), spatial (the Xenium gene sets,
    scored separately because they overlap each other and the identity panels).
    """
    out: Dict[str, Dict[str, List[str]]] = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            role, panel, gene = [(row.get(k) or '').strip() for k in ('role', 'panel', 'gene')]
            if not role or not panel or not gene or gene.startswith('#'):
                continue
            genes = out.setdefault(role, {}).setdefault(panel, [])
            if gene not in genes:
                genes.append(gene)
    if verbose:
        for role, d in out.items():
            print(f'{role:9s} {len(d):2d} panels, {sum(len(v) for v in d.values()):3d} genes | {", ".join(d)}')
    return out


# Which annotated state (ORIGINAL per-arm annotation: taa_l4 for cells, celltype_v3_L3 for nuclei)
# each target programme is expected to concentrate in. `None` = the arm cannot support the target.
EXPECTED_TOP_STATE = {
    'id_perilobular':  {'cells': ['Fb-perilo', 'Fb-perilo_COLEC11'], 'nuclei': ['FB_perilo', 'FB_perilo_COLEC11']},
    'id_interlobular': {'cells': ['Fb-interlo', 'Fb-interlo_PI16'], 'nuclei': ['FB_interlo']},
    'id_medullary':    {'cells': ['Fb-med'], 'nuclei': None},
    'inflammatory':    {'cells': ['Fb-interm'], 'nuclei': ['FB_interm']},
    'adipose_mature':  {'cells': ['Fb-adipo'], 'nuclei': ['FB_adipo']},
    # neither fibrotic programme has a dedicated label in the cells arm -- permissive on purpose
    'ECM_synthesis':      {'cells': ['Fb-interlo', 'Fb-interlo_PI16', 'Fb-perilo', 'Fb-interm'], 'nuclei': ['FB_TGFB1', 'FB_interlo']},
    'myofib_contractile': {'cells': ['Fb-interlo', 'Fb-interm', 'Fb-perilo', 'Fb-perilo_COLEC11'], 'nuclei': ['FB_TGFB1']},
}
# What is jointly recoverable, from the label counts of the object (28 FB_med nuclei; nuclei carry
# the transcribed adipose core but not the terminal lipid-droplet markers)
EXPECTED_ARM_SUPPORT = {
    'id_perilobular': 'both', 'id_interlobular': 'both', 'id_medullary': 'cells_only',
    'inflammatory': 'both', 'adipose_mature': 'both_core_only',
    'ECM_synthesis': 'both', 'myofib_contractile': 'both',
}
# identity programmes must concentrate in one state; activity programmes are used at graded
# intensity across states and are only required to be panel-specific and arm-balanced
PROGRAM_CLASS = {
    'id_perilobular': 'identity', 'id_interlobular': 'identity', 'id_medullary': 'identity',
    'adipose_mature': 'identity', 'inflammatory': 'activity', 'ECM_synthesis': 'activity',
    'myofib_contractile': 'activity',
}
# Fib_TGFb is a target panel in the CSV (it names S1) but was never part of the k-selection
# criterion, so the survival test runs on the seven programmes above.
SURVIVAL_TARGETS = list(PROGRAM_CLASS)


# ---------------------------------------------------------------------------
# 4. scoring primitives used by more than one notebook
# ---------------------------------------------------------------------------
def usage_shares(H: pd.DataFrame) -> pd.DataFrame:
    """Row-normalise raw iNMF usage so each cell's 16 values sum to 1 (a composition)."""
    return H.div(H.sum(axis=1).replace(0, np.nan), axis=0)


def gene_zscores(W: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """(W_norm, Z): each factor normalised to sum 1, then each GENE z-scored across factors.

    Step 1 is only valid for non-negative loadings (iNMF's W always is). Z alone is noisy for
    near-zero genes -- with 16 factors the across-factor SD of an unexpressed gene is tiny -- so
    it is never used on its own to rank single genes (see `distinctiveness`).
    """
    W_norm = W.div(W.sum(axis=1).replace(0, np.nan), axis=0)
    Z = W_norm.sub(W_norm.mean(axis=0), axis=1).div(W_norm.std(axis=0).replace(0, np.nan), axis=1)
    return W_norm, Z


def distinctiveness(W_norm: pd.DataFrame) -> pd.DataFrame:
    """loading x share of that gene's total loading held by this factor.

    High only when a gene both carries mass in the factor AND is concentrated there, so it needs
    no mass gate. This is the column factors are named from and signatures are built from.
    """
    share = W_norm.div(W_norm.sum(axis=0).replace(0, np.nan), axis=1)
    return W_norm * share


def panel_zscores(W: pd.DataFrame, panels: Dict[str, Sequence[str]]) -> pd.DataFrame:
    """panels x factors z-scores: normalise factors, z-score genes across factors, average the
    panel's genes, then z-score each panel across factors so panels of different size compare."""
    _, Z = gene_zscores(W)
    rows = {}
    for name, genes in panels.items():
        g = [x for x in genes if x in Z.columns]
        if g:
            rows[name] = Z[g].mean(axis=1)
    P = pd.DataFrame(rows).T
    return P.sub(P.mean(axis=1), axis=0).div(P.std(axis=1).replace(0, np.nan), axis=0)


def eta_sq(values, groups) -> float:
    """Fraction of variance in `values` explained by the categorical `groups`."""
    v = np.asarray(values, float)
    g = np.asarray(groups)
    gm = v.mean()
    sst = float(((v - gm) ** 2).sum())
    if sst <= 0:
        return np.nan
    ssb = sum(int((g == lev).sum()) * (v[g == lev].mean() - gm) ** 2 for lev in pd.unique(g))
    return float(ssb / sst)


def z_within(values: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    """Standardise every column across cells, separately inside each group (e.g. each arm)."""
    out = values.copy().astype(float)
    g = np.asarray(groups)
    for lev in pd.unique(g):
        m = g == lev
        sub = values[m]
        out.loc[m, :] = ((sub - sub.mean()) / sub.std().replace(0, np.nan)).values
    return out


# ---------------------------------------------------------------------------
# 5. colours
# ---------------------------------------------------------------------------
COLORS = {'blue': '#648FFF', 'purple': '#785EF0', 'magenta': '#DC267F', 'orange': '#FE6100',
          'yellow': '#FFB000', 'teal': '#18868F', 'mint': '#738B79', 'pink': '#C169E8',
          'grass': '#82CD29', 'forest': '#124420'}
GREYS = ['#f7f7f7', '#cccccc', '#969696', '#636363', '#252525']
GREY_BG = GREYS[1]


def tint(color: str, l: float) -> str:
    r, g, b = to_rgb(color)
    return to_hex((1 - l + l * r, 1 - l + l * g, 1 - l + l * b))


# harmonised fibroblast populations (taa_l3_harmonized / taa_l3_cm / taa_l3_cmx vocabulary)
POP_ORDER: List[str] = ['Fb-interlo', 'Fb-perilo', 'Fb-med', 'Fb-interm', 'Fb-TGFB+', 'Fb-adipo', 'Adipocyte']
POP_COLORS: Dict[str, str] = {
    'Fb-interlo': COLORS['orange'], 'Fb-perilo': COLORS['blue'], 'Fb-med': COLORS['magenta'],
    'Fb-interm': COLORS['teal'], 'Fb-TGFB+': COLORS['purple'],
    'Fb-adipo': tint(COLORS['yellow'], 0.55), 'Adipocyte': COLORS['yellow'],
}
ARM_COLORS = {'cells': COLORS['grass'], 'nuclei': GREYS[4]}

CMAP_USAGE = 'RdPu'
CMAP_DIV = LinearSegmentedColormap.from_list('thy_bm', [COLORS['blue'], '#f7f7f7', COLORS['magenta']])

# Xenium
DONOR_AGE = {'U51': 'epaed', 'U48': 'epaed', 'Z16': 'lpaed', 'Z17': 'lpaed', 'A67': 'adult', 'A72': 'adult', 'A70': 'adult'}
DONOR_ORDER = ['U51', 'U48', 'Z16', 'Z17', 'A67', 'A72', 'A70']
AGE_ORDER = ['epaed', 'lpaed', 'adult']
AGE_TITLE = {'epaed': 'early paediatric', 'lpaed': 'late paediatric', 'adult': 'adult'}
AGE_COLORS = {'epaed': '#9ecae1', 'lpaed': '#4292c6', 'adult': '#08306b'}
# the eight H&E morphology niches, in the order used on every spatial figure
NICHE_ORDER = ['cortex_septal', 'cortex', 'medulla', 'medulla_hc', 'adipose', 'fibrosis_adipose', 'fibrosis', 'vessel']
NICHE_ORDER_L3 = ['cortex_septal', 'cortex_outer', 'cortex_dense', 'cortex', 'cortex_residual',
                  'medulla', 'medulla_hc', 'medulla_residual', 'adipose', 'adipose_edge',
                  'fibrosis', 'fibrosis_adipose', 'endothelium_vessel', 'endothelium_edge']
NICHE_ORDERS = {'morpho_niche_l2': NICHE_ORDER, 'morpho_niche_l3': NICHE_ORDER_L3}
REP_NAME = {'patches': 'H&E patches (transcripts)', 'cells': 'segmented cells (all)', 'fbcells': 'fibroblast-compartment cells'}
FB_STATES = ['Fb-perilo', 'Fb-interlo', 'Fb-interm', 'Fb-med', 'Fb-adipo', 'Fb-TGFB+']
FB_STATES_EXT = FB_STATES + ['Adipocyte']


# ---------------------------------------------------------------------------
# 6. figure style and helpers
# ---------------------------------------------------------------------------
MPLSTYLE = f'{REPO}/General_analysis/scripts/plotting/thyAgeing.mplstyle'
ARIAL_TTF = '/nfs/team205/ny1/ThymusSpatialAtlas/software/Arial.ttf'
MM = 1 / 25.4
RASTER_DPI = 600
NATURE_W = {'single': 89.0, 'one-half': 130.0, 'double': 183.0}


def set_style(screen_dpi: int = 110) -> None:
    """House style (Arial, 5-7 pt), Type-42 fonts, fixed-canvas saving, no axis rounding."""
    plt.style.use(MPLSTYLE)
    if os.path.exists(ARIAL_TTF):
        font_manager.fontManager.addfont(ARIAL_TTF)
        plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams.update({
        'figure.dpi': screen_dpi, 'savefig.dpi': RASTER_DPI, 'savefig.format': 'pdf',
        'savefig.bbox': 'standard',        # 'tight' would defeat the fixed millimetre widths
        'savefig.transparent': False, 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'axes.autolimit_mode': 'data',     # the style's round_numbers clips UMAP point clouds
        'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
        'mathtext.bf': 'Arial:bold',
    })
    import logging
    logging.getLogger('fontTools').setLevel(logging.WARNING)


def figsize_mm(width: Union[str, float] = 'double', height_mm: float = 60.0) -> Tuple[float, float]:
    """Figure size in inches from a Nature column width ('single' 89, 'one-half' 130, 'double' 183 mm)."""
    w = NATURE_W.get(width, width) if isinstance(width, str) else float(width)
    return (w * MM, height_mm * MM)


def save(fig, path: str, close: bool = False, **kwargs) -> None:
    """Save on the requested canvas (no bbox cropping) and print the physical size."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with plt.rc_context({'savefig.bbox': 'standard'}):
        fig.savefig(path, **kwargs)
    w, h = fig.get_size_inches()
    print(f'-> {path}   {w / MM:.0f} x {h / MM:.0f} mm')
    if close:
        plt.close(fig)


def char_mm(fontsize: Optional[float] = None) -> float:
    """Approximate width of one Arial character in mm at a given point size (0.60 em)."""
    fs = fontsize or plt.rcParams['ytick.labelsize']
    return 0.60 * fs / 72 * 25.4


def frameless(ax) -> None:
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def heatmap(ax, M: pd.DataFrame, cmap, vmin, vmax, xrot: int = 90, xlabels: bool = True,
            ylabels: bool = True, grid: bool = True, annot: Optional[pd.DataFrame] = None,
            annot_fontsize: float = 4.5):
    """imshow of a DataFrame with its labels, hairline cell separators and optional text."""
    im = ax.imshow(M.values, cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto', interpolation='nearest')
    ax.set_xticks(np.arange(M.shape[1]))
    ax.set_xticklabels(list(M.columns) if xlabels else [], rotation=xrot,
                       ha='center' if xrot in (0, 90) else 'right', va='top')
    ax.set_yticks(np.arange(M.shape[0]))
    ax.set_yticklabels(list(M.index) if ylabels else [])
    ax.tick_params(length=1.5, width=0.4, pad=1, top=False, right=False)
    if not xlabels:
        ax.tick_params(axis='x', length=0)
    if not ylabels:
        ax.tick_params(axis='y', length=0)
    for s in ax.spines.values():
        s.set_linewidth(0.4)
    if grid:   # solid and above the image: the style's dotted grid would chop the cell edges
        ax.set_xticks(np.arange(-.5, M.shape[1], 1), minor=True)
        ax.set_yticks(np.arange(-.5, M.shape[0], 1), minor=True)
        ax.grid(which='minor', color='w', linewidth=0.3, linestyle='-')
        ax.grid(which='major', visible=False)
        ax.set_axisbelow(False)
        ax.tick_params(which='minor', length=0)
    if annot is not None:
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                s = annot.iat[i, j]
                if isinstance(s, str) and s:
                    ax.text(j, i, s, ha='center', va='center', fontsize=annot_fontsize)
    return im


def colorbar(fig, mappable, ax, label: str, width_mm: float = 2.0, pad_mm: float = 1.5,
             height_frac: float = 0.6, ticks=None):
    """A colour bar of fixed physical width to the right of `ax`."""
    fig.canvas.draw()
    box = ax.get_position()
    fw, fh = fig.get_size_inches()
    cax = fig.add_axes([box.x1 + pad_mm * MM / fw, box.y0 + box.height * (1 - height_frac) / 2,
                        width_mm * MM / fw, box.height * height_frac])
    cb = fig.colorbar(mappable, cax=cax, ticks=ticks)
    cb.set_label(label, fontsize=plt.rcParams['axes.labelsize'], labelpad=2)
    cb.ax.tick_params(labelsize=plt.rcParams['ytick.labelsize'], length=1.5, width=0.4, pad=1)
    cb.outline.set_linewidth(0.4)
    return cb


def scatter_umap(ax, xy: np.ndarray, color, size: float = 0.6, order: Optional[np.ndarray] = None, **kw):
    """Rasterised point cloud with a small margin; the axes are stripped and made square."""
    if order is not None:
        xy = xy[order]
        color = color[order] if np.ndim(color) else color
    if len(xy) == 0:
        frameless(ax)
        return
    ax.scatter(xy[:, 0], xy[:, 1], c=color, s=size, linewidths=0, rasterized=True, **kw)
    for lo, hi, setter in [(xy[:, 0].min(), xy[:, 0].max(), ax.set_xlim), (xy[:, 1].min(), xy[:, 1].max(), ax.set_ylim)]:
        pad = 0.04 * (hi - lo)
        setter(lo - pad, hi + pad)
    frameless(ax)
    ax.set_aspect('equal')


def umap_categorical(xy: np.ndarray, values: Sequence[str], palette: Dict[str, str],
                     order: Optional[Sequence[str]] = None, title: str = '', ax=None,
                     point_size: float = 0.35, legend: bool = True, seed: int = 0,
                     mask: Optional[np.ndarray] = None, legend_ncol: int = 1):
    """One categorical UMAP panel; cells outside `mask` are grey background."""
    vals = np.asarray(values).astype(str)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize_mm('single', 70))
    rng = np.random.default_rng(seed)
    if mask is not None:
        scatter_umap(ax, xy[~mask], GREY_BG, point_size * 0.8)
        keep = np.flatnonzero(mask)
    else:
        keep = np.arange(len(xy))
    present = vals[keep]
    levels = [l for l in (order or sorted(set(present))) if l in set(present)]
    keep = keep[rng.permutation(len(keep))]
    cols = np.array([palette.get(v, GREY_BG) for v in vals])
    scatter_umap(ax, xy[keep], cols[keep], point_size)
    scatter_umap(ax, xy, 'none', 0)          # limits from the full embedding
    ax.set_title(title)
    if legend:
        handles = [Line2D([], [], marker='o', linestyle='none', markersize=2.4, markeredgewidth=0,
                          markerfacecolor=palette.get(l, GREY_BG), label=f'{l} ({int((present == l).sum()):,})')
                   for l in levels]
        ax.legend(handles=handles, loc='upper left', frameon=False, ncol=legend_ncol, handletextpad=0.3,
                  columnspacing=0.8, labelspacing=0.3, borderpad=0.2, fontsize=plt.rcParams['legend.fontsize'])
    return ax


def umap_factor_grid(xy: np.ndarray, arms: np.ndarray, U: pd.DataFrame, labels: Dict[str, str],
                     factors: Optional[Sequence[str]] = None, arm: Optional[str] = None,
                     ncols: int = 4, point_size: float = 0.35, pmax: float = 99.0,
                     width: Union[str, float] = 'double', panel_height_mm: float = 40.0,
                     cmap=CMAP_USAGE):
    """Usage share of each factor on the UMAP, one panel per factor, each scaled 0..its p99.

    `arm` colours only that modality and leaves the other as grey background -- the check that a
    programme claimed as shared paints in both arms.
    """
    factors = list(factors or U.columns)
    arms = np.asarray(arms).astype(str)
    keep = np.ones(len(xy), bool) if arm is None else (arms == arm)
    nrows = int(np.ceil(len(factors) / ncols))
    fig = plt.figure(figsize=figsize_mm(width, panel_height_mm * nrows + 4))
    gs = GridSpec(nrows, ncols, figure=fig, hspace=0.16, wspace=0.10, left=0.005, right=0.90, top=0.955, bottom=0.01)
    for i, f in enumerate(factors):
        ax = fig.add_subplot(gs[i // ncols, i % ncols])
        v = U[f].to_numpy(float)
        vmax = float(np.nanpercentile(v[keep], pmax)) or float(np.nanmax(v[keep]))
        if keep.sum() < len(v):
            scatter_umap(ax, xy[~keep], GREY_BG, point_size * 0.8)
        order = np.argsort(v[keep])
        ax.scatter(xy[keep][order, 0], xy[keep][order, 1], c=np.clip(v[keep][order], 0, vmax), s=point_size,
                   cmap=cmap, vmin=0, vmax=vmax, linewidths=0, rasterized=True)
        scatter_umap(ax, xy, 'none', 0)
        ax.set_title(labels.get(f, f))
        ax.text(0.99, 0.01, f'max {vmax:.2f}', transform=ax.transAxes, ha='right', va='bottom', fontsize=4.2, color=GREYS[3])
    cax = fig.add_axes([0.925, 0.35, 0.012, 0.30])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap), cax=cax, ticks=[0, 1])
    cb.ax.set_yticklabels(['0', 'panel max'])
    cb.set_label('usage share', fontsize=plt.rcParams['axes.labelsize'], labelpad=2)
    cb.ax.tick_params(labelsize=plt.rcParams['ytick.labelsize'], length=1.5, width=0.4, pad=1)
    cb.outline.set_linewidth(0.4)
    return fig


def usage_heatmap(U: pd.DataFrame, obs: pd.DataFrame, label_col: str, arms: Sequence[str] = ('cells', 'nuclei'),
                  arm_col: str = 'arm', factors: Optional[Sequence[str]] = None,
                  factor_labels: Optional[Dict[str, str]] = None, pop_order: Optional[Sequence[str]] = None,
                  mode: str = 'z', min_cells: int = 20, cmap=CMAP_USAGE, vlim: float = 2.0,
                  width: Union[str, float] = 'one-half', cell_mm: float = 5.0, suptitle: Optional[str] = None):
    """Factors in ROWS, populations in COLUMNS (labels on top, coloured), one panel per arm.

    mode='z'   each factor z-scored across populations within its arm -- the quantity that is
               comparable between assays (it removes the baseline each assay puts on a factor).
    mode='abs' raw mean usage share on a shared scale.
    Returns (fig, {arm: populations x factors table}).
    """
    factor_labels = factor_labels or {}
    factors = natural_factor_order(factors or U.columns)
    tabs, mats = {}, {}
    for arm in arms:
        idx = obs.index[obs[arm_col].astype(str) == arm].intersection(U.index)
        lab = obs.loc[idx, label_col].astype(str)
        t = U.loc[idx].groupby(lab.values, observed=True).mean()
        n = lab.value_counts()
        keep = [p for p in (pop_order or POP_ORDER) if p in t.index and n.get(p, 0) >= min_cells]
        keep += sorted(p for p in t.index if p not in keep and p not in ('nan', 'None') and n.get(p, 0) >= min_cells)
        tabs[arm] = t.loc[keep]
        M = t.loc[keep].T.loc[[f for f in factors if f in t.columns]]
        if mode == 'z':
            M = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1).replace(0, np.nan), axis=0)
        mats[arm] = M
    if mode == 'z':
        lo = min(float(np.nanmin(M.values)) for M in mats.values())
        vmin, vmax, clab = max(-vlim, np.floor(lo * 2) / 2), vlim, 'usage z across populations'
    else:
        vmin, vmax, clab = 0.0, max(float(np.nanmax(M.values)) for M in mats.values()), 'mean usage share'

    ncol = {a: M.shape[1] for a, M in mats.items()}
    nrow = max(M.shape[0] for M in mats.values())
    fig_w_mm = figsize_mm(width, 10)[0] / MM
    ch = char_mm()
    left_mm = 2.0 + ch * max(len(factor_labels.get(f, f)) for f in factors)
    gap_mm, right_mm, bottom_mm = 6.0, 18.0, 4.0
    top_mm = (12.0 if suptitle else 7.0) + ch * max(len(str(p)) for M in mats.values() for p in M.columns)
    cell = min(cell_mm, (fig_w_mm - left_mm - gap_mm - right_mm) / sum(ncol.values()))
    fig_h_mm = cell * nrow + top_mm + bottom_mm
    fig = plt.figure(figsize=(fig_w_mm * MM, fig_h_mm * MM))
    if suptitle:
        fig.suptitle(suptitle, fontsize=plt.rcParams['axes.titlesize'], y=1 - 1.0 / fig_h_mm, va='top')
    x, axes, im = left_mm, {}, None
    for arm in arms:
        M = mats[arm]
        ax = fig.add_axes([x / fig_w_mm, bottom_mm / fig_h_mm, cell * ncol[arm] / fig_w_mm, cell * M.shape[0] / fig_h_mm])
        im = heatmap(ax, M, cmap, vmin, vmax, xrot=90, xlabels=True, ylabels=(arm == arms[0]))
        ax.set_xticklabels(list(M.columns), rotation=90, ha='center', va='bottom')
        ax.xaxis.set_ticks_position('top')
        ax.tick_params(axis='x', top=True, bottom=False, length=1.5, width=0.4, pad=1)
        for i, pop in enumerate(M.columns):
            ax.get_xticklabels()[i].set_color(POP_COLORS.get(pop, 'k'))
        if arm == arms[0]:
            ax.set_yticklabels([factor_labels.get(f, f) for f in M.index])
        ax.set_xlabel(ARM_NAMES.get(arm, arm), fontsize=plt.rcParams['axes.labelsize'], labelpad=3)
        axes[arm] = ax
        x += cell * ncol[arm] + gap_mm
    colorbar(fig, im, axes[arms[-1]], clab, height_frac=0.55, ticks=[vmin, 0, vmax] if mode == 'z' else None)
    return fig, tabs


def niche_dotplot(frames: Dict[str, pd.DataFrame], factors: Sequence[str], niche_order: Sequence[str],
                  value: str = 'score', niche_col: str = 'group', donor_col: str = 'donor',
                  labels: Optional[Dict[str, str]] = None, centre_within_donor: bool = True,
                  row_annot: Optional[Dict[str, Dict[str, str]]] = None, vlim: Optional[float] = None,
                  cmap=CMAP_DIV, width: Union[str, float] = 'double', cell_mm: float = 6.0,
                  max_dot_mm: float = 2.4, suptitle: Optional[str] = None):
    """Factors x niches, one panel per entry of `frames` (e.g. one per age group).

    Dot COLOUR = mean score over donors (centred within donor across niches first, because sections
    differ several-fold in depth and the question is which niche is high WITHIN a section).
    Dot AREA = fraction of donors on the same side of zero as the mean -- with 3-5 sections that is
    the number that decides whether a niche effect is real. `row_annot[title][factor]` is written
    at the right edge of the row (e.g. Kendall-W stars).
    Returns (fig, {title: (mean, agreement, n_donors)}).
    """
    labels, row_annot = labels or {}, row_annot or {}
    rows = [f for f in natural_factor_order(factors)]
    cols = list(niche_order)
    tables = {}
    for title, df in frames.items():
        d = df[df['source'].isin(rows)].copy()
        d[value] = pd.to_numeric(d[value], errors='coerce')
        if centre_within_donor:
            d[value] = d.groupby(['source', donor_col], observed=True)[value].transform(lambda x: x - x.mean())
        mean = d.pivot_table(index='source', columns=niche_col, values=value, aggfunc='mean', observed=True)
        per_donor = d.pivot_table(index=['source', niche_col], columns=donor_col, values=value, observed=True)
        agree = {}
        for (src, nic), r in per_donor.iterrows():
            v = r.dropna()
            m = mean.loc[src, nic] if (src in mean.index and nic in mean.columns) else np.nan
            agree[(src, nic)] = float((np.sign(v) == np.sign(m)).mean()) if (len(v) and np.isfinite(m) and m != 0) else np.nan
        agree = pd.Series(agree).unstack() if agree else pd.DataFrame()
        cols_here = [c for c in cols if c in mean.columns]
        tables[title] = (mean.reindex(index=rows, columns=cols_here), agree.reindex(index=rows, columns=cols_here),
                         int(per_donor.notna().sum(axis=1).max()) if len(per_donor) else 0)
    v = vlim if vlim is not None else float(np.nanmax([np.nanmax(np.abs(t[0].values)) for t in tables.values()]))
    n_donors = max(t[2] for t in tables.values())

    ch = char_mm()
    left_mm = 2.0 + ch * max(len(labels.get(f, f)) for f in rows)
    right_mm, top_mm = 44.0, (10.0 if suptitle else 6.0)
    bottom_mm = 3.0 + ch * max(len(str(c)) for c in cols)
    gap_mm = 5.0
    n_pan = len(frames)
    fig_w_mm = figsize_mm(width, 10)[0] / MM
    cell = min(cell_mm, (fig_w_mm - left_mm - right_mm - gap_mm * (n_pan - 1)) / max(len(cols) * n_pan, 1))
    fig_h_mm = cell * len(rows) + top_mm + bottom_mm
    fig = plt.figure(figsize=(fig_w_mm * MM, fig_h_mm * MM))
    max_pt = (min(max_dot_mm, cell * 0.85) / 25.4 * 72) ** 2
    x, sc_, ax_last = left_mm, None, None
    for k, (title, (mean, agree, nd)) in enumerate(tables.items()):
        ax = fig.add_axes([x / fig_w_mm, bottom_mm / fig_h_mm, cell * len(cols) / fig_w_mm, cell * len(rows) / fig_h_mm])
        xs, ys, cs, ss = [], [], [], []
        for i, f in enumerate(rows):
            for j, c in enumerate(cols):
                if c not in mean.columns:
                    continue
                a, m = agree.loc[f, c], mean.loc[f, c]
                if not (np.isfinite(a) and np.isfinite(m)):
                    continue
                xs.append(j); ys.append(i); cs.append(m); ss.append(max_pt * a)
        sc_ = ax.scatter(xs, ys, c=cs, s=ss, cmap=cmap, vmin=-v, vmax=v, linewidths=0.15, edgecolors='0.35')
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels(cols, rotation=90, ha='center', va='top')
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([labels.get(f, f) for f in rows] if k == 0 else [])
        ax.set_xlim(-0.6, len(cols) - 0.4); ax.set_ylim(len(rows) - 0.4, -0.6)
        ax.grid(False)
        ax.tick_params(length=1.5, width=0.4, pad=1, top=False, right=False)
        if k > 0:
            ax.tick_params(axis='y', length=0)
        for sp in ax.spines.values():
            sp.set_linewidth(0.4)
        ax.set_title(f'{title} (n = {nd})' if n_pan > 1 or nd else title)
        for i, f in enumerate(rows):
            s = row_annot.get(title, {}).get(f, '')
            if s:
                ax.text(len(cols) - 0.45, i, s, ha='left', va='center', fontsize=4.5, color=GREYS[3])
        x += cell * len(cols) + gap_mm
        ax_last = ax
    ax_right_mm = x - gap_mm
    ax_h_mm = cell * len(rows)
    cbar_h_mm = max(16.0, ax_h_mm * 0.55)
    cax = fig.add_axes([(ax_right_mm + 4.0) / fig_w_mm, (bottom_mm + (ax_h_mm - cbar_h_mm) / 2) / fig_h_mm, 2.0 / fig_w_mm, cbar_h_mm / fig_h_mm])
    cb = fig.colorbar(sc_, cax=cax, ticks=[-v, 0, v])
    cb.set_label('mean score\n(centred within donor)', fontsize=plt.rcParams['axes.labelsize'], labelpad=2)
    cb.ax.tick_params(labelsize=plt.rcParams['ytick.labelsize'], length=1.5, width=0.4, pad=1)
    cb.ax.set_yticklabels([f'{-v:.1f}', '0', f'{v:.1f}'])
    cb.outline.set_linewidth(0.4)
    lax = fig.add_axes([(ax_right_mm + 22.0) / fig_w_mm, bottom_mm / fig_h_mm, 20.0 / fig_w_mm, max(ax_h_mm, 18.0) / fig_h_mm])
    lax.set_axis_off()
    ks = sorted({int(round(k)) for k in np.linspace(1, max(n_donors, 1), min(max(n_donors, 1), 4))})
    handles = [Line2D([], [], marker='o', linestyle='none', markersize=np.sqrt(max_pt * k / max(n_donors, 1)),
                      markerfacecolor='0.55', markeredgewidth=0, label=f'{k / max(n_donors, 1):.0%}') for k in ks]
    lax.legend(handles=handles, loc='center left', frameon=False, title=f'donors agreeing\n(n = {n_donors})',
               title_fontsize=plt.rcParams['axes.labelsize'], fontsize=plt.rcParams['legend.fontsize'],
               handletextpad=0.5, labelspacing=0.8, borderpad=0, borderaxespad=0)
    if suptitle:
        fig.suptitle(suptitle, fontsize=plt.rcParams['axes.titlesize'], y=1 - 0.8 / fig_h_mm, va='top')
    return fig, tables


def niche_heatmap(M: pd.DataFrame, labels: Optional[Dict[str, str]] = None, cmap=CMAP_DIV,
                  vlim: Optional[float] = None, cbar_label: str = 'mean score\n(centred within donor)',
                  width: Union[str, float] = 'single', cell_mm: float = 7.0, suptitle: Optional[str] = None,
                  annot: Optional[pd.DataFrame] = None):
    """Factors (rows) x niches (columns), square cells, diverging around 0."""
    labels = labels or {}
    v = float(np.nanmax(np.abs(M.values))) if vlim is None else float(vlim)
    ch = char_mm()
    right_mm, top_mm = 16.0, (11.0 if suptitle else 6.0)
    left_mm = 2.0 + ch * max(len(labels.get(f, f)) for f in M.index)
    bottom_mm = 3.0 + ch * max(len(str(c)) for c in M.columns)
    fig_w_mm = figsize_mm(width, 10)[0] / MM
    cell = min(cell_mm, (fig_w_mm - left_mm - right_mm) / max(M.shape[1], 1))
    fig_h_mm = cell * M.shape[0] + top_mm + bottom_mm
    fig = plt.figure(figsize=(fig_w_mm * MM, fig_h_mm * MM))
    ax = fig.add_axes([left_mm / fig_w_mm, bottom_mm / fig_h_mm, cell * M.shape[1] / fig_w_mm, cell * M.shape[0] / fig_h_mm])
    im = heatmap(ax, M, cmap, -v, v, annot=annot)
    ax.set_yticklabels([labels.get(f, f) for f in M.index])
    cax = fig.add_axes([(left_mm + cell * M.shape[1] + 3) / fig_w_mm, bottom_mm / fig_h_mm, 2.0 / fig_w_mm, cell * M.shape[0] / fig_h_mm])
    cb = fig.colorbar(im, cax=cax, ticks=[-v, 0, v])
    cb.set_label(cbar_label, fontsize=plt.rcParams['axes.labelsize'], labelpad=2)
    cb.ax.tick_params(labelsize=plt.rcParams['ytick.labelsize'], length=1.5, width=0.4, pad=1)
    cb.ax.set_yticklabels([f'{-v:.1f}', '0', f'{v:.1f}'])
    cb.outline.set_linewidth(0.4)
    if suptitle:
        fig.suptitle(suptitle, fontsize=plt.rcParams['axes.titlesize'], y=1 - 0.8 / fig_h_mm, va='top')
    return fig
