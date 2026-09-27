"""Loaders, alignment guards, pathway sets and a small file cache for the Fig. 1 exploration."""

import json
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

import config as C


# ------------------------------------------------------------------ hydration
def hydrate(path):
    """Force-download a Dropbox online-only placeholder (full size, zero blocks) before reading it."""
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    try:
        if os.stat(path).st_blocks != 0:
            return path
    except AttributeError:
        return path
    warnings.warn(f"hydrating Dropbox placeholder: {path}")
    with open(path, "rb") as fh:
        while fh.read(1 << 24):
            pass
    return path


# ------------------------------------------------------------------ gene sets
def load_panel():
    return np.load(hydrate(C.UTILS_DIR + "shared_genes.npy"), allow_pickle=True).astype(str)


def read_gmt(path):
    """GMT: name <tab> source <tab> gene1 <tab> gene2 ... -> dict name -> (source, [genes])."""
    out = {}
    with open(hydrate(path)) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            out[parts[0]] = (parts[1], [g for g in parts[2:] if g])
    return out


def load_pathways(genes, min_genes=None, extra_sets=None):
    """Gene x pathway indicator (csr) restricted to the panel, plus names, sources and sizes."""
    min_genes = min_genes or C.PATHWAY_MIN_GENES
    sets = read_gmt(C.UTILS_DIR + C.GMT_FILE)
    if extra_sets:
        for k, v in extra_sets.items():
            sets[k] = v
    pos = {g: i for i, g in enumerate(genes)}
    rows, cols, names, sources, sizes = [], [], [], [], []
    for name, (src, members) in sets.items():
        idx = [pos[g] for g in members if g in pos]
        if len(idx) < min_genes:
            continue
        j = len(names)
        rows.extend(idx)
        cols.extend([j] * len(idx))
        names.append(name)
        sources.append(src)
        sizes.append(len(idx))
    M = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(genes), len(names)))
    info = pd.DataFrame({"pathway": names, "source": sources, "n_genes": sizes})
    return M, info


def load_sg_genes(panel):
    """Stress-granule marker genes (fraction in SGs > SG_THR) present on the panel."""
    df = pd.read_excel(hydrate(C.UTILS_DIR + "SG_markers.xlsx"))
    frac = pd.to_numeric(df["Fraction of RNA molecules in SGs"], errors="coerce")
    sg = df.loc[frac > C.SG_THR, "gene"].astype(str)
    return sorted(set(sg) & set(panel))


# ------------------------------------------------------------------ per dataset
def load_tumor_obs(ds):
    """Tumor-cell obs for one dataset, in cell_ids.npy row order (asserted)."""
    p = f"{C.DATA_DIR}{ds}/"
    cell_ids = np.load(hydrate(p + "processed_data/cell_ids.npy"), allow_pickle=True).astype(str)
    adata = sc.read_h5ad(hydrate(p + "intermediate_data/adata.h5ad"), backed="r")
    obs = adata.obs
    is_tum = obs["cell_type_merged"].astype(str).values == C.TUMOR_TYPE
    ot = obs.loc[is_tum].copy()
    if list(ot["cell_id"].astype(str)) != list(cell_ids):
        raise ValueError(f"{ds}: cell_id order mismatch between adata and cell_ids.npy")
    adata.file.close()
    ot = ot.reset_index(drop=True)
    ot["cell_id"] = ot["cell_id"].astype(str)
    for col in ("cell_type_merged", "segmentation_method"):
        if col in ot:
            ot[col] = ot[col].astype(str)
    return ot, cell_ids


def load_all_obs(ds):
    """All cells of the slide: cell_id, coordinates, merged cell type (for neighbourhood geometry)."""
    p = f"{C.DATA_DIR}{ds}/"
    adata = sc.read_h5ad(hydrate(p + "intermediate_data/adata.h5ad"), backed="r")
    obs = adata.obs[["cell_id", "global_x", "global_y", "cell_type_merged"]].copy()
    adata.file.close()
    obs["cell_id"] = obs["cell_id"].astype(str)
    obs["cell_type_merged"] = obs["cell_type_merged"].astype(str)
    return obs.reset_index(drop=True)


def load_compartments(ds, panel=None):
    """Nuclear and cytoplasmic count matrices (tumor cells x panel genes), CSR int32."""
    p = f"{C.DATA_DIR}{ds}/processed_data/"
    Xn = sparse.load_npz(hydrate(p + "nuclear_expression_matrix.npz")).tocsr().astype(np.int32)
    Xc = sparse.load_npz(hydrate(p + "cytoplasmic_expression_matrix.npz")).tocsr().astype(np.int32)
    genes = np.load(hydrate(p + "gene_ids.npy"), allow_pickle=True).astype(str)
    if Xn.shape != Xc.shape:
        raise ValueError(f"{ds}: nuclear {Xn.shape} vs cytoplasmic {Xc.shape}")
    if Xn.shape[1] != len(genes):
        raise ValueError(f"{ds}: {Xn.shape[1]} columns vs {len(genes)} gene_ids")
    if panel is not None and list(genes) != list(panel):
        raise ValueError(f"{ds}: gene_ids.npy does not match shared_genes.npy")
    Xn.sum_duplicates(); Xc.sum_duplicates()
    Xn.eliminate_zeros(); Xc.eliminate_zeros()
    return Xn, Xc, genes


# ------------------------------------------------------------------ output paths and cache
def out_dir(ds=None, sub=None):
    d = C.OUT_DIR if ds is None else C.OUT_DIR + ds + "/"
    if sub:
        d = d + sub + "/"
    os.makedirs(d, exist_ok=True)
    return d


def save(path, obj):
    ext = os.path.splitext(path)[1]
    if ext == ".parquet":
        obj.to_parquet(path, index=False)
    elif ext == ".csv":
        obj.to_csv(path, index=False)
    elif ext == ".json":
        with open(path, "w") as fh:
            json.dump(obj, fh, indent=2, default=_json_default)
    elif ext == ".npz":
        if sparse.issparse(obj):
            sparse.save_npz(path, obj.tocsr())
        else:
            np.savez_compressed(path, **obj)
    elif ext == ".npy":
        np.save(path, obj)
    else:
        raise ValueError(f"unknown extension: {path}")


def load(path):
    ext = os.path.splitext(path)[1]
    if ext == ".parquet":
        return pd.read_parquet(path)
    if ext == ".csv":
        return pd.read_csv(path)
    if ext == ".json":
        with open(path) as fh:
            return json.load(fh)
    if ext == ".npz":
        try:
            return sparse.load_npz(path)
        except Exception:
            with np.load(path, allow_pickle=True) as z:
                return {k: z[k] for k in z.files}
    if ext == ".npy":
        return np.load(path, allow_pickle=True)
    raise ValueError(f"unknown extension: {path}")


def cached(path, fn, force=False, verbose=True):
    """Return the cached object at `path` unless missing or `force`; otherwise compute with fn() and save."""
    if os.path.exists(path) and not force:
        if verbose:
            print(f"    [cache] {os.path.relpath(path, C.OUT_DIR)}")
        return load(path)
    obj = fn()
    save(path, obj)
    return obj


def exists(path):
    return os.path.exists(path)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (pd.Series,)):
        return o.tolist()
    raise TypeError(f"not JSON serialisable: {type(o)}")
