import scanpy as sc
import numpy as np
import gc, resource

def mem():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024

OUT_DIR = "/tmp/claude-0/-home-claude/e7962fd6-be13-56c6-9438-14c50d5a6310/scratchpad/endometriosis"

adata = sc.read_h5ad(f"{OUT_DIR}/adata_norm.h5ad")
print("Loaded:", adata.shape, "peak GB:", mem(), flush=True)

hvg_mask = adata.var["highly_variable"].to_numpy()
adata_hvg = adata[:, hvg_mask].copy()
print("HVG subset:", adata_hvg.shape, "peak GB:", mem(), flush=True)

# free the full-gene object entirely -- we do NOT need it again until a later,
# separate script re-merges results by obs_names
del adata
gc.collect()
print("full-gene object freed. peak GB:", mem(), flush=True)

sc.pp.scale(adata_hvg, max_value=10)
print("scale done. peak GB:", mem(), flush=True)

sc.tl.pca(adata_hvg, n_comps=30, svd_solver="arpack")
print("PCA done. peak GB:", mem(), flush=True)

import harmonypy
ho = harmonypy.run_harmony(adata_hvg.obsm["X_pca"], adata_hvg.obs, ["patient_code"], max_iter_harmony=10)
# harmonypy>=2.0 (C++ backend): Z_corr property already returns (N cells x d), NOT (d x N) like
# the old pure-python harmonypy 0.0.x. Do NOT transpose again here.
Z_corr = ho.Z_corr
print("Z_corr raw shape:", Z_corr.shape, flush=True)
if Z_corr.shape[0] != adata_hvg.n_obs and Z_corr.shape[1] == adata_hvg.n_obs:
    Z_corr = Z_corr.T
adata_hvg.obsm["X_pca_harmony"] = Z_corr
del ho, Z_corr
gc.collect()
print("Harmony done. peak GB:", mem(), flush=True)

sc.pp.neighbors(adata_hvg, use_rep="X_pca_harmony", n_neighbors=15)
print("neighbors done. peak GB:", mem(), flush=True)

sc.tl.leiden(adata_hvg, resolution=1.0, flavor="igraph", n_iterations=2)
print("leiden done. n_clusters:", adata_hvg.obs.leiden.nunique(), "peak GB:", mem(), flush=True)
print(adata_hvg.obs.leiden.value_counts(), flush=True)

# save only the lightweight results needed to re-merge onto the full-gene data
np.save(f"{OUT_DIR}/harmony_pca.npy", adata_hvg.obsm["X_pca_harmony"])
adata_hvg.obs[["leiden"]].to_csv(f"{OUT_DIR}/leiden_labels.csv")
print("Saved harmony_pca.npy + leiden_labels.csv. peak GB:", mem(), flush=True)
