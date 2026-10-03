import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import gc, resource

def mem():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024

OUT_DIR = "/tmp/claude-0/-home-claude/e7962fd6-be13-56c6-9438-14c50d5a6310/scratchpad/endometriosis"

adata = sc.read_h5ad(f"{OUT_DIR}/adata_annotated.h5ad")
print("Loaded:", adata.shape, "peak GB:", mem(), flush=True)

X = adata.X.tocsc()
n_cells, n_genes = X.shape
leiden = adata.obs["leiden"].astype(str).values
cell_type = adata.obs["cell_type"].values
clusters = sorted(adata.obs["leiden"].unique(), key=lambda x: int(x))
gene_names = adata.var_names.to_numpy()

# global mean per gene (dense, 27869-length vector -- cheap)
global_sum = np.asarray(X.sum(axis=0)).ravel()
global_mean = global_sum / n_cells
print("Global mean computed. peak GB:", mem(), flush=True)

rows = []
for cl in clusters:
    mask = leiden == cl
    n = mask.sum()
    Xc = X[mask]
    cl_mean = np.asarray(Xc.sum(axis=0)).ravel() / n
    # simple unbiased score: cluster mean minus global mean (log-space "delta"), proxy for specificity
    delta = cl_mean - global_mean
    top_idx = np.argsort(delta)[::-1][:8]
    top_genes = gene_names[top_idx]
    top_deltas = delta[top_idx]
    ct_call = cell_type[mask][0]
    rows.append({
        "leiden": cl, "cell_type_call": ct_call, "n_cells": int(n),
        "top_unbiased_markers": ", ".join(f"{g}(+{d:.2f})" for g, d in zip(top_genes, top_deltas)),
    })
    del Xc
    gc.collect()

df = pd.DataFrame(rows)
df.to_csv(f"{OUT_DIR}/unbiased_top_markers_per_cluster.csv", index=False)
pd.set_option("display.max_colwidth", 200)
print(df.to_string(index=False), flush=True)
print("\nDONE. peak GB:", mem(), flush=True)
