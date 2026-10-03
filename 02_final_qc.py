import scanpy as sc
import numpy as np
import gc, resource

def mem():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024

OUT_DIR = "/tmp/claude-0/-home-claude/e7962fd6-be13-56c6-9438-14c50d5a6310/scratchpad/endometriosis"

adata = sc.read_h5ad(f"{OUT_DIR}/adata_combined.h5ad")
print("Loaded:", adata.shape, "peak GB:", mem(), flush=True)

X = adata.X  # CSR sparse, float32
n_cells, n_genes = X.shape

# n_genes_by_counts via indptr diff -- zero extra memory
n_genes_by_counts = np.diff(X.indptr)
print("n_genes_by_counts done. peak GB:", mem(), flush=True)

total_counts = np.asarray(X.sum(axis=1)).ravel()
print("total_counts done. peak GB:", mem(), flush=True)

mt_mask = np.array([g.startswith("MT-") for g in adata.var_names])
mt_idx = np.where(mt_mask)[0]
print("n MT genes:", len(mt_idx), flush=True)
mt_counts = np.asarray(X[:, mt_idx].sum(axis=1)).ravel()
print("mt_counts done. peak GB:", mem(), flush=True)

pct_counts_mt = np.divide(mt_counts, total_counts, out=np.zeros_like(mt_counts), where=total_counts > 0) * 100
adata.obs["n_genes_by_counts"] = n_genes_by_counts
adata.obs["total_counts"] = total_counts
adata.obs["pct_counts_mt"] = pct_counts_mt
del mt_counts
gc.collect()
print("QC metrics assigned. peak GB:", mem(), flush=True)

# gene nnz counts via bincount on existing indices array -- no copy of X
gene_nnz = np.bincount(X.indices, minlength=n_genes)
print("gene_nnz done. peak GB:", mem(), flush=True)
gene_keep = gene_nnz >= 3
print(f"Genes kept: {gene_keep.sum()} / {n_genes}", flush=True)

cell_keep = (pct_counts_mt < 20) & (n_genes_by_counts < 8000)
print(f"Cells kept: {cell_keep.sum()} / {n_cells}", flush=True)

np.save(f"{OUT_DIR}/cell_keep.npy", cell_keep)
np.save(f"{OUT_DIR}/gene_keep.npy", gene_keep)
adata.obs.to_csv(f"{OUT_DIR}/obs_with_qc.csv")
print("Saved masks + obs to disk. peak GB:", mem(), flush=True)
print("DONE (subset step deferred to next script to control memory)", flush=True)
