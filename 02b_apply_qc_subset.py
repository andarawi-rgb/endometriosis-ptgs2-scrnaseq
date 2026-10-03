import scanpy as sc
import numpy as np
import gc, resource

def mem():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024

OUT_DIR = "/tmp/claude-0/-home-claude/e7962fd6-be13-56c6-9438-14c50d5a6310/scratchpad/endometriosis"

cell_keep = np.load(f"{OUT_DIR}/cell_keep.npy")
gene_keep = np.load(f"{OUT_DIR}/gene_keep.npy")
print("masks loaded", cell_keep.sum(), gene_keep.sum(), flush=True)

adata = sc.read_h5ad(f"{OUT_DIR}/adata_combined.h5ad")
print("Loaded:", adata.shape, "peak GB:", mem(), flush=True)

adata = adata[cell_keep, :].copy()
gc.collect()
print("after cell subset (materialized):", adata.shape, "peak GB:", mem(), flush=True)

adata.X = adata.X.tocsc()
gc.collect()
print("converted to CSC. peak GB:", mem(), flush=True)

adata = adata[:, gene_keep]
adata.X = adata.X.tocsr()
gc.collect()
print("after gene subset (CSC-sliced, back to CSR):", adata.shape, "peak GB:", mem(), flush=True)

adata = adata.copy()
gc.collect()
print("materialized final copy. peak GB:", mem(), flush=True)

adata.write(f"{OUT_DIR}/adata_qc.h5ad")
print("Saved adata_qc.h5ad. peak GB:", mem(), flush=True)
