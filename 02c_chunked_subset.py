import scanpy as sc
import anndata as ad
import numpy as np
import scipy.sparse as sp
import h5py
import gc, resource

def mem():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024

OUT_DIR = "/tmp/claude-0/-home-claude/e7962fd6-be13-56c6-9438-14c50d5a6310/scratchpad/endometriosis"

cell_keep = np.load(f"{OUT_DIR}/cell_keep.npy")
gene_keep = np.load(f"{OUT_DIR}/gene_keep.npy")
gene_idx_keep = np.where(gene_keep)[0]
print("masks loaded", cell_keep.sum(), gene_keep.sum(), flush=True)

# read obs/var (small) and X in backed mode (does NOT load the sparse matrix into RAM)
adata_backed = sc.read_h5ad(f"{OUT_DIR}/adata_combined.h5ad", backed="r")
print("Backed open:", adata_backed.shape, "peak GB:", mem(), flush=True)

n_cells = adata_backed.n_obs
CHUNK = 8000

data_parts, indices_parts, indptr_list = [], [], [0]
kept_obs_idx = []

for start in range(0, n_cells, CHUNK):
    end = min(start + CHUNK, n_cells)
    block_keep = cell_keep[start:end]
    if not block_keep.any():
        continue
    X_block = adata_backed.X[start:end]  # reads only this row range from disk
    if sp.issparse(X_block):
        X_block = X_block.tocsr()
    else:
        X_block = sp.csr_matrix(X_block)
    X_block = X_block[block_keep][:, gene_idx_keep]
    X_block = X_block.tocsr()

    data_parts.append(X_block.data.astype(np.float32))
    indices_parts.append(X_block.indices.astype(np.int32))
    for i in range(X_block.shape[0]):
        indptr_list.append(indptr_list[-1] + (X_block.indptr[i+1] - X_block.indptr[i]))

    kept_obs_idx.extend(list(np.where(block_keep)[0] + start))
    print(f"  chunk {start}-{end}: kept {X_block.shape[0]} cells, peak GB: {mem():.2f}", flush=True)
    del X_block
    gc.collect()

print("Assembling final sparse matrix...", flush=True)
data = np.concatenate(data_parts)
indices = np.concatenate(indices_parts)
indptr = np.array(indptr_list, dtype=np.int64)
del data_parts, indices_parts
gc.collect()

n_final_cells = len(kept_obs_idx)
n_final_genes = len(gene_idx_keep)
X_final = sp.csr_matrix((data, indices, indptr), shape=(n_final_cells, n_final_genes))
print("Final X:", X_final.shape, X_final.nnz, "peak GB:", mem(), flush=True)

obs_full = adata_backed.obs
var_full = adata_backed.var

obs_final = obs_full.iloc[kept_obs_idx].copy()
var_final = var_full.iloc[gene_idx_keep].copy()

adata_final = ad.AnnData(X=X_final, obs=obs_final, var=var_final)
print("Final AnnData:", adata_final.shape, "peak GB:", mem(), flush=True)

adata_final.write(f"{OUT_DIR}/adata_qc.h5ad")
print("Saved adata_qc.h5ad. peak GB:", mem(), flush=True)
