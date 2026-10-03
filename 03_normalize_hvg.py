import scanpy as sc
import numpy as np
import gc, resource

def mem():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024

OUT_DIR = "/tmp/claude-0/-home-claude/e7962fd6-be13-56c6-9438-14c50d5a6310/scratchpad/endometriosis"

adata = sc.read_h5ad(f"{OUT_DIR}/adata_qc.h5ad")
print("Loaded:", adata.shape, "peak GB:", mem(), flush=True)

sc.pp.normalize_total(adata, target_sum=1e4)
print("normalize_total done. peak GB:", mem(), flush=True)
sc.pp.log1p(adata)
print("log1p done. peak GB:", mem(), flush=True)

# HVG selection (does not alter X by default; adds adata.var['highly_variable'])
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample_id")
print("HVG done. n_hvg:", adata.var.highly_variable.sum(), "peak GB:", mem(), flush=True)

adata.write(f"{OUT_DIR}/adata_norm.h5ad")
print("Saved adata_norm.h5ad. peak GB:", mem(), flush=True)
