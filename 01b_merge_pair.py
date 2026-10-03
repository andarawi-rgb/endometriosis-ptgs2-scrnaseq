import sys, gc
import scanpy as sc
import anndata as ad

f1, f2, out = sys.argv[1], sys.argv[2], sys.argv[3]
print("Loading", f1, flush=True)
a1 = sc.read_h5ad(f1)
print("shape1", a1.shape, flush=True)
print("Loading", f2, flush=True)
a2 = sc.read_h5ad(f2)
print("shape2", a2.shape, flush=True)
print("Concatenating...", flush=True)
merged = ad.concat([a1, a2], join="inner", index_unique=None)
del a1, a2
gc.collect()
print("merged shape", merged.shape, flush=True)
merged.write(out)
print("Saved", out, flush=True)
