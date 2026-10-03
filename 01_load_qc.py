import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import glob, os, re, gc, pickle

sc.settings.verbosity = 1

DATA_DIR = "/mnt/user-data/uploads/Downloads"
OUT_DIR = "/tmp/claude-0/-home-claude/e7962fd6-be13-56c6-9438-14c50d5a6310/scratchpad/endometriosis"

COND_MAP = {
    "Ctrl": "Control",
    "EuE": "Eutopic_Endometrium",
    "EcP": "Ectopic_Peritoneal",
    "EcPA": "Ectopic_Peritoneal_Adjacent",
    "EcPA2": "Ectopic_Peritoneal_Adjacent",
    "EcO": "Ectopic_Ovary",
}

with open(os.path.join(OUT_DIR, "common_genes.pkl"), "rb") as fh:
    REF_GENES = pd.Index(sorted(pickle.load(fh)))
print("Using fixed common gene reference:", len(REF_GENES), flush=True)

files = sorted(glob.glob(os.path.join(DATA_DIR, "GSM*_filtered_feature_bc_matrix.h5")))
print(f"Found {len(files)} files", flush=True)

meta_rows = []
running = None

for i, f in enumerate(files):
    base = os.path.basename(f)
    m = re.match(r"(GSM\d+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)_filtered_feature_bc_matrix\.h5", base)
    if not m:
        print("WARN: could not parse", base)
        continue
    gsm, patient_code, cond_code = m.groups()
    cond = COND_MAP.get(cond_code, cond_code)
    sample_id = f"{patient_code}_{cond_code}"

    a = sc.read_10x_h5(f)
    a.var_names_make_unique()
    a.X = a.X.astype(np.float32)
    a.obs["sample_id"] = sample_id
    a.obs["gsm"] = gsm
    a.obs["patient_code"] = patient_code
    a.obs["condition_code"] = cond_code
    a.obs["condition"] = cond
    a.obs_names = [f"{sample_id}_{bc}" for bc in a.obs_names]

    sc.pp.filter_cells(a, min_genes=200)
    a.var = a.var[[]]
    # subset THIS SMALL sample only, onto the precomputed fixed reference gene list
    a = a[:, REF_GENES].copy()

    meta_rows.append({"gsm": gsm, "sample_id": sample_id, "patient_code": patient_code,
                       "condition_code": cond_code, "condition": cond, "n_cells_after_filter": a.n_obs})
    print(f"[{i+1}/{len(files)}] {sample_id} {cond} shape={a.shape}", flush=True)

    if running is None:
        running = a
    else:
        running = ad.concat([running, a], join="inner", index_unique=None)
        del a
    gc.collect()
    print(f"   running shape={running.shape}", flush=True)

meta = pd.DataFrame(meta_rows)
meta.to_csv(os.path.join(OUT_DIR, "sample_metadata.csv"), index=False)

adata = running
del running
gc.collect()
print("Combined shape (pre-QC):", adata.shape, flush=True)

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

print("Before final QC:", adata.shape, flush=True)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.pct_counts_mt < 20].copy()
adata = adata[adata.obs.n_genes_by_counts < 8000].copy()
print("After final QC:", adata.shape, flush=True)

print(adata.obs.groupby(["condition"]).size(), flush=True)

adata.write(os.path.join(OUT_DIR, "adata_qc.h5ad"))
print("Saved adata_qc.h5ad", flush=True)
