import sys, os, re, gc, pickle
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

sc.settings.verbosity = 0

DATA_DIR = "/mnt/user-data/uploads/Downloads"
OUT_DIR = "/tmp/claude-0/-home-claude/e7962fd6-be13-56c6-9438-14c50d5a6310/scratchpad/endometriosis"

COND_MAP = {
    "Ctrl": "Control", "EuE": "Eutopic_Endometrium", "EcP": "Ectopic_Peritoneal",
    "EcPA": "Ectopic_Peritoneal_Adjacent", "EcPA2": "Ectopic_Peritoneal_Adjacent",
    "EcO": "Ectopic_Ovary",
}

group_idx = int(sys.argv[1])
n_groups = int(sys.argv[2])

with open(os.path.join(OUT_DIR, "common_genes.pkl"), "rb") as fh:
    REF_GENES = pd.Index(sorted(pickle.load(fh)))

files = sorted(__import__("glob").glob(os.path.join(DATA_DIR, "GSM*_filtered_feature_bc_matrix.h5")))
my_files = [f for i, f in enumerate(files) if i % n_groups == group_idx]
print(f"Group {group_idx}: {len(my_files)} files", flush=True)

running = None
meta_rows = []
for f in my_files:
    base = os.path.basename(f)
    m = re.match(r"(GSM\d+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)_filtered_feature_bc_matrix\.h5", base)
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
    a = a[:, REF_GENES].copy()

    meta_rows.append({"gsm": gsm, "sample_id": sample_id, "patient_code": patient_code,
                       "condition_code": cond_code, "condition": cond, "n_cells_after_filter": a.n_obs})
    print(f"  {sample_id} {cond} shape={a.shape}", flush=True)

    if running is None:
        running = a
    else:
        running = ad.concat([running, a], join="inner", index_unique=None)
        del a
    gc.collect()

print(f"Group {group_idx} final shape: {running.shape}", flush=True)
running.write(os.path.join(OUT_DIR, f"group_{group_idx}.h5ad"))
pd.DataFrame(meta_rows).to_csv(os.path.join(OUT_DIR, f"group_{group_idx}_meta.csv"), index=False)
print("Saved.", flush=True)
