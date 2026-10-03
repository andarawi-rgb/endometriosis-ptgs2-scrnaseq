"""
Core scientific analysis: which cell type is the dominant source of PTGS2 (COX-2) /
PTGS1 / PTGES / PTGER2 / PTGER4 expression in the GSE179640 endometriosis atlas,
and does this differ by anatomic/disease condition.

Methodology mirrors the neuroblastoma mPGES-1 project: patient-level pseudobulk
(raw counts summed per patient x cell_type, CPM-normalized against that group's
own total library size, log2(CPM+1)), Wilcoxon signed-rank paired tests across
patients (not per-cell tests, which pseudo-replicate and inflate significance),
Benjamini-Hochberg correction for multiple comparisons, minimum-cell-count
thresholds per pseudobulk unit, and honest reporting of small-N strata.
"""
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
from scipy.stats import wilcoxon, kruskal
import gc, resource, itertools, json

def mem():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024

OUT_DIR = "/tmp/claude-0/-home-claude/e7962fd6-be13-56c6-9438-14c50d5a6310/scratchpad/endometriosis"
GENES = ["PTGS2", "PTGS1", "PTGES", "PTGER2", "PTGER4"]
MIN_CELLS_PER_UNIT = 10  # minimum cells per patient x cell_type pseudobulk unit

def bh_correct(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    return out

# ---------- load raw (post-QC, pre-normalization) counts ----------
adata = sc.read_h5ad(f"{OUT_DIR}/adata_qc.h5ad")
print("Loaded raw QC'd counts:", adata.shape, "peak GB:", mem(), flush=True)

obs_ann = pd.read_csv(f"{OUT_DIR}/obs_annotated.csv", index_col=0)
assert list(obs_ann.index) == list(adata.obs_names), "obs order mismatch"
adata.obs["leiden"] = obs_ann["leiden"].astype(str).values
adata.obs["cell_type"] = obs_ann["cell_type"].values
print("Annotation merged.", flush=True)

X = adata.X.tocsr()
total_counts = np.asarray(X.sum(axis=1)).ravel()  # per-cell library size (all 27,869 genes)
adata.obs["total_counts"] = total_counts
print("Per-cell library size computed. peak GB:", mem(), flush=True)

gene_idx = {g: adata.var_names.get_loc(g) for g in GENES if g in adata.var_names}
print("Gene columns located:", gene_idx, flush=True)

Xcsc = X.tocsc()
gene_counts = {g: np.asarray(Xcsc[:, i].todense()).ravel() for g, i in gene_idx.items()}
del Xcsc
gc.collect()
for g, arr in gene_counts.items():
    adata.obs[f"count_{g}"] = arr
print("Per-cell raw counts for target genes extracted. peak GB:", mem(), flush=True)

obs = adata.obs[["sample_id", "patient_code", "condition_code", "condition", "cell_type",
                  "total_counts"] + [f"count_{g}" for g in GENES]].copy()
del adata, X
gc.collect()
print("Full-gene object freed; working from lightweight obs table only. peak GB:", mem(), flush=True)

obs.to_csv(f"{OUT_DIR}/obs_with_gene_counts.csv")

# ============================================================
# 1) DESCRIPTIVE: % cells expressing + mean raw count, per cell type (pooled)
# ============================================================
desc_rows = []
for ct, sub in obs.groupby("cell_type", observed=True):
    row = {"cell_type": ct, "n_cells": len(sub)}
    for g in GENES:
        c = sub[f"count_{g}"]
        row[f"pct_expr_{g}"] = 100.0 * (c > 0).mean()
        row[f"mean_count_{g}"] = c.mean()
        row[f"mean_cpm_{g}"] = (c / sub["total_counts"] * 1e4).mean()  # mean per-cell CP10K, descriptive only
    desc_rows.append(row)
desc_df = pd.DataFrame(desc_rows).sort_values("pct_expr_PTGS2", ascending=False)
desc_df.to_csv(f"{OUT_DIR}/descriptive_celltype_expression.csv", index=False)
print("\n=== Descriptive: %% cells expressing / mean count per cell type (pooled, all conditions) ===", flush=True)
print(desc_df.round(3).to_string(), flush=True)

# by condition too
desc_cond_rows = []
for (ct, cond), sub in obs.groupby(["cell_type", "condition"], observed=True):
    if len(sub) < MIN_CELLS_PER_UNIT:
        continue
    row = {"cell_type": ct, "condition": cond, "n_cells": len(sub)}
    for g in GENES:
        c = sub[f"count_{g}"]
        row[f"pct_expr_{g}"] = 100.0 * (c > 0).mean()
        row[f"mean_cpm_{g}"] = (c / sub["total_counts"] * 1e4).mean()
    desc_cond_rows.append(row)
desc_cond_df = pd.DataFrame(desc_cond_rows)
desc_cond_df.to_csv(f"{OUT_DIR}/descriptive_celltype_by_condition.csv", index=False)

# ============================================================
# 2) PATIENT-LEVEL PSEUDOBULK (primary rigorous analysis)
# ============================================================
group = obs.groupby(["patient_code", "cell_type"], observed=True)
pb = group.agg(n_cells=("total_counts", "size"), libsize=("total_counts", "sum"),
                **{f"sum_{g}": (f"count_{g}", "sum") for g in GENES}).reset_index()
pb = pb[pb["n_cells"] >= MIN_CELLS_PER_UNIT].copy()
for g in GENES:
    cpm = pb[f"sum_{g}"] / pb["libsize"] * 1e6
    pb[f"log2cpm_{g}"] = np.log2(cpm + 1)
pb.to_csv(f"{OUT_DIR}/pseudobulk_patient_celltype.csv", index=False)
print(f"\nPseudobulk units (patient x cell_type, n_cells>={MIN_CELLS_PER_UNIT}): {len(pb)}", flush=True)

# omnibus Kruskal-Wallis across cell types (pooled pseudobulk units, unpaired) for PTGS2
kw_results = {}
for g in GENES:
    groups = [sub[f"log2cpm_{g}"].values for _, sub in pb.groupby("cell_type", observed=True)]
    groups = [x for x in groups if len(x) >= 3]
    if len(groups) >= 2:
        stat, p = kruskal(*groups)
        kw_results[g] = {"H": stat, "p": p, "n_groups": len(groups)}
print("\n=== Omnibus Kruskal-Wallis across cell types (pseudobulk units, unpaired) ===", flush=True)
for g, r in kw_results.items():
    print(f"  {g}: H={r['H']:.2f}, p={r['p']:.3e}, n_groups={r['n_groups']}", flush=True)

# median log2CPM per cell type per gene -> ranking table
median_table = pb.groupby("cell_type", observed=True)[[f"log2cpm_{g}" for g in GENES]].median()
median_table["n_patients"] = pb.groupby("cell_type", observed=True).size()
median_table = median_table.sort_values("log2cpm_PTGS2", ascending=False)
median_table.to_csv(f"{OUT_DIR}/pseudobulk_median_by_celltype.csv")
print("\n=== Median patient-level pseudobulk log2(CPM+1) per cell type (ranked by PTGS2) ===", flush=True)
print(median_table.round(3).to_string(), flush=True)

# ---- primary hypothesis test: PAIRED Wilcoxon, Macrophage_Myeloid vs every other cell type, for PTGS2 ----
REF_CT = "Macrophage_Myeloid"
cell_types = sorted(pb["cell_type"].unique())
pairwise_rows = []
for g in GENES:
    wide = pb.pivot(index="patient_code", columns="cell_type", values=f"log2cpm_{g}")
    for ct in cell_types:
        if ct == REF_CT:
            continue
        paired = wide[[REF_CT, ct]].dropna()
        if len(paired) < 4:
            continue
        diff = paired[REF_CT] - paired[ct]
        if (diff == 0).all():
            continue
        try:
            stat, p = wilcoxon(paired[REF_CT], paired[ct])
        except ValueError:
            continue
        pairwise_rows.append({
            "gene": g, "ref": REF_CT, "vs": ct, "n_paired_patients": len(paired),
            "median_diff_log2cpm": diff.median(), "wilcoxon_stat": stat, "p_raw": p,
        })
pairwise_df = pd.DataFrame(pairwise_rows)
if len(pairwise_df):
    pairwise_df["p_bh"] = bh_correct(pairwise_df["p_raw"].values)
    pairwise_df = pairwise_df.sort_values(["gene", "p_bh"])
pairwise_df.to_csv(f"{OUT_DIR}/pairwise_wilcoxon_macrophage_vs_others.csv", index=False)
print(f"\n=== Paired Wilcoxon: {REF_CT} vs every other cell type (patient-level pseudobulk log2CPM) ===", flush=True)
print(pairwise_df.round(4).to_string(index=False), flush=True)

# ============================================================
# 3) CONDITION-STRATIFIED (secondary / exploratory — smaller N, report cautiously)
# ============================================================
group_c = obs.groupby(["patient_code", "cell_type", "condition"], observed=True)
pbc = group_c.agg(n_cells=("total_counts", "size"), libsize=("total_counts", "sum"),
                   **{f"sum_{g}": (f"count_{g}", "sum") for g in GENES}).reset_index()
pbc = pbc[pbc["n_cells"] >= MIN_CELLS_PER_UNIT].copy()
for g in GENES:
    cpm = pbc[f"sum_{g}"] / pbc["libsize"] * 1e6
    pbc[f"log2cpm_{g}"] = np.log2(cpm + 1)
pbc.to_csv(f"{OUT_DIR}/pseudobulk_patient_celltype_condition.csv", index=False)

cond_summary = pbc.groupby(["condition", "cell_type"], observed=True).agg(
    n_patients=("patient_code", "nunique"),
    median_log2cpm_PTGS2=("log2cpm_PTGS2", "median"),
).reset_index()
cond_summary.to_csv(f"{OUT_DIR}/pseudobulk_condition_celltype_summary.csv", index=False)
print("\n=== Condition-stratified: n patients & median PTGS2 log2CPM per (condition, cell_type) ===", flush=True)
print(cond_summary.sort_values(["condition", "median_log2cpm_PTGS2"], ascending=[True, False]).round(3).to_string(index=False), flush=True)

# within each condition (>=4 paired patients), Macrophage vs top non-immune competitor (Stromal_Fibroblast, Epithelial)
cond_pairwise_rows = []
for cond, sub_cond in pbc.groupby("condition", observed=True):
    wide = sub_cond.pivot(index="patient_code", columns="cell_type", values="log2cpm_PTGS2")
    for ct in ["Stromal_Fibroblast", "Epithelial", "Endothelial"]:
        if REF_CT not in wide.columns or ct not in wide.columns:
            continue
        paired = wide[[REF_CT, ct]].dropna()
        if len(paired) < 4 or (paired[REF_CT] == paired[ct]).all():
            continue
        try:
            stat, p = wilcoxon(paired[REF_CT], paired[ct])
        except ValueError:
            continue
        cond_pairwise_rows.append({
            "condition": cond, "ref": REF_CT, "vs": ct, "n_paired_patients": len(paired),
            "median_diff_log2cpm_PTGS2": (paired[REF_CT] - paired[ct]).median(),
            "wilcoxon_stat": stat, "p_raw": p,
        })
cond_pairwise_df = pd.DataFrame(cond_pairwise_rows)
if len(cond_pairwise_df):
    cond_pairwise_df["p_bh"] = bh_correct(cond_pairwise_df["p_raw"].values)
cond_pairwise_df.to_csv(f"{OUT_DIR}/pairwise_wilcoxon_by_condition.csv", index=False)
print("\n=== Condition-stratified paired Wilcoxon: Macrophage vs Stromal/Epithelial/Endothelial (PTGS2) ===", flush=True)
print(cond_pairwise_df.round(4).to_string(index=False) if len(cond_pairwise_df) else "(insufficient paired N in any stratum)", flush=True)

print("\nDONE. peak GB:", mem(), flush=True)
