import scanpy as sc
import pandas as pd
import numpy as np
import gc, resource

def mem():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024

OUT_DIR = "/tmp/claude-0/-home-claude/e7962fd6-be13-56c6-9438-14c50d5a6310/scratchpad/endometriosis"

adata = sc.read_h5ad(f"{OUT_DIR}/adata_norm.h5ad")
print("Loaded:", adata.shape, "peak GB:", mem(), flush=True)

leiden = pd.read_csv(f"{OUT_DIR}/leiden_labels.csv", index_col=0)
harmony_pca = np.load(f"{OUT_DIR}/harmony_pca.npy")

assert list(leiden.index) == list(adata.obs_names), "obs_names order mismatch between leiden_labels.csv and adata_norm.h5ad"
assert harmony_pca.shape[0] == adata.n_obs

adata.obs["leiden"] = pd.Categorical(leiden["leiden"].astype(str).values,
                                      categories=[str(i) for i in sorted(leiden["leiden"].unique())])
adata.obsm["X_pca_harmony"] = harmony_pca
print("Merged leiden + harmony_pca onto full-gene object.", "n_clusters:", adata.obs.leiden.nunique(), flush=True)

# ---- canonical marker sets for annotation ----
MARKERS = {
    "Epithelial":         ["EPCAM", "PAX8", "KRT8", "KRT18", "KRT19"],
    "Stromal_Fibroblast":  ["DCN", "PDGFRA", "PDGFRB", "COL1A1", "COL3A1", "LUM"],
    "Endothelial":        ["PECAM1", "VWF", "CDH5", "CLDN5"],
    "Perivascular_Mural":  ["ACTA2", "RGS5", "MYH11", "TAGLN"],
    "Mesothelial":        ["MSLN", "WT1", "UPK3B"],
    "Macrophage_Myeloid":  ["PTPRC", "CD68", "CD14", "ITGAM", "LYZ", "AIF1", "C1QA", "C1QB"],
    "T_NK_Lymphocyte":     ["CD3E", "CD3D", "CD2", "NKG7", "GNLY", "KLRD1"],
    "B_Plasma":            ["CD79A", "MS4A1", "MZB1", "JCHAIN"],
    "Mast":                ["TPSAB1", "TPSB2", "CPA3", "KIT"],
    "Neutrophil_Granulocyte": ["S100A8", "S100A9", "FCGR3B", "CSF3R", "MNDA", "NAMPT"],
}
PROSTAGLANDIN_GENES = ["PTGS2", "PTGS1", "PTGES", "PTGER2", "PTGER4"]

present = {ct: [g for g in gl if g in adata.var_names] for ct, gl in MARKERS.items()}
missing = {ct: [g for g in gl if g not in adata.var_names] for ct, gl in MARKERS.items()}
print("Marker genes present per cell type:", present, flush=True)
print("Marker genes MISSING per cell type:", {k: v for k, v in missing.items() if v}, flush=True)
print("Prostaglandin-pathway genes present:",
      [g for g in PROSTAGLANDIN_GENES if g in adata.var_names],
      "MISSING:", [g for g in PROSTAGLANDIN_GENES if g not in adata.var_names], flush=True)

# mean expression (of log-normalized X) per marker gene per leiden cluster -> score = mean of z-scored marker means
import scipy.sparse as sp
X = adata.X
if sp.issparse(X):
    X = X.tocsc()

def gene_mean_per_cluster(gene):
    gi = adata.var_names.get_loc(gene)
    col = np.asarray(X[:, gi].todense()).ravel() if sp.issparse(X) else X[:, gi]
    return pd.Series(col, index=adata.obs_names).groupby(adata.obs["leiden"], observed=True).mean()

cluster_ids = sorted(adata.obs["leiden"].unique(), key=lambda x: int(x))
score_table = pd.DataFrame(index=cluster_ids)
for ct, genes in present.items():
    if not genes:
        continue
    gene_means = pd.DataFrame({g: gene_mean_per_cluster(g) for g in genes})
    gene_means = gene_means.reindex(cluster_ids)
    # z-score each gene across clusters, then average -> comparable "marker score"
    z = (gene_means - gene_means.mean()) / gene_means.std(ddof=0).replace(0, np.nan)
    score_table[ct] = z.mean(axis=1)

print("\n=== Cell-type marker score per Leiden cluster (z-scored mean expression) ===", flush=True)
print(score_table.round(2).to_string(), flush=True)

cluster_call = score_table.idxmax(axis=1)
print("\n=== Provisional cell-type call per cluster (argmax marker score) ===", flush=True)
print(cluster_call, flush=True)

adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_call).astype(str)
print("\nCell-type composition (n cells):", flush=True)
print(adata.obs["cell_type"].value_counts(), flush=True)

# also dump raw marker gene means per cluster for manual sanity-check / write-up
raw_marker_means = pd.DataFrame(index=cluster_ids)
for ct, genes in present.items():
    for g in genes:
        raw_marker_means[g] = gene_mean_per_cluster(g).reindex(cluster_ids)
raw_marker_means.to_csv(f"{OUT_DIR}/cluster_marker_means.csv")

score_table.to_csv(f"{OUT_DIR}/cluster_celltype_scores.csv")
cluster_call.to_csv(f"{OUT_DIR}/cluster_celltype_calls.csv")
adata.obs[["sample_id", "gsm", "patient_code", "condition_code", "condition", "leiden", "cell_type"]].to_csv(f"{OUT_DIR}/obs_annotated.csv")

adata.write(f"{OUT_DIR}/adata_annotated.h5ad")
print("\nSaved adata_annotated.h5ad, cluster_marker_means.csv, cluster_celltype_scores.csv, cluster_celltype_calls.csv, obs_annotated.csv. peak GB:", mem(), flush=True)
