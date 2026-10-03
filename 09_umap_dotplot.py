import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp
import gc, resource
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def mem():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024

OUT_DIR = "/tmp/claude-0/-home-claude/e7962fd6-be13-56c6-9438-14c50d5a6310/scratchpad/endometriosis"

# ---------- 1) UMAP from the (lightweight) Harmony PCA embedding ----------
harmony_pca = np.load(f"{OUT_DIR}/harmony_pca.npy")
obs_ann = pd.read_csv(f"{OUT_DIR}/obs_annotated.csv", index_col=0)
assert harmony_pca.shape[0] == len(obs_ann)

a = ad.AnnData(X=np.zeros((len(obs_ann), 1), dtype=np.float32), obs=obs_ann.copy())
a.obsm["X_pca_harmony"] = harmony_pca
print("Minimal AnnData built:", a.shape, "peak GB:", mem(), flush=True)

sc.pp.neighbors(a, use_rep="X_pca_harmony", n_neighbors=15)
print("neighbors done. peak GB:", mem(), flush=True)
sc.tl.umap(a)
print("UMAP done. peak GB:", mem(), flush=True)

a.obs[["cell_type"]].assign(UMAP1=a.obsm["X_umap"][:, 0], UMAP2=a.obsm["X_umap"][:, 1]).to_csv(f"{OUT_DIR}/umap_coords.csv")

fig, ax = plt.subplots(figsize=(8, 7))
cell_types = sorted(a.obs["cell_type"].unique())
cmap = plt.get_cmap("tab10")
for i, ct in enumerate(cell_types):
    m = a.obs["cell_type"] == ct
    ax.scatter(a.obsm["X_umap"][m.values, 0], a.obsm["X_umap"][m.values, 1],
               s=1.5, alpha=0.5, label=f"{ct} (n={m.sum()})", color=cmap(i % 10))
ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
ax.set_title("GSE179640 re-analysis: cell-type annotation (Harmony-corrected, n=109,100 cells)")
lgd = ax.legend(markerscale=8, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/fig_umap_celltype.png", dpi=200, bbox_inches="tight")
print("Saved fig_umap_celltype.png", flush=True)

del a
gc.collect()
print("peak GB after UMAP block:", mem(), flush=True)

# ---------- 2) Dot plot: canonical markers + prostaglandin pathway genes, per cell type ----------
MARKER_PANEL = {
    "Epithelial":         ["EPCAM", "PAX8", "KRT8", "KRT18", "KRT19"],
    "Stromal_Fibroblast":  ["DCN", "PDGFRA", "PDGFRB", "COL1A1", "COL3A1", "LUM"],
    "Endothelial":        ["PECAM1", "VWF", "CDH5", "CLDN5"],
    "Perivascular_Mural":  ["ACTA2", "RGS5", "MYH11", "TAGLN"],
    "Mesothelial":        ["MSLN", "WT1", "UPK3B"],
    "Macrophage_Myeloid":  ["PTPRC", "CD68", "CD14", "ITGAM", "LYZ", "AIF1", "C1QA", "C1QB"],
    "Neutrophil_Granulocyte": ["S100A8", "S100A9", "FCGR3B", "CSF3R", "MNDA"],
    "T_NK_Lymphocyte":     ["CD3E", "CD3D", "CD2", "NKG7", "GNLY", "KLRD1"],
    "B_Plasma":            ["CD79A", "MS4A1", "MZB1", "JCHAIN"],
    "Mast":                ["TPSAB1", "TPSB2", "CPA3", "KIT"],
}
PROSTAGLANDIN_GENES = ["PTGS2", "PTGS1", "PTGES", "PTGER2", "PTGER4"]
all_genes = [g for gl in MARKER_PANEL.values() for g in gl] + PROSTAGLANDIN_GENES

adata = sc.read_h5ad(f"{OUT_DIR}/adata_annotated.h5ad")
print("Loaded adata_annotated for dotplot:", adata.shape, "peak GB:", mem(), flush=True)
present_genes = [g for g in all_genes if g in adata.var_names]
adata_small = adata[:, present_genes].copy()
del adata
gc.collect()
print("Subset to marker genes:", adata_small.shape, "peak GB:", mem(), flush=True)

sc.pl.dotplot(adata_small, var_names=present_genes, groupby="cell_type",
              categories_order=["Epithelial", "Stromal_Fibroblast", "Endothelial", "Perivascular_Mural",
                                 "Mesothelial", "Macrophage_Myeloid", "Neutrophil_Granulocyte",
                                 "T_NK_Lymphocyte", "B_Plasma", "Mast"],
              standard_scale="var", show=False, save=None,
              title="Canonical cell-type markers + PTGS2/PTGS1/PTGES/PTGER2/PTGER4 by annotated cell type")
plt.savefig(f"{OUT_DIR}/fig_dotplot_markers.png", dpi=200, bbox_inches="tight")
print("Saved fig_dotplot_markers.png. peak GB:", mem(), flush=True)

print("DONE.", flush=True)
