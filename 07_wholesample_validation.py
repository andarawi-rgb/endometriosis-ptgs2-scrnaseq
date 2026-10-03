"""
Cross-check: whole-tissue (all cell types pooled) patient-level pseudobulk PTGS2,
Control vs Ectopic_Peritoneal, to compare direction/magnitude against the
independent Lin et al. 2026 (BMC Women's Health) reanalysis of this same GSE179640
dataset, which reported PTGS2 downregulated at the whole-tissue pseudobulk level
in EcP vs Control (log2FC = -2.32, p = 0.045, no cell-type resolution).
"""
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, wilcoxon

OUT_DIR = "/tmp/claude-0/-home-claude/e7962fd6-be13-56c6-9438-14c50d5a6310/scratchpad/endometriosis"
GENES = ["PTGS2", "PTGS1", "PTGES", "PTGER2", "PTGER4"]

obs = pd.read_csv(f"{OUT_DIR}/obs_with_gene_counts.csv", index_col=0)

# whole-tissue pseudobulk per patient x condition (all cell types pooled, matching Lin et al.'s approach)
g = obs.groupby(["patient_code", "condition"], observed=True)
pb = g.agg(n_cells=("total_counts", "size"), libsize=("total_counts", "sum"),
           **{f"sum_{x}": (f"count_{x}", "sum") for x in GENES}).reset_index()
for x in GENES:
    cpm = pb[f"sum_{x}"] / pb["libsize"] * 1e6
    pb[f"log2cpm_{x}"] = np.log2(cpm + 1)
pb.to_csv(f"{OUT_DIR}/pseudobulk_wholetissue_patient_condition.csv", index=False)

print("=== Whole-tissue (all cell types pooled) pseudobulk PTGS2 log2(CPM+1), by condition ===")
print(pb.groupby("condition", observed=True)["log2cpm_PTGS2"].agg(["count", "median", "mean", "std"]).round(3).to_string())

ctrl = pb.loc[pb["condition"] == "Control", "log2cpm_PTGS2"]
ecp = pb.loc[pb["condition"] == "Ectopic_Peritoneal", "log2cpm_PTGS2"]
print(f"\nControl (n={len(ctrl)}): median={ctrl.median():.3f}")
print(f"Ectopic_Peritoneal (n={len(ecp)}): median={ecp.median():.3f}")
if len(ctrl) >= 2 and len(ecp) >= 2:
    stat, p = mannwhitneyu(ecp, ctrl, alternative="two-sided")
    log2fc = ecp.median() - ctrl.median()
    print(f"Mann-Whitney U (EcP vs Control, unpaired, different patients): U={stat:.1f}, p={p:.4f}")
    print(f"median log2FC (EcP - Control) = {log2fc:.3f}")
    print("Lin et al. 2026 (BMC Women's Health) reported: log2FC = -2.32, p = 0.045 (EcP vs Control, whole-tissue pseudobulk)")

# paired EuE vs EcP within same patient (many patients have both)
wide = pb.pivot(index="patient_code", columns="condition", values="log2cpm_PTGS2")
if "Eutopic_Endometrium" in wide.columns and "Ectopic_Peritoneal" in wide.columns:
    paired = wide[["Eutopic_Endometrium", "Ectopic_Peritoneal"]].dropna()
    print(f"\nPaired EuE vs EcP within same patient, n={len(paired)}:")
    print(paired.round(3).to_string())
    if len(paired) >= 4 and not (paired["Eutopic_Endometrium"] == paired["Ectopic_Peritoneal"]).all():
        stat, p = wilcoxon(paired["Eutopic_Endometrium"], paired["Ectopic_Peritoneal"])
        diff = (paired["Eutopic_Endometrium"] - paired["Ectopic_Peritoneal"]).median()
        print(f"Wilcoxon signed-rank: stat={stat:.1f}, p={p:.4f}, median(EuE - EcP)={diff:.3f}")
