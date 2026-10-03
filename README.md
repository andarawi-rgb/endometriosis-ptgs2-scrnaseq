# Single-cell PTGS2 (COX-2) cellular-source analysis in endometriosis

Analysis code accompanying the manuscript: *"Single-cell resolution identifies myeloid and epithelial cells, not classical stromal fibroblasts, as the dominant source of PTGS2 (COX-2) transcript in endometriosis"* (Andarawi MO, Elhag OO; manuscript in preparation, target journal BMC Women's Health).

## Data source

Raw, publicly deposited single-cell count matrices for 31 samples from 14 patients (Control, Eutopic Endometrium, Ectopic Peritoneal, Ectopic Peritoneal-Adjacent, Ectopic Ovary), GEO accession **[GSE179640](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179640)** (Tan et al., *Nature Cell Biology* 2022). This repository does not redistribute the raw data — download the 10x Genomics `filtered_feature_bc_matrix.h5` sample files from GEO before running the pipeline.

## Pipeline

Scripts are numbered in the order they were run.

| Script | Purpose |
|---|---|
| `01_load_qc.py`, `01a_load_group.py`, `01b_merge_pair.py` | Load raw per-sample count matrices, establish a common gene panel across all 31 files |
| `02_final_qc.py`, `02b_apply_qc_subset.py`, `02c_chunked_subset.py` | QC filtering (mitochondrial fraction, detected-gene count) and chunked subsetting for memory efficiency |
| `03_normalize_hvg.py` | Normalization, log1p transform, highly-variable-gene selection |
| `04_pca_harmony_cluster.py` | PCA, Harmony batch integration (by patient), Leiden clustering |
| `05_merge_annotate.py` | Marker-panel z-score cell-type annotation (10 cell types), cross-checked by unbiased top-marker-gene discovery |
| `06_pseudobulk_pge2.py` | Patient-level pseudobulk statistics for *PTGS2*/*PTGS1*/*PTGES*/*PTGER2*/*PTGER4* (paired Wilcoxon signed-rank, Benjamini-Hochberg correction, Kruskal-Wallis omnibus test) |
| `07_wholesample_validation.py` | Whole-tissue pseudobulk validation (Control vs. Ectopic Peritoneal; paired Eutopic vs. Ectopic comparison), cross-checked against an independent concurrent reanalysis of the same dataset (Lin et al., *BMC Women's Health* 2026) |
| `08_unbiased_markers.py` | Unbiased top-marker-gene discovery per cluster, independent of the annotation marker panel |
| `09_umap_dotplot.py` | UMAP and marker-panel dot-plot figure generation |

## Environment

```
pip install -r requirements.txt
```

## Citation

If you use this code, please cite the manuscript (citation to be added on publication) and the original data source (Tan et al., *Nature Cell Biology* 2022; GEO: GSE179640).

## License

MIT — see `LICENSE`.

## Use of AI tools

Analysis code was developed with the assistance of Claude (Anthropic), under the authors' direct supervision and independent verification of intermediate outputs and statistical results. See the manuscript's "Use of AI tools" declaration for full detail.
