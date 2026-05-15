import os
import numpy as np
import pandas as pd
import anndata as ad

H5AD_IN = "/in/astrocyte_pilot_rna_merged_with_raw.h5ad"
AUC_CSV = "/work/output/auc_mtx.csv"
OUT_DIR = "/work/final_exports"

OUT_H5AD = os.path.join(OUT_DIR, "astrocyte_pilot_rna_with_pyscenic_auc.h5ad")
OUT_AUC_MATCHED = os.path.join(OUT_DIR, "auc_mtx_matched_to_h5ad.csv")
OUT_REGULON_TXT = os.path.join(OUT_DIR, "regulon_names.txt")
OUT_AUC_MEAN = os.path.join(OUT_DIR, "regulon_auc_mean.csv")
OUT_REPORT = os.path.join(OUT_DIR, "merge_report.txt")

os.makedirs(OUT_DIR, exist_ok=True)

print("[1/5] Reading h5ad ...")
adata = ad.read_h5ad(H5AD_IN)
adata.obs_names = adata.obs_names.astype(str)

print("[2/5] Reading AUC matrix ...")
auc = pd.read_csv(AUC_CSV)
auc = auc.rename(columns={auc.columns[0]: "Cell"})
auc["Cell"] = auc["Cell"].astype(str)
auc = auc.drop_duplicates(subset="Cell").set_index("Cell")

print(f"h5ad cells: {adata.n_obs}")
print(f"h5ad genes: {adata.n_vars}")
print(f"AUC cells: {auc.shape[0]}")
print(f"AUC regulons: {auc.shape[1]}")

shared = adata.obs_names.intersection(auc.index)
missing_in_auc = [x for x in adata.obs_names if x not in auc.index]
missing_in_h5ad = [x for x in auc.index if x not in set(adata.obs_names)]

print(f"shared cells: {len(shared)}")
print(f"missing in AUC: {len(missing_in_auc)}")
print(f"missing in h5ad: {len(missing_in_h5ad)}")

if len(shared) == 0:
    raise ValueError("No shared cells between h5ad and AUC matrix.")

print("[3/5] Reordering and merging ...")
adata2 = adata[shared].copy()
auc2 = auc.loc[shared].copy()

adata2.obsm["X_pyscenic_auc"] = auc2.to_numpy(dtype=np.float32)
adata2.uns["pyscenic_regulon_names"] = auc2.columns.tolist()
adata2.uns["pyscenic_info"] = {
    "auc_csv": AUC_CSV,
    "n_cells": int(auc2.shape[0]),
    "n_regulons": int(auc2.shape[1]),
    "obsm_key": "X_pyscenic_auc",
}

print("[4/5] Writing outputs ...")
adata2.write_h5ad(OUT_H5AD)
auc2.to_csv(OUT_AUC_MATCHED)

with open(OUT_REGULON_TXT, "w", encoding="utf-8") as f:
    for x in auc2.columns:
        f.write(str(x) + "\n")

auc2.mean(axis=0).sort_values(ascending=False).to_csv(OUT_AUC_MEAN, header=["mean_auc"])

with open(OUT_REPORT, "w", encoding="utf-8") as f:
    f.write(f"h5ad_in={H5AD_IN}\n")
    f.write(f"auc_csv={AUC_CSV}\n")
    f.write(f"h5ad_cells={adata.n_obs}\n")
    f.write(f"h5ad_genes={adata.n_vars}\n")
    f.write(f"auc_cells={auc.shape[0]}\n")
    f.write(f"auc_regulons={auc.shape[1]}\n")
    f.write(f"shared_cells={len(shared)}\n")
    f.write(f"missing_in_auc={len(missing_in_auc)}\n")
    f.write(f"missing_in_h5ad={len(missing_in_h5ad)}\n")
    f.write(f"out_h5ad={OUT_H5AD}\n")
    f.write(f"out_auc_matched={OUT_AUC_MATCHED}\n")
    f.write(f"out_regulon_txt={OUT_REGULON_TXT}\n")
    f.write(f"out_auc_mean={OUT_AUC_MEAN}\n")

print("[5/5] Done.")
print("OUT_H5AD =", OUT_H5AD)
print("OUT_AUC_MATCHED =", OUT_AUC_MATCHED)
print("OUT_REGULON_TXT =", OUT_REGULON_TXT)
print("OUT_AUC_MEAN =", OUT_AUC_MEAN)
print("OUT_REPORT =", OUT_REPORT)
