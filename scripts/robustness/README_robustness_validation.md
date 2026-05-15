# 内部稳健性验证运行说明

本目录脚本只复用当前发现队列和已有结果文件，不重跑 pySCENIC 主流程，不重跑 CellOracle 主流程。

## 输入

- `<PROJECT_ROOT_WINDOWS>\final_exports\astrocyte_pilot_rna_with_pyscenic_auc.h5ad`
- `<PROJECT_ROOT_WINDOWS>\analysis_outputs\celloracle_candidates\celloracle_candidate_tf_metrics.csv`
- `<PROJECT_ROOT_WINDOWS>\output\regulons.csv`

## 输出

所有结果统一写入：

- `<PROJECT_ROOT_WINDOWS>\robustness_validation\`

## 直接全跑

```powershell
cd <repository-root>
python scripts\robustness\run_all_robustness_validation.py
```

## 单步运行

```powershell
python scripts\robustness\01_inspect_input_object.py
python scripts\robustness\02_leave_one_donor_out.py
python scripts\robustness\03_pseudobulk_validation.py
python scripts\robustness\04_shortlist_sensitivity.py
python scripts\robustness\05_make_master_summary.py
```

## 依赖

Windows 本地 Python 需要：

```powershell
pip install pandas numpy scipy h5py matplotlib
```

这些脚本使用 `h5py + scipy.sparse` 直接读取 h5ad，不依赖 `anndata`，以避免本机 `torch/anndata` 导入问题。

