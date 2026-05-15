# GO/KEGG + 转录程序汇聚分析运行说明

本目录脚本只复用当前发现队列、pySCENIC regulon、候选 TF 与内部稳健性验证结果。

不会重跑 pySCENIC 主流程，不会重跑 CellOracle 主流程，不会执行药物重定位。

## 输出目录

`<PROJECT_ROOT_WINDOWS>\functional_interpretation\`

## 一键运行

```powershell
cd <repository-root>
python scripts\functional\run_all_functional_interpretation.py
```

## 单步运行

```powershell
python scripts\functional\00_revise_robustness_terms.py
python scripts\functional\01_check_resources.py
python scripts\functional\02_build_tf_gene_sets.py
python scripts\functional\03_run_go_kegg.py
python scripts\functional\04_program_convergence.py
python scripts\functional\05_make_functional_master_summary.py
```

## 依赖

Windows 本地 Python 需要：

```powershell
pip install pandas numpy scipy h5py matplotlib
```

`03_run_go_kegg.py` 不依赖 `gseapy/goatools`。它会优先从 Enrichr 下载并缓存：

- `GO_Biological_Process_2023`
- `KEGG_2021_Human`

若联网失败，脚本会输出清楚说明；前面的 gene set 构建与后面的 program convergence 框架仍可保留。

## 术语修正

当前 discovery 队列中 `donor_id` 只有 2 个 donor，`sample_id` 有 4 个 sample。因此上一轮稳健性验证应统一称为：

- `leave-one-sample-out / 逐样本剔除稳健性分析`
- `sample-level robustness`

不建议写成严格的 `donor-level leave-one-out robustness`。

