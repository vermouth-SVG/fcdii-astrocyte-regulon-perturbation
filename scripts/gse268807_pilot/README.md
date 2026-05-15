# 脚本索引（与设计书阶段对应）

| 脚本 | 作用 |
|------|------|
| `00_check_environment.py` | Python/依赖检查 |
| `01_validate_metadata.py` | 校验 metadata CSV |
| `02_qc_and_downsample_scanpy.py` | 单 h5ad QC + 可选降采样 |
| `03_export_for_cloud.py` | 复制 h5ad + 导出 obs → `export_for_cloud/<GSE>/` |
| `04_ingest_cloud_grn_results.py` | 自 `cloud_inbox` 复制 CSV 到 `integrated/` |
| `05_composition_and_pseudobulk_prep.py` | donor×cell_type 组成长表 |
| `06_recovery_index_postprocess.py` | Recovery Index 排名（设计书 §10） |
| `07_enrichment_enrichr.py` | Enrichr（需 gseapy + 网络） |
| `08_download_geo_gse268807.py` | 自 GEO **HTTPS** 下载 GSE268807 至 `02_raw/GSE268807/`（curl 断点续传） |
| `09_mtx_gz_to_h5ad.py` | 单样本 MTX 三联 → h5ad（可按 feature_type 导出 Gene Expression 或 Peaks） |
| `10_concat_h5ad.py` | 合并多个 QC h5ad，生成本地子集或云端输入 |

**云端**：SCENIC+、pySCENIC、CellOracle — 见 `../00_docs/云端重计算Runbook.md`。

示例（有矩阵后）：

```bash
python scripts/02_qc_and_downsample_scanpy.py --input path/in.h5ad --output 03_processed/GSE268807/qc.h5ad --gse-id GSE268807 --max-cells 25000
python scripts/03_export_for_cloud.py --input 03_processed/GSE268807/qc.h5ad --gse-id GSE268807
```
