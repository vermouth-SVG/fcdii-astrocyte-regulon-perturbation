# pySCENIC/AUCell reproducibility audit

## 1. Executive summary

This report is an updated read-only audit of the project directory plus the WSL Linux-side files under `<WSL_HOME>`. No pySCENIC, ctx, AUCell, GRN, CellOracle, figure, table, or manuscript DOCX analyses were rerun or modified. Only this audit report and the companion CSV were updated.

The WSL search recovered the original discovery-run Docker commands from `<WSL_HOME>/.bash_history`. The discovery pySCENIC run used Docker image `aertslab/pyscenic:0.12.1`, with `--num_workers 4` for GRN, ctx, and AUCell. The exact commands generated `output/grn_adj.tsv`, `output/regulons.csv`, and `output/auc_mtx.csv` under `<PROJECT_ROOT_WSL>`, which is the WSL path corresponding to this Windows project directory.

Additional Docker Desktop/WSL checks support the user's recollection that Docker Desktop was opened and used as the Docker backend. The recovered commands were issued from the Ubuntu-D WSL environment, using `<PROJECT_ROOT_WSL>` bind mounts. Docker Desktop settings record `IntegratedWslDistros: ["Ubuntu-D"]`, the Docker CLI context `desktop-linux` points to `npipe:////./pipe/dockerDesktopLinuxEngine`, and Docker Desktop logs from 2026-04-16 show WSL integration for Ubuntu-D, including Docker context synchronization, mounting of `/var/run/docker.sock`, creation of `/usr/bin/docker`, and startup of the user-distro proxy. These records support Docker Desktop with WSL integration as the runtime environment. No retained per-command Docker Desktop daemon event containing `aertslab/pyscenic` or `<PROJECT_ROOT_WSL>` was found, so the manuscript should describe the analysis as run in Docker using `aertslab/pyscenic:0.12.1`, while the audit can document Docker Desktop/WSL integration as environment provenance.

The recovered GRN command did not explicitly include `--method`; the pySCENIC 0.12.1 CLI source retained in WSL shows that `pyscenic grn --method` defaults to `grnboost2`. Therefore the GRN method can be reported as GRNBoost2 by pySCENIC default, while noting that it was not explicitly written in the original command. The recovered GRN and AUCell commands did not include `--seed`. The pySCENIC 0.12.1 CLI source shows `default=None` for both GRN and AUCell seed arguments and states that the default is to use a random seed. Thus, no fixed random seed was used in the retained discovery commands. WSL contained cloud-run templates with `--seed 777`, but those were templates and do not match the actual local discovery-run Docker commands in `.bash_history`.

The merged discovery object contains 2,322 astrocytes, 36,601 genes, and 105 pySCENIC/AUCell regulons. The AUCell matrix was matched back to the h5ad object without cell loss (`shared_cells=2322`, `missing_in_auc=0`, `missing_in_h5ad=0`). The key pySCENIC resources are `db/allTFs_hg38.txt`, `db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`, and `db/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`.

Recommendation for the current manuscript sentence: delete the internal placeholder sentence "正式英文投稿前，应在本段进一步补入所用 TF list、ranking database、motif annotation 和随机种子等复现信息。". Replace it with the conservative methods text below. The only remaining caution is that the actual discovery run did not set a fixed seed; this should be stated as no fixed seed rather than inventing a seed.

## 2. Confirmed reproducibility items

| Item | Confirmed value | Evidence |
|---|---|---|
| Current mainline pySCENIC output directory | `output/` contains `grn_adj.tsv`, `regulons.csv`, and `auc_mtx.csv` | `RUN_ORDER_CURRENT_MAINLINE.md`, lines 18-26; `README_project_root.md`, line 35 |
| Original discovery Docker image | `aertslab/pyscenic:0.12.1` | `<WSL_HOME>/.bash_history`, lines 822-823 and 1099-1110 |
| Command launch environment | Ubuntu-D WSL, using `<PROJECT_ROOT_WSL>` bind mounts corresponding to this Windows project directory | `<WSL_HOME>/.bash_history`, lines 1099-1110 |
| Docker Desktop/WSL integration | Docker Desktop was configured to integrate the Ubuntu-D WSL distro; `desktop-linux` context targets Docker Desktop's Linux engine endpoint; Docker Desktop host logs show WSL integration setup for Ubuntu-D on 2026-04-16 | `<WINDOWS_USER_HOME>\AppData\Roaming\Docker\settings-store.json`, lines 3 and 6-8; `<WINDOWS_USER_HOME>\.docker\contexts\meta\fe9c6bd7a66301f49ca9b6a70b217107cd1284598bfc254700c989b916da791e\meta.json`, line 1; `<WINDOWS_USER_HOME>\AppData\Local\Docker\log\host\com.docker.backend.exe.log.4`, lines 2426, 2573, 2581, and 2605 |
| Docker Desktop version evidence | Docker Desktop 4.69.0 (build 224084) was installed/existing on 2026-04-16 | `<WINDOWS_USER_HOME>\AppData\Local\Docker\install-log.0.txt`, lines 2 and 14 |
| pySCENIC version | 0.12.1, via Docker image tag `aertslab/pyscenic:0.12.1` | `<WSL_HOME>/.bash_history`, lines 822-823 and 1099-1110 |
| Matrix conversion command | `docker run --rm -v <PROJECT_ROOT_WSL>:/work --entrypoint python aertslab/pyscenic:0.12.1 /work/scripts/mtx_to_csv.py ... /work/input/expression_for_pyscenic.csv` | `<WSL_HOME>/.bash_history`, line 1099 |
| Original GRN command | `pyscenic grn --num_workers 4 -o /data/output/grn_adj.tsv /data/input/expression_for_pyscenic.csv /data/db/allTFs_hg38.txt` inside `aertslab/pyscenic:0.12.1` | `<WSL_HOME>/.bash_history`, line 1102 |
| GRN method | GRNBoost2 by pySCENIC 0.12.1 default; `--method` was not explicitly supplied | pySCENIC CLI source `<WSL_HOME>/miniconda3/envs/scenicplus/lib/python3.11/site-packages/pyscenic/cli/pyscenic.py`, lines 516-521 |
| Original ctx command | `pyscenic ctx /data/output/grn_adj.tsv /data/db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather --annotations_fname /data/db/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname /data/input/expression_for_pyscenic.csv --mode custom_multiprocessing --output /data/output/regulons.csv --num_workers 4` | `<WSL_HOME>/.bash_history`, line 1106 |
| Original AUCell command | `pyscenic aucell /data/input/expression_for_pyscenic.csv /data/output/regulons.csv -o /data/output/auc_mtx.csv --num_workers 4` inside `aertslab/pyscenic:0.12.1` | `<WSL_HOME>/.bash_history`, line 1110 |
| Random seed | No fixed seed was supplied in the recovered discovery GRN or AUCell commands; pySCENIC 0.12.1 defaults to `seed=None`, using a random seed | `<WSL_HOME>/.bash_history`, lines 1102 and 1110; pySCENIC CLI source lines 523-528 and 691-696 |
| Input expression matrix | `input/expression_for_pyscenic.csv`; CSV cell-by-gene format; first column `CellID`; 2,323 lines including header; 36,602 header columns including `CellID`, corresponding to 2,322 cells x 36,601 genes | File dimension audit; first line begins `CellID,MIR1302-2HG,FAM138A,...`; first data row begins `GSE268807_G120_D_FL::AAAGCAAGTTCCGCAC-1,0,0,...` |
| Input h5ad object dimensions | 2,322 cells and 36,601 genes in the object used for AUCell merge | `final_exports/merge_report.txt`, lines 1, 3-4 |
| TF list file | `db/allTFs_hg38.txt`; 1,892 entries; examples include `ZNF354C`, `KLF12`, `ZNF143` | File dimension audit; `db/allTFs_hg38.txt`, first lines |
| TF list source | AertsLab cistarget TF list URL `https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt` | `<WSL_HOME>/ko_epilepsy/.../online_resource_source_candidates.csv`, row `official_hg38_tf_list`; `online_resource_download_log.txt`, line 5 |
| Ranking database | `db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` | File exists in `db/`; original ctx command uses this file; `logs/02_ctx.log`, lines 21-35 show this database loaded in workers |
| Ranking database source | AertsLab gene-based hg38/refseq_r80/mc_v10_clust 10 kb ranking database URL | `<WSL_HOME>/ko_epilepsy/.../online_resource_source_candidates.csv`, row `hg38_gene_ranking_db_10kb`; `online_resource_download_log.txt`, line 6 |
| Motif annotation file | `db/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`; 253,097 lines including header; 13 tab-delimited columns | File dimension audit; header begins `#motif_id motif_name motif_description source_name ... gene_name ...` |
| Motif annotation provenance | Reused from `<WSL_HOME>/ko_epilepsy/resources/scenicplus/motif_annotations/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`; direct URL was not recovered | `<WSL_HOME>/ko_epilepsy/.../online_resource_download_log.txt`, line 4; `.bash_history`, line 348 shows `wget -O ... motifs-v10nr...tbl \` but the URL line is not retained |
| GRN output | `output/grn_adj.tsv`; 6,449,819 lines including header; columns `TF`, `target`, `importance` | File dimension audit; `output/grn_adj.tsv`, first line |
| ctx dropout masking behavior | Dropout masking set to false | `logs/02_ctx.log`, line 8 states default behavior and line 9 states `Dropout masking is currently set to [False].` |
| ctx workers | 4 workers | `<WSL_HOME>/.bash_history`, line 1106; `logs/02_ctx.log`, lines 17 and 19 |
| Regulon output file | `output/regulons.csv`; 257 lines including multi-row header; columns represent `TF`, `MotifID`, and enrichment fields (`AUC`, `NES`, `MotifSimilarityQvalue`, `OrthologousIdentity`, `Annotation`, `Context`, `TargetGenes`, `RankAtMax`) | File dimension audit; `output/regulons.csv`, first 4 lines |
| AUCell matrix | `output/auc_mtx.csv`; 2,323 lines including header; 106 columns including `Cell`, corresponding to 2,322 cells x 105 regulons | File dimension audit; `final_exports/merge_report.txt`, lines 5-6 |
| AUCell-matched matrix | `final_exports/auc_mtx_matched_to_h5ad.csv`; 2,323 lines including header; 106 columns including cell index, corresponding to 2,322 cells x 105 regulons | File dimension audit; `final_exports/merge_report.txt`, lines 7-11 |
| Final regulon count | 105 regulons | `final_exports/merge_report.txt`, line 6; `analysis_outputs/structure_check/structure_summary.json`, line 37 |
| Alignment to 2,322 astrocytes | 2,322 AUC cells, 2,322 h5ad cells, 2,322 shared cells, no missing cells | `final_exports/merge_report.txt`, lines 3, 5, 7-9 |
| h5ad storage keys | AUC stored in `obsm["X_pyscenic_auc"]`; regulon names stored in `uns["pyscenic_regulon_names"]`; info stored in `uns["pyscenic_info"]` | `scripts/merge_auc_to_h5ad.py`, lines 48-54; `analysis_outputs/structure_check/structure_summary.json`, lines 18, 21-22 |
| Discovery group sizes used downstream | `internal_control=1434`, `lesion=888` | `analysis_outputs/group_compare/selected_group_counts.csv`, lines 1-3 |

## 3. Missing or uncertain items

| Item | Status | Detail |
|---|---|---|
| Fixed random seed value | Not fixed in actual run | The recovered discovery commands did not include `--seed`; pySCENIC 0.12.1 defaults to `seed=None` and uses a random seed. Cloud templates containing `--seed 777` were found, but they do not match the actual local Docker commands that produced the retained outputs. |
| Per-command Docker Desktop daemon event | Uncertain | Docker Desktop/WSL integration is strongly supported by settings, contexts, and 2026-04-16 host logs, but no retained daemon event containing `aertslab/pyscenic` or `<PROJECT_ROOT_WSL>` was found that directly records the exact pySCENIC command argv as a container-start event. |
| Direct motif annotation download URL | Not found | WSL records confirm reuse of the local motif annotation file, but the direct URL for `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl` was not recovered. |
| Exact Docker digest | Not found | The image tag `aertslab/pyscenic:0.12.1` is confirmed, but no immutable image digest was found. |

## 4. Evidence table

| Item | Evidence file | Evidence excerpt or dimensions | Confidence |
|---|---|---|---|
| Original discovery commands | `<WSL_HOME>/.bash_history` | Lines 1099, 1102, 1106, 1110 contain matrix conversion, GRN, ctx, and AUCell Docker commands | high |
| Docker image/version | `<WSL_HOME>/.bash_history` | Lines 822-823 show `docker pull aertslab/pyscenic:0.12.1` and `pyscenic -h`; lines 1099-1110 use the same image | high |
| Docker Desktop WSL integration settings | `<WINDOWS_USER_HOME>\AppData\Roaming\Docker\settings-store.json` | Lines 3 and 6-8 show `CustomWslDistroDir: <WINDOWS_ABSOLUTE_PATH>` and `IntegratedWslDistros: ["Ubuntu-D"]` | high |
| Docker Desktop CLI context | `<WINDOWS_USER_HOME>\.docker\contexts\meta\fe9c6bd7a66301f49ca9b6a70b217107cd1284598bfc254700c989b916da791e\meta.json` | `Name: desktop-linux`; `Description: Docker Desktop`; Docker endpoint `npipe:////./pipe/dockerDesktopLinuxEngine` | high |
| Docker Desktop user-distro activity | `<WSL_HOME>/.docker/desktop/log/host/docker-desktop-user-distro.log`; `<WINDOWS_USER_HOME>\AppData\Local\Docker\log\host\com.docker.backend.exe.log.4` | WSL log lines 2-3 show Docker Desktop syncing WSL Docker contexts. Host log lines 2426-2431 add a WSL agent for Ubuntu-D; line 2573 mounts `/var/run/docker.sock`; line 2581 creates `/usr/bin/docker`; line 2605 records `docker-desktop-user-distro proxy has started`; line 2608 starts `docker serve` under Ubuntu-D. | high |
| Docker Desktop installation/version | `<WINDOWS_USER_HOME>\AppData\Local\Docker\install-log.0.txt`; `<WINDOWS_USER_HOME>\AppData\Roaming\Docker Desktop\versions.json` | Docker Desktop version 4.69.0 (build 224084); existing installation found on 2026-04-16 | high |
| GRN method default | pySCENIC 0.12.1 CLI source in WSL | Lines 516-521 show `--method` choices and `default="grnboost2"` | high |
| GRN seed default | pySCENIC 0.12.1 CLI source in WSL | Lines 523-528 show `--seed`, `default=None`, and default random seed behavior | high |
| AUCell seed default | pySCENIC 0.12.1 CLI source in WSL | Lines 691-696 show `--seed`, `default=None`, and default random seed behavior | high |
| Mainline output role | `RUN_ORDER_CURRENT_MAINLINE.md` | `Stores pySCENIC outputs: grn_adj.tsv, regulons.csv, and auc_mtx.csv.` | high |
| Input expression matrix dimensions | `input/expression_for_pyscenic.csv` | 2,323 lines including header; 36,602 header columns including `CellID`; first row is a GSE268807 cell barcode followed by integer counts | high |
| Input h5ad and AUC merge | `final_exports/merge_report.txt` | `h5ad_cells=2322`, `h5ad_genes=36601`, `auc_cells=2322`, `auc_regulons=105`, `shared_cells=2322`, `missing_in_auc=0`, `missing_in_h5ad=0` | high |
| Final object structure | `analysis_outputs/structure_check/structure_summary.json` | `obsm_keys: X_pyscenic_auc`; `uns_keys: pyscenic_info, pyscenic_regulon_names`; `n_regulons: 105` | high |
| TF list | `db/allTFs_hg38.txt` | 1,892 lines; first entries `ZNF354C`, `KLF12`, `ZNF143`; filename includes `hg38` | high |
| TF list source | WSL `online_resource_source_candidates.csv` and `online_resource_download_log.txt` | `https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt`; status OK/content_length 11690 | high |
| Ranking database source | WSL `online_resource_source_candidates.csv` and `online_resource_download_log.txt` | AertsLab URL for `hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`; content_length 311298530 | high |
| Motif annotation | `db/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl` | 253,097 lines including header; 13 tab-delimited columns; filename includes `hgnc` | high |
| GRN adjacency output | `output/grn_adj.tsv` | 6,449,819 lines including header; columns `TF`, `target`, `importance` | high |
| ctx log | `logs/02_ctx.log` | `Creating modules`; `Calculating Pearson correlations`; `Loading databases`; `Calculating regulons`; `Using 4 workers`; database and motif annotations loaded | high |
| ctx dropout masking | `logs/02_ctx.log` | `Dropout masking is currently set to [False].` | high |
| Regulon refinement output | `output/regulons.csv` | 257 lines including multi-row header; fields include `TF`, `MotifID`, `AUC`, `NES`, `MotifSimilarityQvalue`, `Annotation`, `Context`, `TargetGenes` | high |
| AUCell log | `logs/03_aucell.log` | `Loading expression matrix`; `Loading gene signatures`; `Calculating cellular enrichment`; `Writing results to file` | high |
| AUCell matrix | `output/auc_mtx.csv` | 2,323 lines including header; 106 columns including `Cell`; 105 regulon columns | high |
| Matched AUCell matrix | `final_exports/auc_mtx_matched_to_h5ad.csv` | 2,323 lines including header; 106 columns including cell index; written by `scripts/merge_auc_to_h5ad.py` | high |

## 5. Exact commands or scripts identified

The original discovery-run commands were recovered from WSL `<WSL_HOME>/.bash_history`. The commands were issued from the Ubuntu-D WSL environment and used `<PROJECT_ROOT_WSL>` bind mounts. Docker Desktop/WSL integration evidence indicates that Docker Desktop was the runtime backend for this WSL Docker workflow, although no per-command daemon event tying each command line to a container-start record was recovered.

Matrix conversion:

```bash
docker run --rm   -v <PROJECT_ROOT_WSL>:/work   --entrypoint python   aertslab/pyscenic:0.12.1   /work/scripts/mtx_to_csv.py   /work/input/discovery_counts_raw.mtx   /work/input/discovery_barcodes.tsv   /work/input/discovery_genes.tsv   /work/input/expression_for_pyscenic.csv
```

GRN inference:

```bash
time docker run --rm   -v <PROJECT_ROOT_WSL>:/data   aertslab/pyscenic:0.12.1 pyscenic grn   --num_workers 4   -o /data/output/grn_adj.tsv   /data/input/expression_for_pyscenic.csv   /data/db/allTFs_hg38.txt   2>&1 | tee <PROJECT_ROOT_WSL>/logs/01_grn.log
```

ctx/regulon refinement:

```bash
time docker run --rm   -v <PROJECT_ROOT_WSL>:/data   aertslab/pyscenic:0.12.1 pyscenic ctx   /data/output/grn_adj.tsv   /data/db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather   --annotations_fname /data/db/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl   --expression_mtx_fname /data/input/expression_for_pyscenic.csv   --mode custom_multiprocessing   --output /data/output/regulons.csv   --num_workers 4   2>&1 | tee <PROJECT_ROOT_WSL>/logs/02_ctx.log
```

AUCell scoring:

```bash
time docker run --rm   -v <PROJECT_ROOT_WSL>:/data   aertslab/pyscenic:0.12.1 pyscenic aucell   /data/input/expression_for_pyscenic.csv   /data/output/regulons.csv   -o /data/output/auc_mtx.csv   --num_workers 4   2>&1 | tee <PROJECT_ROOT_WSL>/logs/03_aucell.log
```

The script that merged the primary AUCell matrix back to the discovery h5ad object is `scripts/merge_auc_to_h5ad.py`, with the following key paths:

```python
H5AD_IN = "/in/astrocyte_pilot_rna_merged_with_raw.h5ad"
AUC_CSV = "/work/output/auc_mtx.csv"
OUT_H5AD = "/work/final_exports/astrocyte_pilot_rna_with_pyscenic_auc.h5ad"
OUT_AUC_MATCHED = "/work/final_exports/auc_mtx_matched_to_h5ad.csv"
OUT_REGULON_TXT = "/work/final_exports/regulon_names.txt"
OUT_AUC_MEAN = "/work/final_exports/regulon_auc_mean.csv"
```

## 6. Recommended manuscript text

建议删除原句：

> 正式英文投稿前，应在本段进一步补入所用 TF list、ranking database、motif annotation 和随机种子等复现信息。

建议替换为以下中文方法段落：

pySCENIC/AUCell regulon activity analysis：在 GSE268807 astrocyte pilot 分析对象中，将细胞×基因表达矩阵导出并转换为 `input/expression_for_pyscenic.csv`（CSV 格式，第一列为 `CellID`；2,322 个细胞，36,601 个基因）后进行 pySCENIC/AUCell 分析。pySCENIC 分析在 Docker 容器 `aertslab/pyscenic:0.12.1` 中完成。GRN 推断使用 `pyscenic grn`，输入为 `expression_for_pyscenic.csv` 和 TF list `db/allTFs_hg38.txt`（1,892 个 TF 条目，来源记录为 AertsLab cistarget `allTFs_hg38.txt`），输出为 `output/grn_adj.tsv`，并设置 `--num_workers 4`。原始 GRN 命令未显式指定 `--method`；根据 pySCENIC 0.12.1 CLI 源码，`pyscenic grn` 的默认方法为 GRNBoost2。ctx/motif regulon 精炼使用 `output/grn_adj.tsv`、cisTarget ranking database `db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` 和 motif annotation 文件 `db/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`，参数包括 `--expression_mtx_fname /data/input/expression_for_pyscenic.csv`、`--mode custom_multiprocessing` 和 `--num_workers 4`，输出为 `output/regulons.csv`。日志显示 ctx 步骤加载了上述 hg38 ranking database 和 motif annotation，且 dropout masking 为 false。AUCell 活性评分使用 `pyscenic aucell /data/input/expression_for_pyscenic.csv /data/output/regulons.csv -o /data/output/auc_mtx.csv --num_workers 4`。原始 GRN 和 AUCell 命令均未设置 `--seed`；pySCENIC 0.12.1 CLI 记录显示 seed 默认值为 `None`，即未固定随机种子。随后将 AUCell 矩阵合并回 h5ad 对象，保存为 `final_exports/astrocyte_pilot_rna_with_pyscenic_auc.h5ad`，其中 AUCell 活性存储于 `obsm["X_pyscenic_auc"]`，regulon 名称存储于 `uns["pyscenic_regulon_names"]`。合并记录显示 AUCell 矩阵与 h5ad 中的 2,322 个 astrocytes 完全对齐，未出现缺失或额外细胞；最终纳入后续比较的 pySCENIC regulon 为 105 个。本研究使用 pySCENIC/AUCell 的目的，是在单细胞 RNA 层面获得可比较的 TF-regulon activity，而不是构建 enhancer-driven 多组学调控网络。

## 7. English BMC-style manuscript paragraph

### pySCENIC/AUCell regulon activity analysis

Regulon activity was assessed from the GSE268807 astrocyte pilot object using a cell-by-gene expression matrix exported as `input/expression_for_pyscenic.csv` (CSV format with `CellID` as the first column; 2,322 cells and 36,601 genes). The analysis was run in Docker using `aertslab/pyscenic:0.12.1`. The expression matrix was generated from the retained Matrix Market, barcode, and gene files using the project script `mtx_to_csv.py`. GRN inference was performed with `pyscenic grn`, using `expression_for_pyscenic.csv` and the transcription factor list `allTFs_hg38.txt` (1,892 entries; recorded source: AertsLab cistarget TF list), with four workers and output written to `output/grn_adj.tsv`. The recovered GRN command did not explicitly set `--method`; in the pySCENIC 0.12.1 CLI source, the default GRN method is GRNBoost2. Regulon refinement was then performed with `pyscenic ctx` using `output/grn_adj.tsv`, the ranking database `hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`, the motif annotation table `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`, `--expression_mtx_fname /data/input/expression_for_pyscenic.csv`, `--mode custom_multiprocessing`, and four workers. The ctx log confirmed loading of the hg38 ranking database and motif annotations, and recorded dropout masking as false. AUCell scoring was performed with `pyscenic aucell` using `expression_for_pyscenic.csv` and `output/regulons.csv`, with four workers and output written to `output/auc_mtx.csv`. No fixed `--seed` argument was supplied in the recovered GRN or AUCell commands; the pySCENIC 0.12.1 CLI source records the default seed as `None`, corresponding to use of a random seed. The AUCell matrix contained 2,322 cells and 105 regulons and was merged back into the h5ad object as `final_exports/astrocyte_pilot_rna_with_pyscenic_auc.h5ad`, with activity values stored in `obsm["X_pyscenic_auc"]` and regulon names stored in `uns["pyscenic_regulon_names"]`. The merge report confirmed complete alignment between the AUCell matrix and the 2,322 astrocytes in the h5ad object, with no missing or extra cells. This analysis was used to obtain comparable RNA-level TF-regulon activity scores for downstream group comparisons, rather than to construct an enhancer-driven multi-omic gene regulatory network.

## 8. Items that require manual confirmation before submission

1. Decide whether to explicitly state in the manuscript that no fixed random seed was supplied for pySCENIC/AUCell, or to leave this detail in supplementary reproducibility notes.
2. If journal policy requires immutable containers, recover or record the Docker image digest for `aertslab/pyscenic:0.12.1`; only the image tag was found.
3. If possible, recover the direct source URL for `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`; current evidence confirms the local reused file and filename but not the direct URL.
4. If submission requires exact host/container runtime provenance, manually confirm whether to report Docker Desktop 4.69.0 with Ubuntu-D WSL integration. The current project evidence strongly supports this environment, but no per-command Docker Desktop daemon event was recovered.
