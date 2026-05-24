# drug repositioning v2/v2.1 运行说明

本目录是当前合并后的药物重定位脚本目录。active 脚本以 v2/v2.1 为准；v1 脚本已归档到 `archive/legacy_v1_scripts/`。

本流程只重做药物重定位部分，不重跑任何上游分析。

不会重跑：

- pySCENIC
- CellOracle
- CellOracle perturbation
- 外部验证
- 内部稳健性验证
- GO/KEGG
- 转录程序汇聚

## 主线约束

- 主轴: `NFE2L2`, `THRB`
- 补充轴: `BHLHE40`
- 不进入主药物查询: `SOX2`

## 一键运行到 BBB 手工复核前

```powershell
cd <repository-root>
python scripts\drug_repositioning\run_drug_repositioning_until_bbb_v2.py
```

输出目录：

```text
<PROJECT_ROOT_WINDOWS>\drug_repositioning
```

## SwissADME 手工复核

复制这个文件到 SwissADME：

```text
<PROJECT_ROOT_WINDOWS>\drug_repositioning\06_swissadme_copy_paste_input.txt
```

药名对照表：

```text
<PROJECT_ROOT_WINDOWS>\drug_repositioning\06_swissadme_copy_paste_with_names.tsv
```

唯一手工回填文件：

```text
<PROJECT_ROOT_WINDOWS>\drug_repositioning\06_bbb_annotation_manual_review.csv
```

`manual_final_bbb_layer` 只允许：

- `CNS-directed candidates`
- `possible CNS-directed candidates`
- `peripheral/program-modulating candidates`

## 手工 BBB 后继续

```powershell
cd <repository-root>
python scripts\drug_repositioning\run_continue_after_manual_bbb_v2.py
```

## 当前最新 v2.1 after-BBB 结果

当前 canonical 结果已经位于：

```text
<PROJECT_ROOT_WINDOWS>\drug_repositioning
```

最重要的最新文件：

```text
13_integrated_after_bbb_v21.csv
14_primary_leads_after_bbb_v21.csv
15_supportive_after_bbb_v21.csv
16_paper_ready_table_after_bbb_v21.csv
17_after_bbb_summary_cn.txt
18_after_bbb_transition_note_cn.txt
```

## docking

docking 不是当前 v2 版本必需步骤，暂不执行。

## 解释边界

本流程输出的是 exploratory drug-signature enrichment clues / hypothesis-generating mechanism-direction clues，不是治疗候选、治疗建议或经过验证的治疗方案。部分历史文件名仍保留 `lead` 字样，仅作为内部排序标签，不应作为治疗性表述。
