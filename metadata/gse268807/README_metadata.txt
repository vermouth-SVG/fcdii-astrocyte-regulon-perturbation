1. dataset_registry.csv：含 Phase 2 字段（platform、modality_type、n_samples、grouping_summary、has_control、raw_file_types、suited_for_analysis）。modality_type：scRNA/snRNA/multiome/mixed/unclear；has_control：yes/no/partial/unclear；suited_for_analysis 在公开稿中应使用 discovery/supportive/contextual 等保守口径。未定项可写「待核对」「待补充」「待GEO核对」等占位；校验脚本会提示仍待固化的列。
2. sample_metadata.csv：当前仅表头；逐样本信息须从 GEO/论文补充，勿编造。列含义见 列说明_sample_metadata.txt。
3. 校验：在项目根目录运行 python scripts/01_validate_metadata.py
4. 根目录的 *_template.csv 为早期模板，可与本目录正式表并存作参考。
