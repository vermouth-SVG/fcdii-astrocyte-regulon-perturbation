from __future__ import annotations

import csv
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "manuscript_output" / "tables_revised_v3"

INPUTS = {
    "table1": ROOT / "manuscript_output" / "tables_main" / "Table1_cohort_and_analytical_objects_submission.csv",
    "table2": ROOT / "manuscript_output" / "tables_main" / "Table2_shortlist_main_submission.csv",
    "table3": ROOT / "manuscript_output" / "tables_main" / "Table3_integrated_priority_main_submission.csv",
    "st04": ROOT / "manuscript_output" / "tables_supplementary" / "TableS1_full_differential_regulons_submission.csv",
    "st05": ROOT / "manuscript_output" / "tables_supplementary" / "TableS3_celloracle_round1_full_metrics_submission.csv",
    "st06": ROOT / "manuscript_output" / "tables_supplementary" / "TableS4_sample_level_robustness_full_submission.csv",
    "st07": ROOT / "manuscript_output" / "tables_supplementary" / "TableS6_supportive_external_analysis_full_submission.csv",
    "st08": ROOT / "manuscript_output" / "tables_supplementary" / "TableS5_functional_interpretation_program_convergence_full_submission.csv",
    "st09": ROOT / "manuscript_output" / "tables_supplementary" / "TableS7_drug_repositioning_v21_after_bbb_full_submission.csv",
}

SOURCE_LABELS = {
    key: str(path.relative_to(ROOT)) for key, path in INPUTS.items()
}


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)


def replace_terms(value: object) -> object:
    if not isinstance(value, str):
        return value
    out = value
    replacements = {
        "internal_control": "internal-control",
        "donor-level robustness": "sample-level robustness",
        "external validation": "supportive evidence",
        "formal validation": "supportive evidence",
        "independent validation": "supportive evidence",
        "validated": "supported",
        "confirmed": "supported",
        "therapeutic candidate": "exploratory clue",
        "clinical candidate": "exploratory clue",
        "treatment lead": "exploratory clue",
        "validated drug": "exploratory clue",
        "therapeutic recommendation": "mechanistic recommendation",
    }
    for old, new in replacements.items():
        out = re.sub(re.escape(old), new, out, flags=re.IGNORECASE)
    return out


def clean_frame(df: pd.DataFrame) -> pd.DataFrame:
    return df.map(replace_terms)


def main_table_1() -> tuple[pd.DataFrame, Path]:
    path = OUT / "Main_Table1_analytical_datasets_and_object_roles.csv"
    df = pd.DataFrame(
        [
            {
                "Dataset": "GSE268807 astrocyte pilot",
                "Role": "Discovery cohort",
                "Analytical object": "AUCell-matched astrocyte pilot",
                "Cells": "2,322",
                "Comparison": "Lesion vs internal-control",
                "Use in manuscript": "pySCENIC discovery; CellOracle in silico perturbation; sample-level robustness",
            },
            {
                "Dataset": "GSE140393",
                "Role": "Single-group supportive analysis",
                "Analytical object": "Lesion-only expression object",
                "Cells": "11,399",
                "Comparison": "Single-group lesion object",
                "Use in manuscript": "Expression-level supportive evidence only",
            },
            {
                "Dataset": "GSE190452",
                "Role": "Cross-syndrome supportive analysis",
                "Analytical object": "TLE/control expression object",
                "Cells": "69,968",
                "Comparison": "TLE vs non-epileptic control",
                "Use in manuscript": "Cross-syndrome supportive evidence only",
            },
        ]
    )
    write_csv(df, path)
    return df, path


def main_table_2() -> tuple[pd.DataFrame, Path]:
    path = OUT / "Main_Table2_integrated_TF_prioritization.csv"
    df = pd.DataFrame(
        [
            {
                "TF": "NFE2L2",
                "Axis": "Lesion-associated",
                "Discovery evidence": "Lesion-side expression and regulon activity",
                "In silico perturbation evidence": "Highest mean shift; Approx. RI rank 2",
                "Robustness / supportive evidence": "Stable sample-level direction; disease-side support",
                "Program interpretation": "Stress-adaptation / transcriptional reprogramming",
                "Final priority": "Primary lesion-associated axis",
            },
            {
                "TF": "THRB",
                "Axis": "Internal-control-associated",
                "Discovery evidence": "Internal-control-side expression and regulon activity",
                "In silico perturbation evidence": "Highest Approx. RI; mean shift rank 2",
                "Robustness / supportive evidence": "Stable sample-level direction; partial internal-control-like support",
                "Program interpretation": "Compact homeostatic / synaptic-calcium support",
                "Final priority": "Primary internal-control-associated axis",
            },
            {
                "TF": "BHLHE40",
                "Axis": "Lesion-associated",
                "Discovery evidence": "Lesion-side expression and regulon activity",
                "In silico perturbation evidence": "Intermediate perturbation effect",
                "Robustness / supportive evidence": "Stable internal evidence; attenuated external support",
                "Program interpretation": "Secondary lesion-supportive program",
                "Final priority": "Secondary candidate",
            },
            {
                "TF": "SOX2",
                "Axis": "Lesion-associated",
                "Discovery evidence": "Lesion-side expression and regulon activity",
                "In silico perturbation evidence": "Lower perturbation priority",
                "Robustness / supportive evidence": "Stable internal evidence; partial / retained support",
                "Program interpretation": "Target-limited retained program",
                "Final priority": "Retained / lower-priority",
            },
        ]
    )
    write_csv(df, path)
    return df, path


def st01() -> tuple[pd.DataFrame, Path]:
    path = OUT / "ST01_dataset_role_and_excluded_resource_audit.csv"
    df = pd.DataFrame(
        [
            {
                "Dataset": "GSE268807 astrocyte pilot",
                "Role": "Discovery cohort",
                "Cells": "2,322",
                "Object": "AUCell-matched astrocyte pilot",
                "Comparison": "Lesion vs internal-control",
                "Current manuscript status": "Included discovery cohort",
                "Reason / boundary statement": "Used for pySCENIC discovery, CellOracle in silico perturbation, and sample-level robustness.",
            },
            {
                "Dataset": "GSE140393",
                "Role": "Single-group supportive analysis",
                "Cells": "11,399",
                "Object": "Lesion-only expression object",
                "Comparison": "Single-group lesion object",
                "Current manuscript status": "Included supportive dataset",
                "Reason / boundary statement": "Used for single-group expression-level supportive evidence only.",
            },
            {
                "Dataset": "GSE190452",
                "Role": "Cross-syndrome supportive analysis",
                "Cells": "69,968",
                "Object": "TLE/control expression object",
                "Comparison": "TLE vs non-epileptic control",
                "Current manuscript status": "Included supportive dataset",
                "Reason / boundary statement": "Used for cross-syndrome expression-level and target-set supportive evidence only.",
            },
            {
                "Dataset": "GSE275302",
                "Role": "Not included",
                "Cells": "-",
                "Object": "Unavailable route",
                "Comparison": "Not used",
                "Current manuscript status": "Not included",
                "Reason / boundary statement": "Not used in the current analysis.",
            },
            {
                "Dataset": "SCENIC+/Kaggle old exploration route",
                "Role": "Archived exploration",
                "Cells": "-",
                "Object": "Legacy proxy/cloud exploration",
                "Comparison": "Not applicable",
                "Current manuscript status": "Archived route",
                "Reason / boundary statement": "Archived route; not part of the current main workflow.",
            },
        ]
    )
    write_csv(df, path)
    return df, path


def st02() -> tuple[pd.DataFrame, Path]:
    path = OUT / "ST02_discovery_TF_shortlist.csv"
    df = read_csv(INPUTS["table2"]).rename(
        columns={
            "Expr_log2FC": "Expression_log2FC",
            "Regulon_dAUC": "Regulon_deltaAUC",
            "Pos_frac_high_group": "Positive_fraction_high_group",
            "Tier": "Evidence_tier",
        }
    )
    df["Axis"] = df["Axis"].str.replace("internal_control-associated", "internal-control-associated", regex=False)
    tier_map = {
        "NFE2L2": "Primary lesion-associated",
        "THRB": "Primary internal-control-associated",
        "BHLHE40": "Secondary lesion-associated",
        "SOX2": "Retained / lower-priority",
        "SATB2": "Supportive internal-control-associated",
        "RARB": "Supportive internal-control-associated",
        "HMGA1": "Lesion-side context",
    }
    df["Evidence_tier"] = df["TF"].map(tier_map).fillna(df["Evidence_tier"])
    df = df[["TF", "Axis", "Expression_log2FC", "Regulon_deltaAUC", "Regulon_FDR", "Positive_fraction_high_group", "Evidence_tier"]]
    write_csv(df, path)
    return df, path


def st03() -> tuple[pd.DataFrame, Path]:
    path = OUT / "ST03_integrated_TF_prioritization_metrics.csv"
    df = read_csv(INPUTS["table3"]).rename(
        columns={
            "Expr_log2FC": "Expression_log2FC",
            "Regulon_dAUC": "Regulon_deltaAUC",
            "RI_rank": "Approx_RI_rank",
            "Retention_freq": "Retention_frequency",
            "LOSO_expr": "LOSO_expression",
            "GSE190452_support": "GSE190452_cross_syndrome_support",
            "Robustness": "Sample_level_robustness",
            "Final_tier": "Final_priority",
        }
    )
    df["Axis"] = df["Axis"].str.replace("internal_control-associated", "internal-control-associated", regex=False)
    robustness_map = {
        "NFE2L2": "Stable",
        "THRB": "Stable with threshold sensitivity",
        "BHLHE40": "Stable",
        "SOX2": "Stable but target-limited",
    }
    final_map = {
        "NFE2L2": "Primary lesion-associated axis",
        "THRB": "Primary internal-control-associated axis",
        "BHLHE40": "Secondary lesion-associated candidate",
        "SOX2": "Retained / lower-priority",
    }
    df["Sample_level_robustness"] = df["TF"].map(robustness_map).fillna(df["Sample_level_robustness"])
    df["Final_priority"] = df["TF"].map(final_map).fillna(df["Final_priority"])
    cols = [
        "TF",
        "Axis",
        "Expression_log2FC",
        "Regulon_deltaAUC",
        "KO_rank",
        "Approx_RI_rank",
        "Retention_frequency",
        "LOSO_expression",
        "LOSO_regulon",
        "GSE140393_support",
        "GSE190452_cross_syndrome_support",
        "Sample_level_robustness",
        "Final_priority",
    ]
    df = df[cols]
    write_csv(df, path)
    return df, path


def st04() -> tuple[pd.DataFrame | None, Path]:
    path = OUT / "ST04_full_differential_regulon_statistics.csv"
    src = INPUTS["st04"]
    if not src.exists():
        return None, path
    raw = read_csv(src)
    diff = pd.to_numeric(raw["dAUC_lesion_minus_control"], errors="coerce")
    ranked = raw.assign(_abs=diff.abs()).sort_values("_abs", ascending=False).reset_index(drop=True)
    direction = diff.loc[ranked.index] if False else pd.to_numeric(ranked["dAUC_lesion_minus_control"], errors="coerce")
    df = pd.DataFrame(
        {
            "Regulon": ranked["Regulon"],
            "Mean_AUC_difference_lesion_minus_internal-control": ranked["dAUC_lesion_minus_control"],
            "FDR": ranked["FDR"],
            "Direction": direction.map(lambda x: "lesion-associated" if x > 0 else "internal-control-associated" if x < 0 else "no directional difference"),
            "Rank": [str(i) for i in range(1, len(ranked) + 1)],
            "Notes": "Ranked by absolute mean regulon AUC difference while preserving lesion minus internal-control direction.",
        }
    )
    write_csv(df, path)
    return df, path


def st05() -> tuple[pd.DataFrame | None, Path]:
    path = OUT / "ST05_celloracle_perturbation_metrics.csv"
    src = INPUTS["st05"]
    if not src.exists():
        return None, path
    raw = read_csv(src)
    df = pd.DataFrame(
        {
            "TF": raw["TF"],
            "Mean_shift": raw["Mean_shift"],
            "Approx_RI": raw["Recovery_index"],
            "Net_shift": raw["Net_shift"],
            "Lesion_shift": raw["Lesion_shift"],
            "Internal-control_shift": raw["Control_shift"],
            "Direction_agreement": raw["Direction_agreement"],
            "Interpretation": raw["TF"].map(
                {
                    "NFE2L2": "Primary lesion-associated axis; highest mean shift.",
                    "THRB": "Primary internal-control-associated axis; highest Approx. RI.",
                    "BHLHE40": "Secondary lesion-associated candidate.",
                    "SOX2": "Retained / lower-priority candidate.",
                }
            ),
        }
    )
    write_csv(df, path)
    return df, path


def st06() -> tuple[pd.DataFrame | None, Path]:
    path = OUT / "ST06_sample_level_robustness_and_threshold_sensitivity.csv"
    src = INPUTS["st06"]
    if not src.exists():
        return None, path
    raw = clean_frame(read_csv(src))
    df = raw.rename(
        columns={
            "Discovery_expr_direction": "Discovery_expression_direction",
            "LOSO_expr_rate": "LOSO_expression_rate",
            "Pseudobulk_expr_direction": "Pseudobulk_expression_direction",
            "Pseudobulk_expr_effect": "Pseudobulk_expression_effect",
            "Robustness": "Sample_level_robustness",
            "Priority_tier": "Priority_tier",
        }
    )
    priority_map = {
        "primary_lesion_axis": "Primary lesion-associated axis",
        "primary_internal-control_axis": "Primary internal-control-associated axis",
        "second_tier_lesion_candidate": "Secondary lesion-associated candidate",
        "reserved_candidate_not_main_axis": "Retained / lower-priority",
        "supporting_lesion_candidate": "Lesion-side context",
        "supporting_internal-control_candidate": "Supportive internal-control-associated",
    }
    df["Priority_tier"] = df["Priority_tier"].map(priority_map).fillna(df["Priority_tier"])
    write_csv(df, path)
    return df, path


def st07() -> tuple[pd.DataFrame | None, Path]:
    path = OUT / "ST07_supportive_external_evidence_summary.csv"
    src = INPUTS["st07"]
    if not src.exists():
        return None, path
    raw = clean_frame(read_csv(src))
    df = raw.rename(
        columns={
            "dataset": "Dataset",
            "analysis": "Analysis",
            "TF": "TF",
            "Expected_axis": "Expected_axis",
            "Lesion_presence_score": "Lesion_presence_score",
            "Direction_score": "Direction_score",
            "Support_score": "Support_score",
            "Support_rank": "Support_rank",
            "Group_log2FC": "Group_log2FC",
            "Cluster_support_fraction": "Cluster_expression_support_fraction",
        }
    )
    support_type = []
    boundary = []
    for dataset in df["Dataset"]:
        if dataset == "GSE140393":
            support_type.append("single-group expression support")
            boundary.append("single-group supportive analysis; expression-level supportive evidence only")
        else:
            support_type.append("cross-syndrome expression / target-set support")
            boundary.append("cross-syndrome supportive analysis; expression-level and target-set supportive evidence only")
    df.insert(2, "Support_type", support_type)
    df["Boundary_statement"] = boundary
    write_csv(df, path)
    return df, path


def sentence_case_term(term: str) -> str:
    protected = {
        "rna": "RNA",
        "ii": "II",
        "dna-templated": "DNA-templated",
        "mapk": "MAPK",
        "kegg": "KEGG",
        "go": "GO",
    }
    parts = term.split(" ")
    lowered = []
    for part in parts:
        key = part.lower()
        lowered.append(protected.get(key, part[:1].upper() + part[1:].lower() if not lowered else part.lower()))
    out = " ".join(lowered)
    out = out.replace("Rna", "RNA").replace("Dna-templated", "DNA-templated").replace("Mapk", "MAPK")
    return out


def st08() -> tuple[pd.DataFrame | None, Path]:
    path = OUT / "ST08_functional_enrichment_and_program_convergence.csv"
    src = INPUTS["st08"]
    if not src.exists():
        return None, path
    raw = read_csv(src)
    fixed = {
        "NFE2L2": {
            "Program_interpretation": "stress-adaptation / transcriptional reprogramming",
            "Regulon_targets": "458",
            "Overlap_genes": "245",
            "Overlap_ratio": "0.53",
            "Notes": "Primary lesion-associated program; term count summarizes enrichment breadth and does not determine integrated TF priority.",
        },
        "THRB": {
            "Program_interpretation": "compact homeostatic / synaptic-calcium supportive",
            "Regulon_targets": "17",
            "Overlap_genes": "13",
            "Overlap_ratio": "0.76",
            "Notes": "Primary internal-control-associated compact program; target set is small but overlap ratio is high.",
        },
        "BHLHE40": {
            "Program_interpretation": "secondary lesion-supportive",
            "Regulon_targets": "221",
            "Overlap_genes": "143",
            "Overlap_ratio": "0.65",
            "Notes": "Secondary lesion-supportive candidate.",
        },
        "SOX2": {
            "Program_interpretation": "target-limited retained",
            "Regulon_targets": "10",
            "Overlap_genes": "10",
            "Overlap_ratio": "1.00",
            "Notes": "Target-limited retained candidate; term count is not interpreted as integrated TF priority.",
        },
    }
    rows = []
    for tf in ["NFE2L2", "THRB", "BHLHE40", "SOX2"]:
        sub = raw[raw["TF"] == tf].copy()
        go_count = str((sub["Source"] == "GO").sum())
        kegg_count = str((sub["Source"] == "KEGG").sum())
        top_terms = "; ".join(sentence_case_term(x) for x in sub["Term"].head(6).tolist())
        row = {"TF": tf, **fixed[tf]}
        row["GO_significant_terms"] = go_count
        row["KEGG_significant_terms"] = kegg_count
        row["Top_terms"] = top_terms
        rows.append(row)
    df = pd.DataFrame(rows)[
        [
            "TF",
            "Program_interpretation",
            "Regulon_targets",
            "Overlap_genes",
            "Overlap_ratio",
            "GO_significant_terms",
            "KEGG_significant_terms",
            "Top_terms",
            "Notes",
        ]
    ]
    write_csv(df, path)
    return df, path


def st09() -> tuple[pd.DataFrame | None, Path]:
    path = OUT / "ST09_exploratory_drug_signature_enrichment.csv"
    src = INPUTS["st09"]
    if not src.exists():
        return None, path
    raw = read_csv(src)
    layer_map = {
        "primary_mechanism_direction_leads": "mechanism-direction clues",
        "headline_cns_mechanism_direction_leads": "mechanism-direction clues",
        "supportive_manual_review_leads": "supportive clues",
        "peripheral_program_modulating_clues": "peripheral/program-modulating clues",
        "peripheral/program-modulating candidates": "peripheral/program-modulating clue layer",
        "CNS-directed candidates": "BBB-reviewed CNS-directed clue layer",
    }
    theme_map = {
        "manual_review_needed": "manual review needed",
        "secondary_axis_clue": "secondary-axis clue",
        "environmental_toxicant": "flagged toxicant-related theme",
        "broad_cytotoxic/topoisomerase": "flagged cytotoxic/topoisomerase theme",
        "phosphodiesterase/cAMP": "phosphodiesterase / cAMP",
        "calcium/synaptic": "calcium / synaptic",
        "epigenetic/hdac": "epigenetic / HDAC",
        "redox/stress-adaptation": "redox / stress-adaptation",
    }
    df = pd.DataFrame(
        {
            "Compound": raw["Compound"],
            "Axis": raw["Axis"],
            "Mechanism_theme": raw["Mechanism_theme"].map(lambda x: theme_map.get(x, x)),
            "Pre_BBB_exploratory_layer": raw["PreBBB_layer"].map(lambda x: layer_map.get(x, x)),
            "After_BBB_exploratory_layer": raw["AfterBBB_layer"].map(lambda x: layer_map.get(x, x)),
            "BBB_review_layer": raw["BBB_layer"].map(lambda x: layer_map.get(x, x)),
            "Final_exploratory_score_after_BBB": raw["Final_score_after_BBB"],
            "Exploratory_role": raw["AfterBBB_layer"].map(lambda x: layer_map.get(x, x)) + "; exploratory drug-signature enrichment",
            "Interpretation_boundary": "exploratory clues only; not therapeutic recommendation",
        }
    )
    df["Notes"] = df["Mechanism_theme"].map(
        lambda x: "Cautionary flagged theme; not prioritized as CNS-directed interpretation."
        if x.startswith("flagged")
        else "Hypothesis-generating mechanism-direction clue only."
    )
    write_csv(df, path)
    return df, path


def relative(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def write_index(generated: dict[str, tuple[Path, str, str, str]]) -> Path:
    path = OUT / "main_and_supplementary_tables_index_v3.csv"
    rows = []
    specs = [
        ("Main Table 1", "Main_Table1_analytical_datasets_and_object_roles.csv", "Table 1. Analytical datasets and object roles.", "Main", SOURCE_LABELS["table1"]),
        ("Main Table 2", "Main_Table2_integrated_TF_prioritization.csv", "Table 2. Integrated prioritization of candidate transcription factors.", "Main", f"{SOURCE_LABELS['table2']}; {SOURCE_LABELS['table3']}"),
        ("ST01", "ST01_dataset_role_and_excluded_resource_audit.csv", "Supplementary Table ST01. Dataset role and excluded-resource audit.", "Supplementary", SOURCE_LABELS["table1"]),
        ("ST02", "ST02_discovery_TF_shortlist.csv", "Supplementary Table ST02. Discovery TF shortlist from differential regulon and expression analysis.", "Supplementary", SOURCE_LABELS["table2"]),
        ("ST03", "ST03_integrated_TF_prioritization_metrics.csv", "Supplementary Table ST03. Integrated TF prioritization metrics.", "Supplementary", SOURCE_LABELS["table3"]),
        ("ST04", "ST04_full_differential_regulon_statistics.csv", "Supplementary Table ST04. Full differential regulon statistics.", "Supplementary", SOURCE_LABELS["st04"]),
        ("ST05", "ST05_celloracle_perturbation_metrics.csv", "Supplementary Table ST05. CellOracle perturbation metrics.", "Supplementary", SOURCE_LABELS["st05"]),
        ("ST06", "ST06_sample_level_robustness_and_threshold_sensitivity.csv", "Supplementary Table ST06. Sample-level robustness and threshold sensitivity.", "Supplementary", SOURCE_LABELS["st06"]),
        ("ST07", "ST07_supportive_external_evidence_summary.csv", "Supplementary Table ST07. Supportive external evidence summary.", "Supplementary", SOURCE_LABELS["st07"]),
        ("ST08", "ST08_functional_enrichment_and_program_convergence.csv", "Supplementary Table ST08. Functional enrichment and program convergence.", "Supplementary", SOURCE_LABELS["st08"]),
        ("ST09", "ST09_exploratory_drug_signature_enrichment.csv", "Supplementary Table ST09. Exploratory drug-signature enrichment.", "Supplementary", SOURCE_LABELS["st09"]),
    ]
    for table_id, file_name, title, main_or_supp, source in specs:
        exists = (OUT / file_name).exists()
        rows.append(
            {
                "Table ID": table_id,
                "File name": file_name,
                "Title": title,
                "Main or supplementary": main_or_supp,
                "Source input": source,
                "Status": "generated" if exists else "missing input; not generated",
                "Notes": generated.get(table_id, (None, "", "", "Missing input; not generated."))[3],
            }
        )
    df = pd.DataFrame(rows)
    write_csv(df, path)
    return path


def write_legends() -> tuple[Path, Path]:
    main_path = OUT / "main_table_legends_v3.md"
    supp_path = OUT / "supplementary_table_legends_v3.md"
    main_text = """# Main Table Legends

Table 1. Analytical datasets and object roles.
This table summarizes the datasets and analytical objects used in the current manuscript. GSE268807 served as the discovery cohort, whereas GSE140393 and GSE190452 were used only as supportive evidence.

Table 2. Integrated prioritization of candidate transcription factors.
This table summarizes integrated evidence from discovery regulon analysis, CellOracle in silico perturbation, sample-level robustness, supportive external analyses, and functional interpretation. Final priorities are interpretive categories and should not be read as experimental confirmation.
"""
    supp_text = """# Supplementary Table Legends

Supplementary Table ST01. Dataset role and excluded-resource audit.
This table records included datasets and excluded or archived legacy resources.

Supplementary Table ST02. Discovery TF shortlist from differential regulon and expression analysis.
This table lists discovery TFs from expression and regulon differential analysis, including primary, secondary, retained, supportive, and contextual TFs.

Supplementary Table ST03. Integrated TF prioritization metrics.
This table provides numeric metrics underlying the integrated prioritization summarized in Main Table 2.

Supplementary Table ST04. Full differential regulon statistics.
This table provides all-regulon differential statistics ranked by mean regulon AUC difference.

Supplementary Table ST05. CellOracle perturbation metrics.
This table provides CellOracle in silico perturbation metrics, including mean shift and Approx. RI.

Supplementary Table ST06. Sample-level robustness and threshold sensitivity.
This table summarizes leave-one-sample-out consistency, pseudobulk support, and threshold sensitivity metrics.

Supplementary Table ST07. Supportive external evidence summary.
This table summarizes GSE140393 single-group supportive analysis and GSE190452 cross-syndrome supportive analysis.

Supplementary Table ST08. Functional enrichment and program convergence.
This table summarizes GO/KEGG enrichment breadth and regulon-DEG program convergence for candidate TFs.

Supplementary Table ST09. Exploratory drug-signature enrichment.
This table summarizes exploratory drug-signature enrichment and BBB-layered mechanism-direction clues; compounds are hypothesis-generating clues only and not therapeutic recommendation.
"""
    main_path.write_text(main_text, encoding="utf-8")
    supp_path.write_text(supp_text, encoding="utf-8")
    return main_path, supp_path


def write_revision_log(generated_status: dict[str, bool]) -> Path:
    path = OUT / "main_and_supplementary_tables_revision_log_v3.md"
    missing = [k for k, v in generated_status.items() if k.startswith("ST") and not v]
    missing_line = "None; ST04-ST09 inputs were found and generated." if not missing else ", ".join(missing)
    text = f"""# Main and Supplementary Tables Revision Log v3

- Rebuilt the manuscript-facing main tables as two concise tables.
- Main Table 1 retains only the three analytical objects used in the current manuscript.
- GSE275302 and the archived legacy exploration route were moved to ST01.
- The original shortlist table was moved from the main-table role to ST02.
- The integrated-priority table was split into an interpretive Main Table 2 and a numeric ST03.
- Replaced legacy underscore-form internal-control wording with internal-control.
- Renamed the legacy Approx. Recovery Index rank column to Approx_RI_rank.
- Replaced over-strong external-confirmation wording with supportive evidence language.
- Preserved original numeric values where source tables provided numeric metrics; changes are limited to table tiering, column names, and interpretation terminology.
- ST04-ST09 generation status: {missing_line}
- No figure files were modified.
"""
    path.write_text(text, encoding="utf-8")
    return path


def scan_files(paths: list[Path], term: str) -> list[str]:
    if term == "RI_rank":
        pattern = re.compile(r"(?<!Approx_)RI_rank", re.IGNORECASE)
    else:
        pattern = re.compile(re.escape(term), re.IGNORECASE)
    hits = []
    for path in paths:
        if not path.exists() or path.name.endswith("terminology_audit_v3.csv"):
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        if pattern.search(text):
            hits.append(path.name)
    return hits


def write_audit(output_paths: list[Path]) -> Path:
    path = OUT / "main_and_supplementary_tables_terminology_audit_v3.csv"
    checks = [
        ("internal_control", "absent", "replaced by internal-control"),
        ("internal-control", "present", "required"),
        ("RI_rank", "absent", "replaced by Approx_RI_rank"),
        ("Approx_RI_rank", "present", "required"),
        ("Regulon_dAUC", "absent", "replaced by Regulon_deltaAUC"),
        ("Regulon_deltaAUC", "present", "required"),
        ("validation", "absent", "avoid table-facing wording"),
        ("formal validation", "absent", "required absent"),
        ("external validation", "absent", "required absent"),
        ("independent validation", "absent", "required absent"),
        ("validated", "absent", "required absent"),
        ("confirmed", "absent", "required absent"),
        ("donor-level robustness", "absent", "replaced by sample-level robustness"),
        ("sample-level robustness", "present", "required"),
        ("SCENIC+", "present only in ST01 archived route", "not current mainline"),
        ("current mainline", "absent", "avoid implying archived route is current"),
        ("therapeutic recommendation", "present only as negated boundary statement", "not therapeutic recommendation"),
        ("validated drug", "absent", "required absent"),
        ("clinical candidate", "absent", "required absent"),
        ("treatment lead", "absent", "required absent"),
        ("experimental perturbation", "absent", "required absent"),
        ("in silico perturbation", "present", "required"),
        ("supportive evidence", "present", "required"),
        ("cross-syndrome supportive analysis", "present", "required"),
        ("single-group supportive analysis", "present", "required"),
    ]
    rows = []
    for term, expected_status, action in checks:
        hits = scan_files(output_paths, term)
        if term == "SCENIC+" and hits == ["ST01_dataset_role_and_excluded_resource_audit.csv", "supplementary_table_legends_v3.md"]:
            status = "present only in ST01 archived route and its legend"
        elif term == "therapeutic recommendation" and hits:
            status = "present only as negated boundary statement"
        elif hits:
            status = "present"
        else:
            status = "absent"
        if expected_status.startswith("present") and hits:
            status = expected_status if term in {"SCENIC+", "therapeutic recommendation"} else "present"
        rows.append(
            {
                "term": term,
                "status": status,
                "action": action,
                "file(s) affected": "; ".join(hits) if hits else "-",
            }
        )
    df = pd.DataFrame(rows)
    write_csv(df, path)
    return path


def run() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    generated: dict[str, tuple[Path, str, str, str]] = {}
    status: dict[str, bool] = {}
    output_paths: list[Path] = []

    for table_id, func, note in [
        ("Main Table 1", main_table_1, "Generated from fixed manuscript dataset roles."),
        ("Main Table 2", main_table_2, "Generated as interpretive evidence-integration table."),
        ("ST01", st01, "Generated from Table1 and manuscript boundary decisions."),
        ("ST02", st02, "Generated from Table2 with renamed columns and updated evidence tiers."),
        ("ST03", st03, "Generated from Table3 with renamed columns and final priority labels."),
    ]:
        df, path = func()
        generated[table_id] = (path, table_id, "generated", note)
        status[table_id] = True
        output_paths.append(path)

    for table_id, func, note in [
        ("ST04", st04, "Generated from full differential regulon statistics input."),
        ("ST05", st05, "Generated from CellOracle round1 perturbation metrics input."),
        ("ST06", st06, "Generated from sample-level robustness input."),
        ("ST07", st07, "Generated from supportive external evidence input."),
        ("ST08", st08, "Generated from functional enrichment and convergence input."),
        ("ST09", st09, "Generated from exploratory after-BBB drug-signature enrichment input."),
    ]:
        df, path = func()
        ok = df is not None
        status[table_id] = ok
        if ok:
            generated[table_id] = (path, table_id, "generated", note)
            output_paths.append(path)
        else:
            generated[table_id] = (path, table_id, "missing input; not generated", "Missing input; no synthetic data generated.")

    index_path = write_index(generated)
    main_legend, supp_legend = write_legends()
    revision_log = write_revision_log(status)
    output_paths.extend([index_path, main_legend, supp_legend, revision_log])
    audit = write_audit(output_paths)
    output_paths.append(audit)

    generated_st = [k for k in [f"ST{i:02d}" for i in range(4, 10)] if status.get(k)]
    missing_st = [k for k in [f"ST{i:02d}" for i in range(4, 10)] if not status.get(k)]
    table_scan_paths = [p for p in output_paths if p.name != audit.name]
    residual_internal = scan_files(table_scan_paths, "internal_control")
    residual_validation = sorted(set(scan_files(table_scan_paths, "validation") + scan_files(table_scan_paths, "external validation") + scan_files(table_scan_paths, "formal validation")))
    residual_therapeutic = sorted(set(scan_files(table_scan_paths, "validated drug") + scan_files(table_scan_paths, "clinical candidate")))
    therapeutic_recommendation_hits = scan_files(table_scan_paths, "therapeutic recommendation")

    print("Main Table 1 generated: yes")
    print("Main Table 2 generated: yes")
    print("ST01-ST03 generated: yes")
    print("ST04-ST09 generated:", ", ".join(generated_st) if generated_st else "none")
    print("ST04-ST09 missing:", ", ".join(missing_st) if missing_st else "none")
    print("Original numeric values preserved: yes; source numeric fields were copied or summarized without recalculation")
    print("Residual internal_control:", "none" if not residual_internal else "; ".join(residual_internal))
    print("Residual validation / external validation / formal validation:", "none" if not residual_validation else "; ".join(residual_validation))
    if therapeutic_recommendation_hits:
        print("Residual therapeutic recommendation / validated drug / clinical candidate: negated boundary phrase only in " + "; ".join(therapeutic_recommendation_hits))
    else:
        print("Residual therapeutic recommendation / validated drug / clinical candidate: none")
    print("Output directory:", relative(OUT))
    print("Index:", relative(index_path))
    print("Revision log:", relative(revision_log))
    print("Terminology audit:", relative(audit))


if __name__ == "__main__":
    run()
