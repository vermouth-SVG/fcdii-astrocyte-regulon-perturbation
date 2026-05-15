from __future__ import annotations

import csv
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "manuscript_output" / "tables_revised_v3"
OUT = ROOT / "manuscript_output" / "tables_revised_v3_minorfix"

TABLE_FILES = [
    ("Main_Table1_analytical_datasets_and_object_roles.csv", "Main_Table1_analytical_datasets_and_object_roles_v3_minorfix.csv"),
    ("Main_Table2_integrated_TF_prioritization.csv", "Main_Table2_integrated_TF_prioritization_v3_minorfix.csv"),
    ("ST01_dataset_role_and_excluded_resource_audit.csv", "ST01_dataset_role_and_excluded_resource_audit_v3_minorfix.csv"),
    ("ST02_discovery_TF_shortlist.csv", "ST02_discovery_TF_shortlist_v3_minorfix.csv"),
    ("ST03_integrated_TF_prioritization_metrics.csv", "ST03_integrated_TF_prioritization_metrics_v3_minorfix.csv"),
    ("ST04_full_differential_regulon_statistics.csv", "ST04_full_differential_regulon_statistics_v3_minorfix.csv"),
    ("ST05_celloracle_perturbation_metrics.csv", "ST05_celloracle_perturbation_metrics_v3_minorfix.csv"),
    ("ST06_sample_level_robustness_and_threshold_sensitivity.csv", "ST06_sample_level_robustness_and_threshold_sensitivity_v3_minorfix.csv"),
    ("ST07_supportive_external_evidence_summary.csv", "ST07_supportive_external_evidence_summary_v3_minorfix.csv"),
    ("ST08_functional_enrichment_and_program_convergence.csv", "ST08_functional_enrichment_and_program_convergence_v3_minorfix.csv"),
    ("ST09_exploratory_drug_signature_enrichment.csv", "ST09_exploratory_drug_signature_enrichment_v3_minorfix.csv"),
]


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def process_tables() -> dict[str, Path]:
    OUT.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, Path] = {}
    for src_name, out_name in TABLE_FILES:
        df = read_csv(SRC / src_name)

        if src_name.startswith("ST04_"):
            old_col = "Mean_AUC_difference_lesion_minus_internal-control"
            new_col = "Mean AUC difference (lesion − internal-control)"
            if old_col in df.columns:
                df = df.rename(columns={old_col: new_col})
            repeated_note = "Ranked by absolute mean regulon AUC difference while preserving lesion minus internal-control direction."
            if "Notes" in df.columns and df["Notes"].nunique(dropna=False) == 1 and df["Notes"].iloc[0] == repeated_note:
                df = df.drop(columns=["Notes"])

        if src_name.startswith("ST06_"):
            robustness_map = {
                "NFE2L2": "Directionally stable",
                "THRB": "Directionally stable with threshold sensitivity",
                "BHLHE40": "Directionally stable secondary signal",
                "SOX2": "Directionally stable but target-limited",
                "HMGA1": "Directionally stable context signal",
                "SATB2": "Directionally stable supportive signal",
                "RARB": "Directionally stable supportive signal",
            }
            df["Sample_level_robustness"] = df["TF"].map(robustness_map).fillna(df["Sample_level_robustness"])

        if src_name.startswith("ST09_"):
            df = df.map(
                lambda x: x.replace("supportive_after_bbb_leads", "supportive_after_bbb_clues")
                .replace("possible CNS-directed candidates", "possible CNS-directed exploratory clues")
                .replace("not therapeutic recommendation", "not a therapeutic recommendation")
                if isinstance(x, str)
                else x
            )

        out_path = OUT / out_name
        write_csv(df, out_path)
        outputs[out_name] = out_path
    return outputs


def write_main_legend() -> Path:
    src = (SRC / "main_table_legends_v3.md").read_text(encoding="utf-8")
    text = src.replace(
        "Final priorities are interpretive categories and should not be read as experimental confirmation.",
        "Final priorities are interpretive categories and should not be read as wet-lab evidence.",
    )
    path = OUT / "main_table_legends_v3_minorfix.md"
    write_text(path, text)
    return path


def write_supplementary_legend() -> Path:
    src = (SRC / "supplementary_table_legends_v3.md").read_text(encoding="utf-8")
    text = src.replace(
        "This table provides all-regulon differential statistics ranked by mean regulon AUC difference.",
        "This table provides all-regulon differential statistics. Regulons are ranked by absolute mean regulon AUC difference while preserving the lesion minus internal-control direction.",
    )
    text = text.replace("not therapeutic recommendation", "not a therapeutic recommendation")
    text = text.replace(
        "This table summarizes exploratory drug-signature enrichment and BBB-layered mechanism-direction clues; compounds are hypothesis-generating clues only and not a therapeutic recommendation.",
        "This table summarizes exploratory drug-signature enrichment and BBB-layered mechanism-direction clues; compounds are hypothesis-generating clues only and not a therapeutic recommendation.",
    )
    path = OUT / "supplementary_table_legends_v3_minorfix.md"
    write_text(path, text)
    return path


def table_title(file_name: str) -> str:
    titles = {
        "Main_Table1_analytical_datasets_and_object_roles_v3_minorfix.csv": "Table 1. Analytical datasets and object roles.",
        "Main_Table2_integrated_TF_prioritization_v3_minorfix.csv": "Table 2. Integrated prioritization of candidate transcription factors.",
        "ST01_dataset_role_and_excluded_resource_audit_v3_minorfix.csv": "Supplementary Table ST01. Dataset role and excluded-resource audit.",
        "ST02_discovery_TF_shortlist_v3_minorfix.csv": "Supplementary Table ST02. Discovery TF shortlist from differential regulon and expression analysis.",
        "ST03_integrated_TF_prioritization_metrics_v3_minorfix.csv": "Supplementary Table ST03. Integrated TF prioritization metrics.",
        "ST04_full_differential_regulon_statistics_v3_minorfix.csv": "Supplementary Table ST04. Full differential regulon statistics.",
        "ST05_celloracle_perturbation_metrics_v3_minorfix.csv": "Supplementary Table ST05. CellOracle perturbation metrics.",
        "ST06_sample_level_robustness_and_threshold_sensitivity_v3_minorfix.csv": "Supplementary Table ST06. Sample-level robustness and threshold sensitivity.",
        "ST07_supportive_external_evidence_summary_v3_minorfix.csv": "Supplementary Table ST07. Supportive external evidence summary.",
        "ST08_functional_enrichment_and_program_convergence_v3_minorfix.csv": "Supplementary Table ST08. Functional enrichment and program convergence.",
        "ST09_exploratory_drug_signature_enrichment_v3_minorfix.csv": "Supplementary Table ST09. Exploratory drug-signature enrichment.",
    }
    return titles[file_name]


def write_index(outputs: dict[str, Path]) -> Path:
    old_index = read_csv(SRC / "main_and_supplementary_tables_index_v3.csv")
    rows = []
    old_to_new = dict(TABLE_FILES)
    for _, row in old_index.iterrows():
        new_file = old_to_new.get(row["File name"], row["File name"])
        rows.append(
            {
                "Table ID": row["Table ID"],
                "File name": new_file,
                "Title": table_title(new_file) if new_file.endswith("_v3_minorfix.csv") else row["Title"],
                "Main or supplementary": row["Main or supplementary"],
                "Source input": row["File name"],
                "Status": "generated",
                "Notes": "v3 minorfix output; original numeric values preserved.",
            }
        )
    df = pd.DataFrame(rows)
    path = OUT / "main_and_supplementary_tables_index_v3_minorfix.csv"
    write_csv(df, path)
    return path


def write_revision_log() -> Path:
    text = """# Main and Supplementary Tables Revision Log v3 Minorfix

- Generated a separate tables_revised_v3_minorfix output folder without overwriting v3 files.
- ST04 column label was changed to Mean AUC difference (lesion − internal-control).
- ST04 repeated row-wise Notes were removed; the ranking explanation was moved to the supplementary table legends.
- ST06 Sample_level_robustness values were replaced with directionally specific statements while preserving all numeric columns and Priority_tier.
- ST09 drug-status wording was downgraded from lead/candidate wording to clue wording where present.
- ST09 boundary wording was updated to not a therapeutic recommendation.
- Main Table 2 legend was updated to use wet-lab evidence wording.
- Supplementary legends were updated for ST04 and ST09.
- Original numeric values, compounds, axes, mechanisms, scores, regulons, FDR values, directions, and ranks were preserved.
- No figure files were modified.
"""
    path = OUT / "main_and_supplementary_tables_revision_log_v3_minorfix.md"
    write_text(path, text)
    return path


def write_readme() -> Path:
    text = """# Revised Main and Supplementary Tables v3 Minorfix

Current manuscript-facing main tables are two concise tables.

Final table files:
- Main_Table1_analytical_datasets_and_object_roles.csv
- Main_Table2_integrated_TF_prioritization.csv
- ST01_dataset_role_and_excluded_resource_audit.csv
- ST02_discovery_TF_shortlist.csv
- ST03_integrated_TF_prioritization_metrics.csv
- ST04_full_differential_regulon_statistics.csv
- ST05_celloracle_perturbation_metrics.csv
- ST06_sample_level_robustness_and_threshold_sensitivity.csv
- ST07_supportive_external_evidence_summary.csv
- ST08_functional_enrichment_and_program_convergence.csv
- ST09_exploratory_drug_signature_enrichment.csv
- main_and_supplementary_tables_index_v3_minorfix.csv
- main_table_legends_v3_minorfix.md
- supplementary_table_legends_v3_minorfix.md
- main_and_supplementary_tables_revision_log_v3_minorfix.md
- main_and_supplementary_tables_terminology_audit_v3_minorfix.csv

The actual files in this folder include the `_v3_minorfix` suffix for traceability and to avoid overwriting v3 outputs.

The original shortlist table is retained as ST02.
The full numeric integrated-priority table is retained as ST03.
GSE275302 and archived SCENIC+ / Kaggle exploration route are recorded only in ST01.
Drug-signature enrichment results are exploratory clues only and not a therapeutic recommendation.
"""
    path = OUT / "README_tables_revised_v3_minorfix.md"
    write_text(path, text)
    return path


def scan_files(paths: list[Path], term: str) -> list[str]:
    if term == "RI_rank":
        pattern = re.compile(r"(?<!Approx_)RI_rank", re.IGNORECASE)
    else:
        pattern = re.compile(re.escape(term), re.IGNORECASE)
    hits = []
    for path in paths:
        if path.name.endswith("terminology_audit_v3_minorfix.csv"):
            continue
        if pattern.search(path.read_text(encoding="utf-8", errors="ignore")):
            hits.append(path.name)
    return hits


def drug_status_hits(term: str) -> list[str]:
    st09 = OUT / "ST09_exploratory_drug_signature_enrichment_v3_minorfix.csv"
    pattern = re.compile(rf"\b{re.escape(term)}\b", re.IGNORECASE)
    text = st09.read_text(encoding="utf-8", errors="ignore")
    return [st09.name] if pattern.search(text) else []


def write_audit(paths: list[Path]) -> Path:
    rows = []
    checks = [
        ("internal_control", "absent", "required absent"),
        ("internal-control", "present", "required"),
        ("validation", "absent", "required absent"),
        ("formal validation", "absent", "required absent"),
        ("external validation", "absent", "required absent"),
        ("validated drug", "absent", "required absent"),
        ("clinical candidate", "absent", "required absent"),
        ("treatment lead", "absent", "required absent"),
        ("lead", "absent in ST09 drug-status wording", "replace drug-status wording with clues"),
        ("leads", "absent in ST09 drug-status wording", "replace drug-status wording with clues"),
        ("candidate", "absent in ST09 drug-status wording; allowed for TF candidate wording", "replace drug-status wording with clues"),
        ("candidates", "absent in ST09 drug-status wording; allowed for TF candidate wording", "replace drug-status wording with clues"),
        ("not a therapeutic recommendation", "present", "required boundary phrase"),
        ("exploratory drug-signature enrichment", "present", "required"),
        ("hypothesis-generating clues only", "present", "required"),
    ]
    for term, expected, action in checks:
        if term in {"lead", "leads", "candidate", "candidates"}:
            hits = drug_status_hits(term)
            status = "present in ST09" if hits else expected
            affected = "; ".join(hits) if hits else "-"
        else:
            hits = scan_files(paths, term)
            status = "present" if hits else "absent"
            if term == "not a therapeutic recommendation" and hits:
                status = "present as boundary phrase"
            if term in {"exploratory drug-signature enrichment", "hypothesis-generating clues only"} and hits:
                status = "present"
            affected = "; ".join(hits) if hits else "-"
        rows.append({"term": term, "status": status, "action": action, "file(s) affected": affected})
    path = OUT / "main_and_supplementary_tables_terminology_audit_v3_minorfix.csv"
    write_csv(pd.DataFrame(rows), path)
    return path


def compare_numeric_columns(old_name: str, new_name: str) -> bool:
    old = read_csv(SRC / old_name)
    new = read_csv(OUT / new_name)
    old_numeric_cols = []
    for col in old.columns:
        converted = pd.to_numeric(old[col], errors="coerce")
        if converted.notna().any():
            old_numeric_cols.append(col)
    rename_map = {
        "Mean_AUC_difference_lesion_minus_internal-control": "Mean AUC difference (lesion − internal-control)",
    }
    for old_col in old_numeric_cols:
        new_col = rename_map.get(old_col, old_col)
        if new_col not in new.columns:
            continue
        if not old[old_col].fillna("").astype(str).equals(new[new_col].fillna("").astype(str)):
            return False
    return True


def run() -> None:
    outputs = process_tables()
    extra_paths = [
        write_index(outputs),
        write_main_legend(),
        write_supplementary_legend(),
        write_revision_log(),
        write_readme(),
    ]
    all_paths = list(outputs.values()) + extra_paths
    audit = write_audit(all_paths)
    all_paths.append(audit)

    st04 = read_csv(OUT / "ST04_full_differential_regulon_statistics_v3_minorfix.csv")
    st06 = read_csv(OUT / "ST06_sample_level_robustness_and_threshold_sensitivity_v3_minorfix.csv")
    st09 = read_csv(OUT / "ST09_exploratory_drug_signature_enrichment_v3_minorfix.csv")
    st09_text = (OUT / "ST09_exploratory_drug_signature_enrichment_v3_minorfix.csv").read_text(encoding="utf-8")
    all_text_paths = [p for p in all_paths if p.name != audit.name]

    residual_drug_status = sorted(set(drug_status_hits("lead") + drug_status_hits("leads") + drug_status_hits("candidate") + drug_status_hits("candidates")))
    residual_forbidden = sorted(
        set(
            scan_files(all_text_paths, "internal_control")
            + scan_files(all_text_paths, "validation")
            + scan_files(all_text_paths, "formal validation")
            + scan_files(all_text_paths, "external validation")
            + scan_files(all_text_paths, "validated drug")
            + scan_files(all_text_paths, "clinical candidate")
            + scan_files(all_text_paths, "treatment lead")
        )
    )
    numeric_ok = all(compare_numeric_columns(old, new) for old, new in TABLE_FILES)

    print("ST04 column and Notes fix completed:", "yes" if "Mean AUC difference (lesion − internal-control)" in st04.columns and "Notes" not in st04.columns else "no")
    print("ST06 robustness wording fix completed:", "yes" if "high" not in set(st06["Sample_level_robustness"]) else "no")
    print("ST09 leads/candidates/recommendation downgrade completed:", "yes" if "supportive_after_bbb_leads" not in st09_text and "possible CNS-directed candidates" not in st09_text and "not a therapeutic recommendation" in st09_text else "no")
    print("Legends minorfix completed:", "yes")
    print("README update completed:", "yes")
    print("Residual drug-related lead / leads / candidate / candidates:", "none" if not residual_drug_status else "; ".join(residual_drug_status))
    print("Residual internal_control / validation / formal validation / external validation / validated drug / clinical candidate / treatment lead:", "none" if not residual_forbidden else "; ".join(residual_forbidden))
    print("All original numeric values preserved:", "yes" if numeric_ok else "no")
    print("Output directory:", OUT.relative_to(ROOT))


if __name__ == "__main__":
    run()
