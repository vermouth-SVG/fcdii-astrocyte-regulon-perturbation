import argparse
import ast
from pathlib import Path
from typing import Iterable, List, Sequence

import pandas as pd


DEFAULT_INPUT = Path(r"<PROJECT_ROOT_WINDOWS>\output\regulons.csv")
DEFAULT_GMT = Path(
    r"<PROJECT_ROOT_WINDOWS>\external_validation_round2\pyscenic_recalc\output\discovery_regulons_4tf.gmt"
)
DEFAULT_SELECTED = Path(
    r"<PROJECT_ROOT_WINDOWS>\external_validation_round2\pyscenic_recalc\output\discovery_regulons_4tf_selected_rows.csv"
)
DEFAULT_TFS = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]


def normalize_level_name(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.startswith("Unnamed:"):
        return ""
    return text


def column_has_name(column: Sequence[object], target: str) -> bool:
    target_norm = target.strip().lower()
    return any(normalize_level_name(level).lower() == target_norm for level in column)


def find_column(columns: Iterable[Sequence[object]], target: str):
    for column in columns:
        if column_has_name(column, target):
            return column
    raise KeyError(f"Could not find column '{target}' in regulons.csv multi-level header.")


def parse_target_genes(raw_value: object) -> List[str]:
    if raw_value is None or (isinstance(raw_value, float) and pd.isna(raw_value)):
        return []
    text = str(raw_value).strip()
    if not text:
        return []
    parsed = ast.literal_eval(text)
    genes: List[str] = []
    if isinstance(parsed, list):
        for item in parsed:
            if isinstance(item, (list, tuple)) and item:
                gene = str(item[0]).strip()
                if gene:
                    genes.append(gene)
            elif item is not None:
                gene = str(item).strip()
                if gene:
                    genes.append(gene)
    else:
        gene = str(parsed).strip()
        if gene:
            genes.append(gene)
    deduped: List[str] = []
    seen = set()
    for gene in genes:
        if gene not in seen:
            deduped.append(gene)
            seen.add(gene)
    return deduped


def flatten_regulons_table(df: pd.DataFrame) -> pd.DataFrame:
    tf_col = find_column(df.columns, "TF")
    motif_col = find_column(df.columns, "MotifID")
    auc_col = find_column(df.columns, "AUC")
    nes_col = find_column(df.columns, "NES")
    annotation_col = find_column(df.columns, "Annotation")
    context_col = find_column(df.columns, "Context")
    target_genes_col = find_column(df.columns, "TargetGenes")

    flat = pd.DataFrame(
        {
            "TF": df[tf_col].astype(str).str.strip(),
            "MotifID": df[motif_col].astype(str).str.strip(),
            "AUC": pd.to_numeric(df[auc_col], errors="coerce"),
            "NES": pd.to_numeric(df[nes_col], errors="coerce"),
            "Annotation": df[annotation_col].astype(str),
            "Context": df[context_col].astype(str),
            "TargetGenes": df[target_genes_col].astype(str),
        }
    )
    return flat


def select_best_row_for_tf(tf_rows: pd.DataFrame) -> pd.Series:
    ranked = tf_rows.copy()
    annotation_lower = ranked["Annotation"].str.lower()
    context_lower = ranked["Context"].str.lower()
    ranked["annotation_direct"] = annotation_lower.str.contains("directly annotated") | annotation_lower.str.contains(
        "direct annotation"
    )
    ranked["context_activating"] = context_lower.str.contains("activating")
    ranked = ranked.sort_values(
        by=["annotation_direct", "context_activating", "NES", "AUC"],
        ascending=[False, False, False, False],
        na_position="last",
    ).reset_index(drop=True)
    return ranked.iloc[0]


def build_gmt_line(tf: str, selected: pd.Series, genes: List[str]) -> str:
    sign = "+" if bool(selected["context_activating"]) else "?"
    regulon_name = f"{tf}({sign})"
    description = (
        f"selected_from_discovery_regulons|motif={selected['MotifID']}|NES={selected['NES']:.6f}|AUC={selected['AUC']:.6f}"
    )
    return "\t".join([regulon_name, description, *genes])


def main() -> None:
    parser = argparse.ArgumentParser(description="Build 4-TF GMT from discovery regulons.csv with multi-level header.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--output-gmt", default=str(DEFAULT_GMT))
    parser.add_argument("--output-selected", default=str(DEFAULT_SELECTED))
    parser.add_argument("--tfs", nargs="+", default=DEFAULT_TFS)
    args = parser.parse_args()

    input_path = Path(args.input)
    output_gmt = Path(args.output_gmt)
    output_selected = Path(args.output_selected)
    output_gmt.parent.mkdir(parents=True, exist_ok=True)
    output_selected.parent.mkdir(parents=True, exist_ok=True)

    raw_df = pd.read_csv(input_path, header=[0, 1, 2])
    flat_df = flatten_regulons_table(raw_df)
    flat_df["TF_upper"] = flat_df["TF"].str.upper()

    selected_rows = []
    gmt_lines: List[str] = []
    missing_tfs: List[str] = []

    for tf in args.tfs:
        tf_upper = tf.upper()
        tf_rows = flat_df.loc[flat_df["TF_upper"] == tf_upper, :].copy()
        if tf_rows.empty:
            missing_tfs.append(tf)
            continue
        selected = select_best_row_for_tf(tf_rows)
        genes = parse_target_genes(selected["TargetGenes"])
        if not genes:
            raise ValueError(f"No target genes could be parsed for {tf} from motif row {selected['MotifID']}.")

        gmt_lines.append(build_gmt_line(tf, selected, genes))
        selected_rows.append(
            {
                "TF": tf,
                "regulon_name": f"{tf}(+)" if bool(selected["context_activating"]) else f"{tf}(?)",
                "MotifID": selected["MotifID"],
                "AUC": selected["AUC"],
                "NES": selected["NES"],
                "Annotation": selected["Annotation"],
                "Context": selected["Context"],
                "annotation_direct": bool(selected["annotation_direct"]),
                "context_activating": bool(selected["context_activating"]),
                "target_gene_count": len(genes),
                "target_genes_preview": ",".join(genes[:20]),
            }
        )

    if missing_tfs:
        raise ValueError(f"Missing TF rows in regulons.csv for: {', '.join(missing_tfs)}")

    with output_gmt.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(gmt_lines) + "\n")

    selected_df = pd.DataFrame(selected_rows)
    selected_df.to_csv(output_selected, index=False)

    print(f"Input regulons.csv: {input_path}")
    print(f"Selected TF count: {len(selected_rows)}")
    print(f"GMT written to: {output_gmt}")
    print(f"Selected rows table written to: {output_selected}")
    for row in selected_rows:
        print(
            f"{row['TF']}: {row['MotifID']} | direct={row['annotation_direct']} | "
            f"activating={row['context_activating']} | NES={row['NES']:.6f} | "
            f"AUC={row['AUC']:.6f} | targets={row['target_gene_count']}"
        )


if __name__ == "__main__":
    main()
