#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "external_validation_round2" / "pyscenic_recalc" / "input"
OUTPUT_DIR = ROOT / "external_validation_round2" / "pyscenic_recalc" / "output"
LOG_DIR = ROOT / "external_validation_round2" / "pyscenic_recalc" / "logs"
DB_DIR = ROOT / "db"
DOCKER_IMAGE = "aertslab/pyscenic:0.12.1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the round2 pySCENIC pipeline inside the existing Docker route.")
    parser.add_argument("--expression-csv", default=str(INPUT_DIR / "expression_for_pyscenic.csv"))
    parser.add_argument("--tf-list", default=str(DB_DIR / "allTFs_hg38.txt"))
    parser.add_argument(
        "--ranking-db",
        default=str(DB_DIR / "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"),
    )
    parser.add_argument(
        "--motif-annotation",
        default=str(DB_DIR / "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"),
    )
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--image", default=DOCKER_IMAGE)
    parser.add_argument("--force", action="store_true", help="Re-run even if output files already exist.")
    parser.add_argument("--dry-run", action="store_true", help="Only generate commands and reports, do not execute.")
    return parser.parse_args()


def to_mount_spec(path: Path) -> str:
    posix = path.resolve().as_posix()
    if len(posix) >= 2 and posix[1] == ":":
        posix = posix[0].lower() + posix[1:]
    return f"{posix}:/work"


def to_container_path(path: Path) -> str:
    relative = path.resolve().relative_to(ROOT.resolve())
    return f"/work/{relative.as_posix()}"


def ensure_parent_dirs() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)


def build_commands(args: argparse.Namespace) -> list[tuple[str, list[str], Path]]:
    expression_csv = Path(args.expression_csv)
    tf_list = Path(args.tf_list)
    ranking_db = Path(args.ranking_db)
    motif_annotation = Path(args.motif_annotation)

    grn_out = OUTPUT_DIR / "grn_adj.tsv"
    ctx_out = OUTPUT_DIR / "regulons.csv"
    auc_out = OUTPUT_DIR / "auc_mtx.csv"

    volume = to_mount_spec(ROOT)
    prefix = ["docker", "run", "--rm", "-v", volume, "-w", "/work", args.image]

    grn_cmd = prefix + [
        "pyscenic",
        "grn",
        to_container_path(expression_csv),
        to_container_path(tf_list),
        "-o",
        to_container_path(grn_out),
        "--num_workers",
        str(args.workers),
        "--method",
        "grnboost2",
        "--sparse",
    ]
    ctx_cmd = prefix + [
        "pyscenic",
        "ctx",
        to_container_path(grn_out),
        to_container_path(ranking_db),
        "--annotations_fname",
        to_container_path(motif_annotation),
        "--expression_mtx_fname",
        to_container_path(expression_csv),
        "--num_workers",
        str(args.workers),
        "--mode",
        "custom_multiprocessing",
        "-o",
        to_container_path(ctx_out),
    ]
    auc_cmd = prefix + [
        "pyscenic",
        "aucell",
        to_container_path(expression_csv),
        to_container_path(ctx_out),
        "-o",
        to_container_path(auc_out),
        "--num_workers",
        str(args.workers),
    ]
    return [
        ("01_grn", grn_cmd, grn_out),
        ("02_ctx", ctx_cmd, ctx_out),
        ("03_aucell", auc_cmd, auc_out),
    ]


def write_command_report(commands: list[tuple[str, list[str], Path]]) -> Path:
    report = LOG_DIR / "pyscenic_commands.txt"
    lines = ["round2 pySCENIC commands"]
    for step, cmd, expected in commands:
        lines.append(f"[{step}]")
        lines.append(" ".join(shlex.quote(x) for x in cmd))
        lines.append(f"expected_output: {expected}")
        lines.append("")
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


def run_step(step: str, cmd: list[str], expected_output: Path, force: bool) -> None:
    log_path = LOG_DIR / f"{step}.log"
    if expected_output.exists() and not force:
        log_path.write_text(f"Skipped {step} because output exists: {expected_output}\n", encoding="utf-8")
        print(f"Skip {step}: {expected_output.name} already exists")
        return
    with log_path.open("w", encoding="utf-8") as handle:
        handle.write("COMMAND:\n")
        handle.write(" ".join(shlex.quote(x) for x in cmd) + "\n\n")
        handle.flush()
        completed = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"{step} failed with returncode={completed.returncode}. See {log_path}")
    if not expected_output.exists():
        raise FileNotFoundError(f"{step} completed but expected output was not found: {expected_output}")
    print(f"Finished {step}: {expected_output}")


def main() -> None:
    args = parse_args()
    ensure_parent_dirs()

    expression_csv = Path(args.expression_csv)
    for resource in [Path(args.tf_list), Path(args.ranking_db), Path(args.motif_annotation)]:
        if not resource.exists():
            raise FileNotFoundError(f"Missing pySCENIC resource: {resource}")

    commands = build_commands(args)
    command_report = write_command_report(commands)
    print(f"Wrote {command_report}")

    if args.dry_run:
        if not expression_csv.exists():
            print(f"Dry-run note: expression matrix does not exist yet: {expression_csv}")
        print("Dry-run only; commands were generated but not executed.")
        return

    if not expression_csv.exists():
        raise FileNotFoundError(
            f"Missing expression matrix: {expression_csv}\n"
            "Run export_round2_for_pyscenic.py first."
        )

    for step, cmd, expected in commands:
        run_step(step, cmd, expected, args.force)

    summary_path = LOG_DIR / "pipeline_summary.txt"
    summary_path.write_text(
        "\n".join(
            [
                "round2 pySCENIC pipeline completed",
                f"grn_adj: {OUTPUT_DIR / 'grn_adj.tsv'}",
                f"regulons: {OUTPUT_DIR / 'regulons.csv'}",
                f"auc_mtx: {OUTPUT_DIR / 'auc_mtx.csv'}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
