#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXPRESSION = ROOT / "external_validation_round2" / "pyscenic_recalc" / "input" / "expression_for_pyscenic_4tf_targets.csv"
DEFAULT_GMT = ROOT / "external_validation_round2" / "pyscenic_recalc" / "output" / "discovery_regulons_4tf.gmt"
DEFAULT_OUTPUT = ROOT / "external_validation_round2" / "pyscenic_recalc" / "output" / "auc_mtx_4tf_from_discovery.csv"
DEFAULT_LOG = ROOT / "external_validation_round2" / "pyscenic_recalc" / "logs" / "03_aucell_4tf_from_discovery.log"
DEFAULT_REPORT = ROOT / "external_validation_round2" / "pyscenic_recalc" / "logs" / "auc_projection_report.txt"
DEFAULT_IMAGE = "aertslab/pyscenic:0.12.1"


def to_docker_mount_path(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    tail = resolved.as_posix().split(":/", 1)[1]
    return f"{drive}:/{tail}"


def to_container_path(path: Path) -> str:
    return "/work/" + path.resolve().relative_to(ROOT.resolve()).as_posix()


def read_csv_shape(path: Path) -> tuple[int, int]:
    rows = 0
    cols = 0
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, [])
        cols = max(len(header) - 1, 0)
        for _ in reader:
            rows += 1
    return rows, cols


def main() -> None:
    parser = argparse.ArgumentParser(description="Run AUCell projection for 4 discovery regulons on minimal round2 matrix.")
    parser.add_argument("--expression", default=str(DEFAULT_EXPRESSION))
    parser.add_argument("--gmt", default=str(DEFAULT_GMT))
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--log", default=str(DEFAULT_LOG))
    parser.add_argument("--report", default=str(DEFAULT_REPORT))
    parser.add_argument("--image", default=DEFAULT_IMAGE)
    parser.add_argument("--num-workers", type=int, default=1)
    args = parser.parse_args()

    expression = Path(args.expression)
    gmt = Path(args.gmt)
    output = Path(args.output)
    log = Path(args.log)
    report = Path(args.report)

    output.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    report.parent.mkdir(parents=True, exist_ok=True)

    if not expression.exists():
        raise FileNotFoundError(f"Missing minimal expression matrix: {expression}")
    if not gmt.exists():
        raise FileNotFoundError(f"Missing discovery GMT: {gmt}")
    if shutil.which("docker") is None:
        raise RuntimeError("Docker CLI not found in PATH. Run this script from the Windows host with Docker Desktop available.")

    if output.exists():
        output.unlink()

    cmd = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{to_docker_mount_path(ROOT)}:/work",
        "-w",
        "/work",
        args.image,
        "pyscenic",
        "aucell",
        to_container_path(expression),
        to_container_path(gmt),
        "-o",
        to_container_path(output),
        "--num_workers",
        str(args.num_workers),
    ]

    with log.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("Command:\n")
        handle.write(" ".join(cmd) + "\n\n")
        process = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, text=True, check=False)

    success = output.exists() and output.stat().st_size > 0 and process.returncode == 0
    output_rows = 0
    output_regulons = 0
    if output.exists() and output.stat().st_size > 0:
        output_rows, output_regulons = read_csv_shape(output)

    report_lines = [
        f"Expression input: {expression}",
        f"GMT input: {gmt}",
        f"Docker image: {args.image}",
        f"Workers: {args.num_workers}",
        f"Return code: {process.returncode}",
        f"Output exists: {output.exists()}",
        f"Output non-empty: {output.exists() and output.stat().st_size > 0 if output.exists() else False}",
        f"Output path: {output}",
        f"Output size bytes: {output.stat().st_size if output.exists() else 0}",
        f"Output cells: {output_rows}",
        f"Output regulons: {output_regulons}",
        f"Log path: {log}",
        f"Projection success: {success}",
    ]
    report.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    if not success:
        raise SystemExit(
            f"AUCell projection failed or output is empty.\n"
            f"See log: {log}\n"
            f"See report: {report}"
        )

    print(f"Wrote {output}")
    print(f"Wrote {log}")
    print(f"Wrote {report}")
    print(f"Cells={output_rows}, regulons={output_regulons}")


if __name__ == "__main__":
    main()
