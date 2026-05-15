#!/usr/bin/env python3
from __future__ import annotations

import gzip
import json
import re
import tarfile
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "external_validation_round2"
REPORT_TXT = OUT_DIR / "candidate_assessment_report.txt"
REPORT_JSON = OUT_DIR / "candidate_assessment_report.json"


@dataclass
class CandidateResult:
    candidate_name: str
    source_path: str
    source_type: str
    available: bool
    has_group: bool
    has_sample_or_donor: bool
    has_celltype: bool
    can_generate_cluster: bool
    recommended: bool
    summary: str
    extra: dict[str, object]


def safe_bool(value: bool) -> str:
    return "yes" if value else "no"


def fetch_geo_status(accession: str) -> dict[str, object]:
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
    result: dict[str, object] = {"accession": accession, "url": url, "reachable": False}
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            html = response.read().decode("utf-8", errors="ignore")
    except urllib.error.URLError as exc:
        result["error"] = str(exc)
        return result

    result["reachable"] = True
    private_match = re.search(
        r'Accession "(GSE\d+)" is currently private and is scheduled to be released on ([A-Za-z]{3} \d{1,2}, \d{4})',
        html,
    )
    if private_match:
        result["status"] = "private"
        result["release_date"] = private_match.group(2)
        return result

    title_match = re.search(r"!Series_title = ([^\n<]+)", html)
    if title_match:
        result["status"] = "public"
        result["series_title"] = title_match.group(1).strip()
    elif f"Series&nbsp;{accession}" in html:
        result["status"] = "public"
    else:
        result["status"] = "unknown"
    return result


def inspect_existing_h5ad(path: Path, candidate_name: str) -> CandidateResult:
    if not path.exists():
        return CandidateResult(
            candidate_name=candidate_name,
            source_path=str(path),
            source_type="h5ad",
            available=False,
            has_group=False,
            has_sample_or_donor=False,
            has_celltype=False,
            can_generate_cluster=False,
            recommended=False,
            summary="对象不存在。",
            extra={},
        )

    adata = ad.read_h5ad(path, backed="r")
    obs_cols = list(adata.obs.columns)
    extra: dict[str, object] = {
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "obs_columns": obs_cols,
        "obsm_keys": list(adata.obsm.keys()),
        "uns_keys": list(adata.uns.keys()),
    }

    group_col = next((c for c in ["group", "disease", "condition", "diagnosis", "status"] if c in obs_cols), None)
    sample_col = next((c for c in ["sample", "sample_id", "donor", "donor_id", "patient", "patient_id"] if c in obs_cols), None)
    celltype_col = next((c for c in ["celltype", "cell_type", "broad_celltype", "cluster", "leiden"] if c in obs_cols), None)

    has_group = False
    if group_col is not None:
        groups = adata.obs[group_col].astype("string").dropna().unique().tolist()
        extra["group_column"] = group_col
        extra["group_levels"] = groups
        has_group = len(groups) >= 2
    has_sample = sample_col is not None
    if sample_col is not None:
        extra["sample_or_donor_column"] = sample_col
        extra["sample_or_donor_nunique"] = int(adata.obs[sample_col].astype("string").nunique())
    has_celltype = False
    if celltype_col is not None:
        extra["celltype_or_cluster_column"] = celltype_col
        extra["celltype_or_cluster_nunique"] = int(adata.obs[celltype_col].astype("string").nunique())
        has_celltype = extra["celltype_or_cluster_nunique"] >= 2

    summary = (
        f"{candidate_name}: 已存在标准 h5ad；"
        f"group={group_col or 'none'}，sample/donor={sample_col or 'none'}，celltype/cluster={celltype_col or 'none'}。"
    )
    recommended = has_group and has_sample and has_celltype
    adata.file.close()
    return CandidateResult(
        candidate_name=candidate_name,
        source_path=str(path),
        source_type="h5ad",
        available=True,
        has_group=has_group,
        has_sample_or_donor=has_sample,
        has_celltype=has_celltype,
        can_generate_cluster=False,
        recommended=recommended,
        summary=summary,
        extra=extra,
    )


def parse_soft_metadata(soft_path: Path) -> pd.DataFrame:
    records: list[dict[str, str]] = []
    current: dict[str, object] | None = None
    with gzip.open(soft_path, "rt", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    records.append(_finalize_soft(current))
                current = {"sample_id": line.split("=", 1)[1].strip(), "characteristics": [], "title": ""}
                continue
            if current is None:
                continue
            if line.startswith("!Sample_title = "):
                current["title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1 = "):
                current["characteristics"].append(line.split("=", 1)[1].strip())
    if current is not None:
        records.append(_finalize_soft(current))
    return pd.DataFrame(records)


def _finalize_soft(record: dict[str, object]) -> dict[str, str]:
    char_map: dict[str, str] = {}
    for item in record.get("characteristics", []):
        if ":" in item:
            key, value = item.split(":", 1)
            char_map[key.strip().lower()] = value.strip()
    diagnosis = char_map.get("diagnosis", "")
    group = "unknown"
    disease = diagnosis
    lowered = diagnosis.lower()
    if "control" in lowered or "non-epileptic" in lowered or "normal" in lowered:
        group = "internal_control"
        disease = "control"
    elif "tle" in lowered or "epilep" in lowered:
        group = "lesion"
        disease = "TLE"
    sample_name = str(record.get("title", "")).replace("single-cell_", "").strip()
    if not sample_name:
        sample_name = str(record.get("sample_id", ""))
    return {
        "sample_id": str(record.get("sample_id", "")),
        "title": str(record.get("title", "")),
        "diagnosis": diagnosis,
        "group": group,
        "disease": disease,
        "sample": sample_name,
        "donor": sample_name,
    }


def inspect_raw_series(series_dir: Path, accession: str, candidate_name: str) -> CandidateResult:
    soft_path = series_dir / f"{accession}_family.soft.gz"
    extract_dir = series_dir / "extracted"
    raw_tar = series_dir / f"{accession}_RAW.tar"
    if not soft_path.exists():
        return CandidateResult(
            candidate_name=candidate_name,
            source_path=str(series_dir),
            source_type="geo_raw",
            available=False,
            has_group=False,
            has_sample_or_donor=False,
            has_celltype=False,
            can_generate_cluster=False,
            recommended=False,
            summary=f"{candidate_name}: 缺少 {soft_path.name}。",
            extra={},
        )

    meta = parse_soft_metadata(soft_path)
    has_group = bool(meta["group"].nunique() >= 2) if not meta.empty and "group" in meta.columns else False
    has_sample = bool("sample" in meta.columns and meta["sample"].nunique() >= 2)
    has_celltype = False

    matrix_files = sorted(extract_dir.glob("*matrix.mtx.gz")) if extract_dir.exists() else []
    feature_files = sorted(extract_dir.glob("*features.tsv.gz")) if extract_dir.exists() else []
    barcode_files = sorted(extract_dir.glob("*barcodes.tsv.gz")) if extract_dir.exists() else []
    can_generate_cluster = bool(matrix_files and feature_files and barcode_files)

    sample_groups = {}
    if not meta.empty:
        sample_groups = meta[["sample", "group", "disease"]].to_dict(orient="records")

    summary = (
        f"{candidate_name}: "
        f"样本数={int(meta.shape[0]) if meta is not None else 0}；"
        f"group={safe_bool(has_group)}；sample/donor={safe_bool(has_sample)}；"
        f"原始 celltype 注释={safe_bool(has_celltype)}；"
        f"可自动生成 cluster={safe_bool(can_generate_cluster)}。"
    )
    recommended = has_group and has_sample and can_generate_cluster
    extra = {
        "soft_exists": soft_path.exists(),
        "raw_tar_exists": raw_tar.exists(),
        "extracted_matrix_files": [p.name for p in matrix_files],
        "sample_metadata": sample_groups,
        "group_levels": sorted(meta["group"].dropna().unique().tolist()) if "group" in meta.columns else [],
    }
    return CandidateResult(
        candidate_name=candidate_name,
        source_path=str(series_dir),
        source_type="geo_raw",
        available=True,
        has_group=has_group,
        has_sample_or_donor=has_sample,
        has_celltype=has_celltype,
        can_generate_cluster=can_generate_cluster,
        recommended=recommended,
        summary=summary,
        extra=extra,
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results: list[CandidateResult] = []

    geo275302 = fetch_geo_status("GSE275302")
    if geo275302.get("status") == "private":
        results.append(
            CandidateResult(
                candidate_name="GSE275302",
                source_path=str(geo275302.get("url", "")),
                source_type="geo_remote",
                available=False,
                has_group=False,
                has_sample_or_donor=False,
                has_celltype=False,
                can_generate_cluster=False,
                recommended=False,
                summary=(
                    f"GSE275302 当前不可用于正式外部验证：该 GEO accession 处于 private 状态，"
                    f"计划公开日期为 {geo275302.get('release_date', 'unknown')}。"
                ),
                extra=geo275302,
            )
        )
    else:
        results.append(
            CandidateResult(
                candidate_name="GSE275302",
                source_path=str(geo275302.get("url", "")),
                source_type="geo_remote",
                available=bool(geo275302.get("reachable")),
                has_group=False,
                has_sample_or_donor=False,
                has_celltype=False,
                can_generate_cluster=False,
                recommended=False,
                summary=f"GSE275302 远程状态：{geo275302.get('status', 'unknown')}",
                extra=geo275302,
            )
        )

    results.append(
        inspect_existing_h5ad(
            ROOT / "external_validation" / "input" / "external_validation.h5ad",
            "GSE140393-derived external_validation.h5ad",
        )
    )
    results.append(
        inspect_raw_series(
            ROOT / "external_validation_round2" / "source_gse190452",
            accession="GSE190452",
            candidate_name="GSE190452",
        )
    )

    recommended = [r for r in results if r.recommended]
    preferred = next((r for r in recommended if r.candidate_name == "GSE190452"), recommended[0] if recommended else None)

    lines: list[str] = []
    lines.append("第二个正式外部验证队列候选评估报告")
    lines.append("")
    for item in results:
        lines.append(f"[{item.candidate_name}]")
        lines.append(f"source_path: {item.source_path}")
        lines.append(f"source_type: {item.source_type}")
        lines.append(f"available: {safe_bool(item.available)}")
        lines.append(f"has_group: {safe_bool(item.has_group)}")
        lines.append(f"has_sample_or_donor: {safe_bool(item.has_sample_or_donor)}")
        lines.append(f"has_celltype: {safe_bool(item.has_celltype)}")
        lines.append(f"can_generate_cluster: {safe_bool(item.can_generate_cluster)}")
        lines.append(f"recommended: {safe_bool(item.recommended)}")
        lines.append(f"summary: {item.summary}")
        if item.extra:
            lines.append("extra:")
            lines.append(json.dumps(item.extra, ensure_ascii=False, indent=2))
        lines.append("")

    lines.append("[结论]")
    if preferred is not None:
        lines.append(
            f"推荐将 {preferred.candidate_name} 作为第二个正式外部验证队列。"
        )
        if preferred.candidate_name == "GSE190452":
            lines.append(
                "推荐理由：该队列在当前工作区已具备可直接读取的 filtered 10x matrix，"
                "同时包含 Normal 与 TLE 双组样本，样本层面可保留 sample/donor 信息；"
                "虽然原始数据没有现成 celltype 注释，但可以通过最小 Scanpy 流程自动生成 cluster，用于后续 group/celltype 双层验证。"
            )
    else:
        lines.append("当前未发现满足正式外部验证最低要求的候选队列。")

    REPORT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    REPORT_JSON.write_text(
        json.dumps(
            {
                "results": [item.__dict__ for item in results],
                "recommended_candidate": preferred.candidate_name if preferred is not None else None,
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {REPORT_TXT}")
    print(f"Wrote {REPORT_JSON}")
    if preferred is not None:
        print(f"Recommended candidate: {preferred.candidate_name}")


if __name__ == "__main__":
    main()
