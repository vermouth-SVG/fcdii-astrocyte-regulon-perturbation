#!/usr/bin/env python3
from __future__ import annotations

import argparse
import ast
import json
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


ROUND1_KO_TFS = ["BHLHE40", "NFE2L2", "SOX2", "THRB"]


def project_root_from_file(script_file: str | Path) -> Path:
    return Path(script_file).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    root = project_root_from_file(__file__)
    parser = argparse.ArgumentParser(
        description="Prepare first-round CellOracle analysis objects using the current pySCENIC shortlist."
    )
    parser.add_argument(
        "--h5ad",
        default=str(root / "final_exports" / "astrocyte_pilot_rna_with_pyscenic_auc.h5ad"),
        help="Merged h5ad with pySCENIC AUC already added.",
    )
    parser.add_argument(
        "--candidate-dir",
        default=str(root / "analysis_outputs" / "celloracle_candidates"),
        help="Directory containing shortlisted TF results.",
    )
    parser.add_argument(
        "--regulons-csv",
        default=str(root / "output" / "regulons.csv"),
        help="pySCENIC regulons.csv file.",
    )
    parser.add_argument(
        "--grn-adj",
        default=str(root / "output" / "grn_adj.tsv"),
        help="pySCENIC grn adjacency file.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(root / "celloracle_run"),
        help="Output directory root for CellOracle round 1.",
    )
    parser.add_argument(
        "--cluster-column",
        default="group",
        help="obs column used as the CellOracle cluster unit.",
    )
    parser.add_argument(
        "--embedding-key",
        default="X_umap",
        help="Embedding key to generate or reuse.",
    )
    parser.add_argument(
        "--hvg-n-top",
        type=int,
        default=2500,
        help="Number of HVGs to keep before union with base-GRN genes.",
    )
    parser.add_argument(
        "--max-targets-per-tf",
        type=int,
        default=300,
        help="Maximum number of target genes retained per TF for the first-round base GRN.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=10.0,
        help="Regularization strength for CellOracle get_links.",
    )
    parser.add_argument(
        "--bagging-number",
        type=int,
        default=20,
        help="Bagging number for GRN inference.",
    )
    parser.add_argument(
        "--edge-threshold",
        type=int,
        default=2000,
        help="Top edge number retained in filtered links per cluster.",
    )
    parser.add_argument(
        "--edge-pvalue",
        type=float,
        default=0.001,
        help="P-value threshold used in link filtering.",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=-1,
        help="CPU cores used in CellOracle GRN inference.",
    )
    parser.add_argument(
        "--ko-tfs",
        default=",".join(ROUND1_KO_TFS),
        help="Comma-separated TFs for round-1 perturbation simulation.",
    )
    return parser.parse_args()


def load_shortlist_tfs(candidate_dir: Path) -> list[str]:
    shortlist_path = candidate_dir / "celloracle_tf_shortlist.csv"
    if not shortlist_path.exists():
        raise FileNotFoundError(f"Missing shortlist file: {shortlist_path}")
    df = pd.read_csv(shortlist_path)
    if "tf" not in df.columns:
        raise ValueError(f"'tf' column not found in {shortlist_path}")
    return df["tf"].astype(str).dropna().tolist()


def _parse_target_gene_field(raw_value: str) -> list[tuple[str, float]]:
    parsed = ast.literal_eval(raw_value)
    output: list[tuple[str, float]] = []
    if isinstance(parsed, (list, tuple)):
        for item in parsed:
            if isinstance(item, (list, tuple)) and item:
                gene = str(item[0])
                weight = float(item[1]) if len(item) > 1 else np.nan
                output.append((gene, weight))
    return output


def load_edges_from_regulons(regulons_csv: Path, tf_list: list[str], max_targets_per_tf: int) -> pd.DataFrame:
    df = pd.read_csv(regulons_csv, skiprows=1)
    if df.shape[1] < 9:
        raise ValueError(f"Unexpected regulons.csv format: {regulons_csv}")
    df = df.iloc[:, :10].copy()
    df.columns = [
        "TF",
        "MotifID",
        "AUC",
        "NES",
        "MotifSimilarityQvalue",
        "OrthologousIdentity",
        "Annotation",
        "Context",
        "TargetGenes",
        "RankAtMax",
    ]
    df = df.loc[df["TF"].astype(str).isin(tf_list)].copy()

    rows: list[dict[str, object]] = []
    for _, row in df.iterrows():
        tf = str(row["TF"])
        try:
            targets = _parse_target_gene_field(str(row["TargetGenes"]))
        except Exception:
            targets = []
        for gene, weight in targets:
            rows.append(
                {
                    "tf": tf,
                    "target": gene,
                    "weight": float(weight) if np.isfinite(weight) else np.nan,
                    "source": "regulons_csv",
                }
            )

    if not rows:
        return pd.DataFrame(columns=["tf", "target", "weight", "source"])

    edges = pd.DataFrame(rows)
    edges = edges.sort_values(by=["tf", "weight"], ascending=[True, False])
    edges = edges.drop_duplicates(subset=["tf", "target"], keep="first")
    edges = edges.groupby("tf", group_keys=False).head(max_targets_per_tf).reset_index(drop=True)
    return edges


def load_edges_from_grn_adj(
    grn_adj_path: Path,
    tf_list: list[str],
    existing_edges: pd.DataFrame,
    max_targets_per_tf: int,
) -> pd.DataFrame:
    existing_targets: dict[str, set[str]] = defaultdict(set)
    if not existing_edges.empty:
        for _, row in existing_edges.iterrows():
            existing_targets[str(row["tf"])].add(str(row["target"]))

    collected: list[pd.DataFrame] = []
    usecols = ["TF", "target", "importance"]
    for chunk in pd.read_csv(grn_adj_path, sep="\t", usecols=usecols, chunksize=200000):
        sub = chunk.loc[chunk["TF"].astype(str).isin(tf_list)].copy()
        if sub.empty:
            continue
        sub.columns = ["tf", "target", "weight"]
        collected.append(sub)

    if not collected:
        return pd.DataFrame(columns=["tf", "target", "weight", "source"])

    fallback = pd.concat(collected, axis=0, ignore_index=True)
    fallback["tf"] = fallback["tf"].astype(str)
    fallback["target"] = fallback["target"].astype(str)
    fallback = fallback.sort_values(by=["tf", "weight"], ascending=[True, False])

    rows: list[dict[str, object]] = []
    per_tf_count = {tf: len(existing_targets.get(tf, set())) for tf in tf_list}
    for _, row in fallback.iterrows():
        tf = str(row["tf"])
        target = str(row["target"])
        if target in existing_targets.get(tf, set()):
            continue
        if per_tf_count[tf] >= max_targets_per_tf:
            continue
        rows.append(
            {
                "tf": tf,
                "target": target,
                "weight": float(row["weight"]),
                "source": "grn_adj_fallback",
            }
        )
        existing_targets[tf].add(target)
        per_tf_count[tf] += 1

    return pd.DataFrame(rows)


def build_tg_to_tf_dict(edges: pd.DataFrame) -> dict[str, list[str]]:
    tg_to_tf: dict[str, set[str]] = defaultdict(set)
    for _, row in edges.iterrows():
        tg_to_tf[str(row["target"])].add(str(row["tf"]))
    return {target: sorted(list(tfs)) for target, tfs in tg_to_tf.items()}


def choose_pca_dims(adata) -> int:
    n_obs = int(adata.n_obs)
    n_vars = int(adata.n_vars)
    upper = max(2, min(50, n_obs - 1, n_vars - 1))
    return max(10, min(upper, 30))


def choose_knn_k(n_cells: int) -> int:
    return max(20, min(80, int(round(n_cells * 0.025))))


def main() -> None:
    args = parse_args()
    root = project_root_from_file(__file__)

    output_dir = Path(args.output_dir).resolve()
    checks_dir = output_dir / "checks"
    prepared_dir = output_dir / "prepared_data"
    network_dir = output_dir / "network"
    for path in [checks_dir, prepared_dir, network_dir]:
        path.mkdir(parents=True, exist_ok=True)

    print("[1/7] Importing CellOracle dependencies")
    try:
        import anndata as ad
        import celloracle as co
        import scanpy as sc
        from scipy import sparse
    except Exception as exc:
        raise RuntimeError(
            "CellOracle preparation requires an environment with anndata, scanpy, and celloracle. "
            "Please run this script inside the CellOracle Docker image."
        ) from exc

    warnings.filterwarnings("ignore", category=FutureWarning)

    shortlist_tfs = load_shortlist_tfs(Path(args.candidate_dir))
    ko_tfs = [item.strip() for item in args.ko_tfs.split(",") if item.strip()]
    missing_ko = [tf for tf in ko_tfs if tf not in shortlist_tfs]
    if missing_ko:
        raise ValueError(f"Round-1 perturbation TFs are not all in the shortlist: {missing_ko}")

    print("[2/7] Loading base-GRN edges from current pySCENIC results")
    regulon_edges = load_edges_from_regulons(Path(args.regulons_csv), shortlist_tfs, args.max_targets_per_tf)
    fallback_edges = load_edges_from_grn_adj(Path(args.grn_adj), shortlist_tfs, regulon_edges, args.max_targets_per_tf)
    edges = pd.concat([regulon_edges, fallback_edges], axis=0, ignore_index=True)
    if edges.empty:
        raise ValueError("No TF-target edges were recovered from regulons.csv or grn_adj.tsv.")
    edges["tf"] = edges["tf"].astype(str)
    edges["target"] = edges["target"].astype(str)
    edges = edges.sort_values(by=["tf", "weight"], ascending=[True, False]).reset_index(drop=True)

    print("[3/7] Reading h5ad and building minimal CellOracle-ready embedding")
    adata = sc.read_h5ad(args.h5ad)
    if adata.raw is not None:
        adata = adata.raw.to_adata()
    else:
        adata = adata.copy()

    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    adata.var_names_make_unique()

    if args.cluster_column not in adata.obs.columns:
        raise KeyError(f"Cluster column not found in obs: {args.cluster_column}")
    adata.obs["celloracle_cluster"] = adata.obs[args.cluster_column].astype(str)
    if "group" in adata.obs.columns:
        adata.obs["celloracle_group"] = adata.obs["group"].astype(str)

    counts = adata.X.copy()
    if sparse.issparse(counts):
        counts_data = counts.data
        if counts_data.size and not np.allclose(counts_data, np.round(counts_data)):
            raise ValueError("Current counts matrix is not integer-like.")
    else:
        if not np.allclose(np.asarray(counts), np.round(np.asarray(counts))):
            raise ValueError("Current counts matrix is not integer-like.")
    adata.layers["raw_count"] = counts.copy()
    adata.layers["counts"] = counts.copy()

    sc.pp.filter_genes(adata, min_counts=1)

    base_genes = set(shortlist_tfs) | set(edges["target"].astype(str).tolist())
    base_genes = {gene for gene in base_genes if gene in adata.var_names}

    adata_hvg = adata.copy()
    try:
        sc.pp.highly_variable_genes(
            adata_hvg,
            layer="counts",
            flavor="seurat_v3",
            n_top_genes=args.hvg_n_top,
            subset=True,
        )
        hvg_genes = set(adata_hvg.var_names)
    except Exception as exc:
        print(f"[3/7] seurat_v3 HVG failed, fallback to seurat: {exc}")
        adata_hvg = adata.copy()
        sc.pp.normalize_total(adata_hvg, target_sum=1e4)
        sc.pp.log1p(adata_hvg)
        sc.pp.highly_variable_genes(
            adata_hvg,
            flavor="seurat",
            n_top_genes=args.hvg_n_top,
            subset=True,
        )
        hvg_genes = set(adata_hvg.var_names)

    genes_to_keep = hvg_genes | base_genes
    keep_order = [gene for gene in adata.var_names if gene in genes_to_keep]
    adata = adata[:, keep_order].copy()
    adata.layers["raw_count"] = adata.layers["raw_count"].copy()
    adata.layers["counts"] = adata.layers["counts"].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.scale(adata, max_value=10)
    n_pcs = choose_pca_dims(adata)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=n_pcs)
    sc.pp.neighbors(adata, n_neighbors=min(30, max(10, adata.n_obs - 1)), n_pcs=min(n_pcs, 30))
    sc.tl.umap(adata)
    try:
        sc.tl.leiden(adata, key_added="celloracle_leiden", resolution=0.6)
    except Exception:
        pass

    prepared_h5ad_path = prepared_dir / "astrocyte_pilot_celloracle_round1_input.h5ad"
    adata.write_h5ad(prepared_h5ad_path)

    # Filter edges to genes retained in the final object.
    edges = edges.loc[edges["target"].isin(adata.var_names)].copy().reset_index(drop=True)
    tg_to_tf = build_tg_to_tf_dict(edges)
    if not tg_to_tf:
        raise ValueError("No base-GRN target genes remain after gene subsetting.")

    shortlist_df = pd.DataFrame({"tf": shortlist_tfs})
    shortlist_df.to_csv(prepared_dir / "celloracle_round1_shortlist_tfs.csv", index=False)
    pd.DataFrame({"tf": ko_tfs}).to_csv(prepared_dir / "celloracle_round1_ko_tfs.csv", index=False)
    edges.to_csv(prepared_dir / "celloracle_round1_base_grn_edges.csv", index=False)
    with open(prepared_dir / "celloracle_round1_tg_to_tf.json", "w", encoding="utf-8") as handle:
        json.dump(tg_to_tf, handle, indent=2, ensure_ascii=False)
    adata.obs.to_csv(prepared_dir / "celloracle_round1_metadata.csv")

    print("[4/7] Creating Oracle object")
    adata_for_oracle = adata.copy()
    adata_for_oracle.X = adata_for_oracle.layers["raw_count"].copy()

    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=adata_for_oracle,
        cluster_column_name="celloracle_cluster",
        embedding_name=args.embedding_key,
    )
    oracle.import_TF_data(TFdict=tg_to_tf)
    oracle.perform_PCA()

    k = choose_knn_k(int(adata.n_obs))
    b_sight = min(int(max(k * 8, 40)), max(adata.n_obs - 1, 1))
    b_maxl = min(int(max(k * 4, 20)), max(adata.n_obs - 1, 1))
    try:
        oracle.knn_imputation(
            n_pca_dims=n_pcs,
            k=k,
            balanced=True,
            b_sight=b_sight,
            b_maxl=b_maxl,
            n_jobs=args.n_jobs,
        )
    except TypeError:
        oracle.knn_imputation(
            n_pca_dims=n_pcs,
            k=k,
            balanced=True,
            b_sight=b_sight,
            b_maxl=b_maxl,
        )

    oracle_path = network_dir / "oracle_round1_base.celloracle.oracle"
    oracle.to_hdf5(str(oracle_path))

    print("[5/7] Inferring cluster-specific GRNs")
    links = oracle.get_links(
        cluster_name_for_GRN_unit="celloracle_cluster",
        alpha=args.alpha,
        bagging_number=args.bagging_number,
        verbose_level=1,
        model_method="bagging_ridge",
        n_jobs=args.n_jobs,
    )
    links.filter_links(
        p=0.001,
        weight="coef_abs",
        threshold_number=10000,
    )
    links.get_network_score()
    links_path = network_dir / "links_round1_raw.celloracle.links"
    links.to_hdf5(str(links_path))

    print("[6/7] Filtering links and writing summaries")
    links.filter_links(p=args.edge_pvalue, weight="coef_abs", threshold_number=args.edge_threshold)
    filtered_links_path = network_dir / "links_round1_filtered.celloracle.links"
    links.to_hdf5(str(filtered_links_path))

    if hasattr(links, "merged_score") and links.merged_score is not None:
        links.merged_score.to_csv(network_dir / "links_round1_merged_score.csv", index=False)

    filtered_rows: list[pd.DataFrame] = []
    if hasattr(links, "filtered_links") and isinstance(links.filtered_links, dict):
        for cluster_name, df_cluster in links.filtered_links.items():
            if df_cluster is None or len(df_cluster) == 0:
                continue
            tmp = df_cluster.copy()
            tmp["cluster"] = str(cluster_name)
            filtered_rows.append(tmp)
    if filtered_rows:
        pd.concat(filtered_rows, axis=0, ignore_index=True).to_csv(
            network_dir / "links_round1_filtered_long.csv",
            index=False,
        )

    prep_config = {
        "project_root": str(root),
        "prepared_h5ad": str(prepared_h5ad_path),
        "oracle_path": str(oracle_path),
        "links_raw_path": str(links_path),
        "links_filtered_path": str(filtered_links_path),
        "base_grn_edges_csv": str(prepared_dir / "celloracle_round1_base_grn_edges.csv"),
        "base_grn_tg_to_tf_json": str(prepared_dir / "celloracle_round1_tg_to_tf.json"),
        "shortlist_tfs_csv": str(prepared_dir / "celloracle_round1_shortlist_tfs.csv"),
        "ko_tfs_csv": str(prepared_dir / "celloracle_round1_ko_tfs.csv"),
        "cluster_column": "celloracle_cluster",
        "embedding_key": args.embedding_key,
        "alpha": args.alpha,
        "bagging_number": args.bagging_number,
        "edge_threshold": args.edge_threshold,
        "edge_pvalue": args.edge_pvalue,
        "n_pcs": n_pcs,
        "knn_k": k,
    }
    with open(output_dir / "prep_config_round1.json", "w", encoding="utf-8") as handle:
        json.dump(prep_config, handle, indent=2, ensure_ascii=False)

    summary_lines = [
        f"prepared_h5ad={prepared_h5ad_path}",
        f"oracle_path={oracle_path}",
        f"links_filtered_path={filtered_links_path}",
        f"cluster_column=celloracle_cluster",
        f"embedding_key={args.embedding_key}",
        f"n_cells={adata.n_obs}",
        f"n_genes_after_subset={adata.n_vars}",
        f"shortlist_tfs={len(shortlist_tfs)}",
        f"ko_tfs={','.join(ko_tfs)}",
        f"base_grn_edges={edges.shape[0]}",
        f"base_grn_targets={len(tg_to_tf)}",
        f"n_pcs={n_pcs}",
        f"knn_k={k}",
    ]
    with open(output_dir / "prep_summary_round1.txt", "w", encoding="utf-8") as handle:
        handle.write("\n".join(summary_lines) + "\n")

    print("[7/7] Done")
    print(f"  prepared_h5ad={prepared_h5ad_path}")
    print(f"  oracle_path={oracle_path}")
    print(f"  links_filtered_path={filtered_links_path}")
    print(f"  outputs={output_dir}")


if __name__ == "__main__":
    main()
