#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_KO_TFS = ["BHLHE40", "NFE2L2", "SOX2", "THRB"]


def project_root_from_file(script_file: str | Path) -> Path:
    return Path(script_file).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    root = project_root_from_file(__file__)
    parser = argparse.ArgumentParser(
        description="Run first-round CellOracle in silico perturbation for four prioritized TFs."
    )
    parser.add_argument(
        "--config",
        default=str(root / "celloracle_run" / "prep_config_round1.json"),
        help="Config JSON produced by prepare_celloracle_first_round.py.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(root / "celloracle_run" / "ko_round1"),
        help="Output directory for round-1 perturbation results.",
    )
    parser.add_argument(
        "--candidate-dir",
        default=str(root / "analysis_outputs" / "celloracle_candidates"),
        help="Candidate TF screening output directory.",
    )
    parser.add_argument(
        "--tf-list",
        default=",".join(DEFAULT_KO_TFS),
        help="Comma-separated TF list to perturb in round 1.",
    )
    parser.add_argument(
        "--simulation-alpha",
        type=float,
        default=10.0,
        help="Alpha used in fit_GRN_for_simulation.",
    )
    parser.add_argument(
        "--n-propagation",
        type=int,
        default=3,
        help="Propagation steps for simulate_shift.",
    )
    parser.add_argument(
        "--transition-neighbors",
        type=int,
        default=200,
        help="Neighbors used in estimate_transition_prob and p-mass calculation.",
    )
    parser.add_argument(
        "--sigma-corr",
        type=float,
        default=0.05,
        help="Sigma for calculate_embedding_shift.",
    )
    parser.add_argument(
        "--n-grid",
        type=int,
        default=40,
        help="Grid number for vector-field visualization.",
    )
    parser.add_argument(
        "--min-mass",
        type=float,
        default=0.01,
        help="Mass threshold for vector-field visualization.",
    )
    parser.add_argument(
        "--quiver-scale",
        type=float,
        default=25.0,
        help="Scale used in plot_quiver.",
    )
    parser.add_argument(
        "--grid-scale",
        type=float,
        default=0.5,
        help="Scale used in plot_simulation_flow_on_grid.",
    )
    parser.add_argument(
        "--markov-steps",
        type=int,
        default=500,
        help="Steps for Markov simulation.",
    )
    parser.add_argument(
        "--markov-duplication",
        type=int,
        default=5,
        help="Duplication number for Markov simulation.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=123,
        help="Random seed for Markov simulation.",
    )
    parser.add_argument(
        "--force-refit",
        action="store_true",
        help="Force rebuilding the simulation-ready Oracle object.",
    )
    return parser.parse_args()


def safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in str(value))


def load_candidate_metrics(candidate_dir: Path) -> pd.DataFrame:
    path = candidate_dir / "celloracle_candidate_tf_metrics.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    df = pd.read_csv(path)
    if "tf" not in df.columns:
        raise ValueError(f"'tf' column not found in {path}")
    return df


def ensure_simulation_ready_oracle(co, config: dict[str, object], output_dir: Path, alpha: float, force_refit: bool):
    sim_ready_dir = output_dir / "simulation_ready"
    sim_ready_dir.mkdir(parents=True, exist_ok=True)
    sim_ready_path = sim_ready_dir / "oracle_round1_simulation_ready.celloracle.oracle"
    if sim_ready_path.exists() and not force_refit:
        return sim_ready_path

    oracle = co.load_hdf5(str(config["oracle_path"]))
    links = co.load_hdf5(str(config["links_filtered_path"]))
    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    oracle.fit_GRN_for_simulation(alpha=alpha, use_cluster_specific_TFdict=True, verbose_level=1)
    oracle.to_hdf5(str(sim_ready_path))
    return sim_ready_path


def extract_embedding_and_shift(oracle, embedding_key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    embedding = np.asarray(oracle.adata.obsm[embedding_key], dtype=float)
    delta = None
    delta_random = None
    candidate_attrs = ["delta_embedding", "delta_embedding_use", "embedding_shift"]
    for attr in candidate_attrs:
        if hasattr(oracle, attr):
            value = getattr(oracle, attr)
            if value is None:
                continue
            arr = np.asarray(value, dtype=float)
            if arr.ndim == 2 and arr.shape[0] == embedding.shape[0]:
                delta = arr[:, :2]
                break
    random_attrs = ["delta_embedding_random", "delta_embedding_randomized", "embedding_shift_random"]
    for attr in random_attrs:
        if hasattr(oracle, attr):
            value = getattr(oracle, attr)
            if value is None:
                continue
            arr = np.asarray(value, dtype=float)
            if arr.ndim == 2 and arr.shape[0] == embedding.shape[0]:
                delta_random = arr[:, :2]
                break
    if delta is None:
        raise AttributeError("Unable to locate embedding shift matrix after calculate_embedding_shift().")
    return embedding[:, :2], delta, delta_random


def make_score_table(oracle, tf: str, embedding_key: str, cluster_column: str) -> pd.DataFrame:
    embedding, delta, delta_random = extract_embedding_and_shift(oracle, embedding_key)
    obs = oracle.adata.obs.copy()
    obs.index = obs.index.astype(str)
    score_df = pd.DataFrame(
        {
            "cell_id": obs.index,
            "tf": tf,
            "cluster": obs[cluster_column].astype(str).values,
            "embedding_x": embedding[:, 0],
            "embedding_y": embedding[:, 1],
            "delta_x": delta[:, 0],
            "delta_y": delta[:, 1],
        }
    )
    if "group" in obs.columns:
        score_df["group"] = obs["group"].astype(str).values
    score_df["shift_length"] = np.sqrt(score_df["delta_x"] ** 2 + score_df["delta_y"] ** 2)
    score_df["shift_signed_x"] = score_df["delta_x"]
    score_df["shift_signed_y"] = score_df["delta_y"]
    if delta_random is not None:
        score_df["delta_random_x"] = delta_random[:, 0]
        score_df["delta_random_y"] = delta_random[:, 1]
        score_df["shift_length_random"] = np.sqrt(score_df["delta_random_x"] ** 2 + score_df["delta_random_y"] ** 2)
    return score_df


def summarize_score_table(score_df: pd.DataFrame, tf: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    group_col = "group" if "group" in score_df.columns else "cluster"
    summary = (
        score_df.groupby(group_col, dropna=False)
        .agg(
            n_cells=("cell_id", "size"),
            mean_shift_length=("shift_length", "mean"),
            median_shift_length=("shift_length", "median"),
            p90_shift_length=("shift_length", lambda x: float(np.quantile(x, 0.9))),
            mean_delta_x=("delta_x", "mean"),
            mean_delta_y=("delta_y", "mean"),
        )
        .reset_index()
        .rename(columns={group_col: "group"})
    )
    summary.insert(0, "tf", tf)

    overall = pd.DataFrame(
        [
            {
                "tf": tf,
                "n_cells": int(score_df.shape[0]),
                "mean_shift_length": float(score_df["shift_length"].mean()),
                "median_shift_length": float(score_df["shift_length"].median()),
                "p90_shift_length": float(np.quantile(score_df["shift_length"], 0.9)),
                "max_shift_length": float(score_df["shift_length"].max()),
            }
        ]
    )
    return summary, overall


def save_quiver_plots(oracle, tf: str, tf_dir: Path, quiver_scale: float) -> None:
    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    oracle.plot_quiver(scale=quiver_scale, ax=ax[0])
    ax[0].set_title(f"{tf} perturbation: simulated shift")
    oracle.plot_quiver_random(scale=quiver_scale, ax=ax[1])
    ax[1].set_title(f"{tf} perturbation: randomized control")
    fig.tight_layout()
    fig.savefig(tf_dir / f"{tf}_quiver.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def save_shift_scatter(score_df: pd.DataFrame, tf: str, tf_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    sc = ax.scatter(
        score_df["embedding_x"],
        score_df["embedding_y"],
        c=score_df["shift_length"],
        s=8,
        cmap="viridis",
        linewidths=0,
    )
    ax.set_title(f"{tf} perturbation: shift magnitude on embedding")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    cbar = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.03)
    cbar.set_label("Shift length")
    fig.tight_layout()
    fig.savefig(tf_dir / f"{tf}_shift_magnitude_embedding.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def save_shift_boxplot(score_df: pd.DataFrame, tf: str, tf_dir: Path) -> None:
    group_col = "group" if "group" in score_df.columns else "cluster"
    groups = list(pd.Series(score_df[group_col]).dropna().unique())
    data = [score_df.loc[score_df[group_col] == group, "shift_length"].to_numpy(dtype=float) for group in groups]
    fig, ax = plt.subplots(figsize=(7, 5))
    box = ax.boxplot(data, labels=groups, patch_artist=True, showfliers=False)
    colors = plt.cm.Set2(np.linspace(0.0, 1.0, len(groups)))
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.85)
    ax.set_ylabel("Shift length")
    ax.set_title(f"{tf} perturbation: shift magnitude by {group_col}")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(tf_dir / f"{tf}_shift_magnitude_boxplot.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def save_grid_flow_plot(oracle, tf: str, tf_dir: Path, n_grid: int, min_mass: float, neighbors: int, grid_scale: float) -> None:
    try:
        oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=neighbors)
        oracle.calculate_mass_filter(min_mass=min_mass, plot=False)
        fig, ax = plt.subplots(1, 2, figsize=(13, 6))
        oracle.plot_simulation_flow_on_grid(scale=grid_scale, ax=ax[0])
        ax[0].set_title(f"{tf} perturbation: flow on grid")
        oracle.plot_simulation_flow_random_on_grid(scale=grid_scale, ax=ax[1])
        ax[1].set_title(f"{tf} perturbation: randomized flow")
        fig.tight_layout()
        fig.savefig(tf_dir / f"{tf}_grid_flow.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

        fig2, ax2 = plt.subplots(figsize=(8, 8))
        oracle.plot_cluster_whole(ax=ax2, s=10)
        oracle.plot_simulation_flow_on_grid(scale=grid_scale, ax=ax2, show_background=False)
        ax2.set_title(f"{tf} perturbation: cluster background + flow")
        fig2.tight_layout()
        fig2.savefig(tf_dir / f"{tf}_grid_flow_with_clusters.png", dpi=200, bbox_inches="tight")
        plt.close(fig2)
    except Exception as exc:
        with open(tf_dir / f"{tf}_grid_flow_skipped.txt", "w", encoding="utf-8") as handle:
            handle.write(f"Grid-flow plot was skipped.\n{repr(exc)}\n")


def run_markov_and_save(oracle, tf: str, tf_dir: Path, cluster_column: str, n_steps: int, n_duplication: int, seed: int) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    try:
        oracle.prepare_markov_simulation(verbose=False)
        oracle.run_markov_chain_simulation(
            n_steps=n_steps,
            n_duplication=n_duplication,
            seed=seed,
            calculate_randomized=True,
        )

        transition_df = oracle.get_markov_simulation_cell_transition_table(
            cluster_column_name=cluster_column,
            end=-1,
            return_df=True,
        )
        if isinstance(transition_df, pd.DataFrame):
            transition_df.to_csv(tf_dir / f"{tf}_markov_transition_table.csv")

        summary_df = oracle.summarize_mc_results_by_cluster(cluster_use=cluster_column, random=False)
        if isinstance(summary_df, pd.DataFrame):
            summary_df.to_csv(tf_dir / f"{tf}_markov_summary_by_cluster.csv", index=False)

        try:
            summary_random_df = oracle.summarize_mc_results_by_cluster(cluster_use=cluster_column, random=True)
            if isinstance(summary_random_df, pd.DataFrame):
                summary_random_df.to_csv(tf_dir / f"{tf}_markov_summary_by_cluster_random.csv", index=False)
        except Exception:
            pass
        return transition_df if isinstance(transition_df, pd.DataFrame) else None, summary_df if isinstance(summary_df, pd.DataFrame) else None
    except Exception as exc:
        with open(tf_dir / f"{tf}_markov_skipped.txt", "w", encoding="utf-8") as handle:
            handle.write(f"Markov simulation was skipped.\n{repr(exc)}\n")
        return None, None


def main() -> None:
    args = parse_args()
    config_path = Path(args.config).resolve()
    output_dir = Path(args.output_dir).resolve()
    candidate_dir = Path(args.candidate_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print("[1/6] Importing CellOracle runtime")
    try:
        import celloracle as co
    except Exception as exc:
        raise RuntimeError(
            "This script requires a CellOracle runtime. Please run it inside the CellOracle Docker image."
        ) from exc

    with open(config_path, "r", encoding="utf-8") as handle:
        config = json.load(handle)

    embedding_key = str(config["embedding_key"])
    cluster_column = str(config["cluster_column"])
    tf_list = [item.strip() for item in args.tf_list.split(",") if item.strip()]

    candidate_metrics = load_candidate_metrics(candidate_dir)
    missing_tf = [tf for tf in tf_list if tf not in candidate_metrics["tf"].astype(str).tolist()]
    if missing_tf:
        raise ValueError(f"TFs missing from candidate metrics table: {missing_tf}")

    print("[2/6] Preparing simulation-ready Oracle object")
    sim_ready_path = ensure_simulation_ready_oracle(
        co=co,
        config=config,
        output_dir=output_dir,
        alpha=args.simulation_alpha,
        force_refit=args.force_refit,
    )

    combined_overall_rows: list[pd.DataFrame] = []
    combined_group_rows: list[pd.DataFrame] = []

    print("[3/6] Running in silico perturbation for prioritized TFs")
    for tf in tf_list:
        tf_dir = output_dir / safe_name(tf)
        tf_dir.mkdir(parents=True, exist_ok=True)

        oracle = co.load_hdf5(str(sim_ready_path))
        if tf not in oracle.adata.var_names:
            with open(tf_dir / f"{tf}_skipped.txt", "w", encoding="utf-8") as handle:
                handle.write(f"{tf} not found in Oracle gene space.\n")
            continue

        oracle.simulate_shift(
            perturb_condition={tf: 0.0},
            GRN_unit="cluster",
            n_propagation=args.n_propagation,
        )
        oracle.estimate_transition_prob(
            n_neighbors=min(args.transition_neighbors, max(int(oracle.adata.n_obs) - 1, 1)),
            knn_random=True,
            sampled_fraction=1,
        )
        oracle.calculate_embedding_shift(sigma_corr=args.sigma_corr)

        score_df = make_score_table(oracle, tf=tf, embedding_key=embedding_key, cluster_column=cluster_column)
        score_df.to_csv(tf_dir / f"{tf}_state_shift_scores.csv", index=False)

        group_summary_df, overall_df = summarize_score_table(score_df, tf=tf)
        group_summary_df.to_csv(tf_dir / f"{tf}_state_shift_summary_by_group.csv", index=False)
        overall_df.to_csv(tf_dir / f"{tf}_state_shift_summary_overall.csv", index=False)

        markov_table_df, markov_summary_df = run_markov_and_save(
            oracle=oracle,
            tf=tf,
            tf_dir=tf_dir,
            cluster_column=cluster_column,
            n_steps=args.markov_steps,
            n_duplication=args.markov_duplication,
            seed=args.seed,
        )

        save_quiver_plots(oracle, tf=tf, tf_dir=tf_dir, quiver_scale=args.quiver_scale)
        save_shift_scatter(score_df, tf=tf, tf_dir=tf_dir)
        save_shift_boxplot(score_df, tf=tf, tf_dir=tf_dir)
        save_grid_flow_plot(
            oracle=oracle,
            tf=tf,
            tf_dir=tf_dir,
            n_grid=args.n_grid,
            min_mass=args.min_mass,
            neighbors=min(args.transition_neighbors, max(int(oracle.adata.n_obs) - 1, 1)),
            grid_scale=args.grid_scale,
        )

        meta_row = candidate_metrics.loc[candidate_metrics["tf"].astype(str) == tf].iloc[0].to_dict()
        overall_df["expected_regulon_effect_group"] = meta_row.get("regulon_effect_group")
        overall_df["expected_expr_higher_group"] = meta_row.get("expr_higher_group")
        overall_df["regulon_fdr"] = meta_row.get("regulon_fdr")
        overall_df["expr_wilcoxon_fdr"] = meta_row.get("expr_wilcoxon_fdr")
        combined_overall_rows.append(overall_df)

        group_summary_df["expected_regulon_effect_group"] = meta_row.get("regulon_effect_group")
        group_summary_df["expected_expr_higher_group"] = meta_row.get("expr_higher_group")
        combined_group_rows.append(group_summary_df)

    print("[4/6] Writing combined summaries")
    if combined_overall_rows:
        pd.concat(combined_overall_rows, axis=0, ignore_index=True).to_csv(
            output_dir / "round1_ko_summary_overall.csv",
            index=False,
        )
    if combined_group_rows:
        pd.concat(combined_group_rows, axis=0, ignore_index=True).to_csv(
            output_dir / "round1_ko_summary_by_group.csv",
            index=False,
        )

    print("[5/6] Writing run metadata")
    run_meta = {
        "config": str(config_path),
        "simulation_ready_oracle": str(sim_ready_path),
        "tf_list": tf_list,
        "n_propagation": args.n_propagation,
        "transition_neighbors": args.transition_neighbors,
        "sigma_corr": args.sigma_corr,
        "markov_steps": args.markov_steps,
        "markov_duplication": args.markov_duplication,
        "output_dir": str(output_dir),
    }
    with open(output_dir / "round1_ko_run_config.json", "w", encoding="utf-8") as handle:
        json.dump(run_meta, handle, indent=2, ensure_ascii=False)

    print("[6/6] Done")
    print(f"  simulation_ready_oracle={sim_ready_path}")
    print(f"  outputs={output_dir}")


if __name__ == "__main__":
    main()
