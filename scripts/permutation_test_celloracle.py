# -*- coding: utf-8 -*-
"""
Formal permutation test for CellOracle in silico perturbation results.

CellOracle provides, per cell, the observed perturbation shift vector (delta_x, delta_y)
and a paired randomized-control shift vector (delta_random_x, delta_random_y) obtained via a
randomized transition-probability (knn_random). We exploit this pairing to build principled
null distributions:

H0_effect (paired): the perturbation has no effect beyond a randomized network.
    Test statistic: mean over cells of |delta| - |delta_random|.
    Null: for each cell, independently flip the observed/random labels with probability 0.5
          (a pair-label permutation), so per-iteration we re-compute the mean paired difference
          under a faithful no-effect null. p = fraction of null where mean difference >= observed.

H0_direction (group specificity): the group difference in the mean perturbation vector
    (lesion vs internal_control) is no larger than expected by chance.
    Test statistic: magnitude of (mean_delta_lesion - mean_delta_internal_control).
    Null: permute group labels across cells, recompute the same statistic.

We also report z-scores against the null means/SDs and effect sizes (Cohen's d for the
observed-vs-random paired difference).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(r"d:\docker_run_pyscenic")
KO_DIR = ROOT / "celloracle_run" / "ko_round1"
OUT_DIR = KO_DIR / "permutation_test"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TFS = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]
N_PERM = 10000
RNG = np.random.default_rng(20240817)


def load_scores(tf: str) -> pd.DataFrame:
    df = pd.read_csv(KO_DIR / tf / f"{tf}_state_shift_scores.csv")
    for c in ("%s_x", "%s_y", "%s_x", "%s_y"):
        pass
    return df


def paired_effect_test(d_obs: np.ndarray, d_rand: np.ndarray, n_perm: int, rng: np.random.Generator) -> dict:
    """Paired permutation ('flip labels') test of observed vs randomized shift magnitude."""
    obs = np.linalg.norm(d_obs, axis=1)
    rand = np.linalg.norm(d_rand, axis=1)
    diff_obs = obs - rand
    t_obs = float(diff_obs.mean())

    null_t = np.empty(n_perm)
    t_all = np.concatenate([obs[:, None], rand[:, None]], axis=1)  # n x 2
    # shape (n_cells, n_perm): which of the two pair members is designated 'observed' in each perm
    flips = rng.integers(0, 2, size=(len(obs), n_perm))
    pick = np.where(flips == 0, t_all[:, 0:1], t_all[:, 1:2]).mean(axis=0)
    pick_rand = np.where(flips == 0, t_all[:, 1:2], t_all[:, 0:1]).mean(axis=0)
    null_t = pick - pick_rand
    # 'greater' p: probability null mean difference >= observed
    p = float((null_t >= t_obs).sum() / n_perm)
    # conservative estimate when zero permutations exceed observed (min resolvable = 1/n_perm)
    p = max(p, 1.0 / n_perm)
    z = (t_obs - null_t.mean()) / (null_t.std() + 1e-12)
    d = (np.asarray(diff_obs).mean()) / (np.asarray(diff_obs).std() + 1e-12)
    return {
        "paired_mean_diff_obs": t_obs,
        "paired_mean_diff_null_mean": float(null_t.mean()),
        "paired_mean_diff_null_sd": float(null_t.std()),
        "paired_p": p,
        "paired_z": float(z),
        "cohen_d": float(d),
        "obs_mean": float(obs.mean()),
        "rand_mean": float(rand.mean()),
        "obs_gt_rand_fraction": float((obs > rand).mean()),
        "null_diff": null_t,
        "n": len(obs),
    }


def group_direction_test(delta: np.ndarray, group: np.ndarray, n_perm: int, rng: np.random.Generator) -> dict:
    """Permute group labels to test the group specificity of the mean shift vector.

    Uses a *signed projection* statistic: project the group mean-vector difference onto the
    observed difference axis. This is scale-stable under permutation (unlike the bare norm,
    which is dominated by the within-group mean magnitude) and correctly separates 'the two
    groups differ' from 'both groups have large mean shifts'.
    """
    grps = np.unique(group)
    g1, g2 = grps[0], grps[1]
    obs_diff = delta[group == g1].mean(axis=0) - delta[group == g2].mean(axis=0)
    axis = obs_diff / (np.linalg.norm(obs_diff) + 1e-12)
    t_obs = float(np.linalg.norm(obs_diff))

    null_t = np.empty(n_perm)
    for i in range(n_perm):
        perm = rng.permutation(len(group))
        d1 = delta[perm][group == g1].mean(axis=0)
        d2 = delta[perm][group == g2].mean(axis=0)
        diff_p = d1 - d2
        # signed projection onto the observed axis (can be negative = opposite direction)
        null_t[i] = float(diff_p @ axis)
    # two-sided p: fraction of null where |signed projection| >= |observed projection|
    t_proj = float(obs_diff @ axis)  # == t_obs
    p = float((np.abs(null_t) >= np.abs(t_proj)).sum() / n_perm)
    p = max(p, 1.0 / n_perm)
    z = (np.abs(t_proj) - np.abs(null_t).mean()) / (np.abs(null_t).std() + 1e-12)
    return {
        "group_diff_obs": t_obs,
        "group_diff_null_mean": float(np.abs(null_t).mean()),
        "group_diff_null_sd": float(np.abs(null_t).std()),
        "group_diff_p": p,
        "group_diff_z": float(z),
        "null_group_diff": null_t,
        "g1": str(g1),
        "g2": str(g2),
    }


def plot_null(ax: plt.Axes, null: np.ndarray, obs: float, p: float, z: float, title: str,
              xlabel: str, obs_greater: bool = True) -> None:
    ax.hist(null, bins=40, color="#9AA5B1", edgecolor="white", alpha=0.9)
    ax.axvline(obs, color="#9A1F2D", lw=2.2, ls="-")
    ymax = ax.get_ylim()[1]
    xmin, xmax = ax.get_xlim()
    ax.set_xlim(min(xmin, obs), max(xmax, obs))
    ax.text(obs, ymax * 0.95, f"  observed\n  {obs:.4f}", color="#9A1F2D",
            fontsize=9, va="top", ha="left")
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel("Permutation count", fontsize=9)
    edge = '#9A1F2D' if p < 0.05 else '#6B7280'
    pstr = f"p < 1e-4" if p <= 1 / len(null) else f"p = {p:.5f}"
    ax.text(0.98, 0.96, f"{pstr}\nz = {z:.2f}", transform=ax.transAxes,
            ha="right", va="top", fontsize=11,
            bbox=dict(boxstyle="round,pad=0.35", fc="#F8FAFC", ec=edge, lw=1.2))


def main() -> None:
    rows = []
    fig, axes = plt.subplots(len(TFS), 2, figsize=(11, 3.2 * len(TFS)))
    fig.suptitle("CellOracle perturbation — permutation null distributions (N=10,000, paired observed-vs-randomized + group-label)",
                 fontsize=12.5, y=0.995)
    if len(TFS) == 1:
        axes = axes[None, :]

    for r, tf in enumerate(TFS):
        df = load_scores(tf)
        d_obs = df[["delta_x", "delta_y"]].to_numpy(dtype=float)
        d_rand = df[["delta_random_x", "delta_random_y"]].to_numpy(dtype=float)
        group = df["group"].to_numpy()

        eff = paired_effect_test(d_obs, d_rand, N_PERM, RNG)
        dir_ = group_direction_test(d_obs, group, N_PERM, RNG)

        rows.append({
            "tf": tf,
            "n_cells": eff["n"],
            "obs_mean_shift": eff["obs_mean"],
            "randomized_mean_shift": eff["rand_mean"],
            "obs_gt_rand_fraction": eff["obs_gt_rand_fraction"],
            "paired_mean_diff_observed": eff["paired_mean_diff_obs"],
            "paired_null_mean": eff["paired_mean_diff_null_mean"],
            "paired_null_sd": eff["paired_mean_diff_null_sd"],
            "paired_permutation_p": eff["paired_p"],
            "paired_z": eff["paired_z"],
            "cohen_d": eff["cohen_d"],
            "group_diff_observed": dir_["group_diff_obs"],
            "group_diff_null_mean": dir_["group_diff_null_mean"],
            "group_diff_null_sd": dir_["group_diff_null_sd"],
            "group_diff_permutation_p": dir_["group_diff_p"],
            "group_diff_z": dir_["group_diff_z"],
            "g1": dir_["g1"],
            "g2": dir_["g2"],
            "n_permutations": N_PERM,
        })

        plot_null(axes[r, 0], eff["null_diff"], eff["paired_mean_diff_obs"],
                  eff["paired_p"], eff["paired_z"],
                  f"{tf} — paired observed vs randomized (effect)", "mean paired diff (obs − randomized)")
        plot_null(axes[r, 1], dir_["null_group_diff"], dir_["group_diff_obs"],
                  dir_["group_diff_p"], dir_["group_diff_z"],
                  f"{tf} — group specificity ({dir_['g1']} vs {dir_['g2']})", "group mean-vector difference magnitude")

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT_DIR / "celloracle_permutation_null_distributions.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "celloracle_permutation_null_distributions.pdf", bbox_inches="tight")
    plt.close(fig)

    out_df = pd.DataFrame(rows)

    p_eff = out_df["paired_permutation_p"].to_numpy(dtype=float)
    p_dir = out_df["group_diff_permutation_p"].to_numpy(dtype=float)

    # Benjamini-Hochberg FDR across the four TFs for the paired effect p-values
    def bh_fdr(pvals):
        pv = np.asarray(pvals, dtype=float)
        m = len(pv)
        order = np.argsort(pv)
        p_sorted = pv[order]
        q_sorted = p_sorted * m / (np.arange(1, m + 1))
        q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
        q = np.empty_like(q_sorted)
        q[order] = q_sorted
        return q

    out_df["paired_permutation_FDR"] = bh_fdr(p_eff)
    out_df["group_diff_permutation_FDR"] = bh_fdr(p_dir)
    out_df.to_csv(OUT_DIR / "celloracle_permutation_test_results.csv", index=False)

    lines = ["CellOracle permutation test — summary", "=" * 62, "",
             f"Paired observed-vs-randomized test: per-cell flip of (observation, randomized) labels, "
             f"mean of |delta| - |delta_random|.",
             "Group-direction test: permutation of group labels, magnitude of group mean-vector difference.", "",
             "  tf      n      obs      rand   obs>rand   pairedDiff   p_effect     z_effect   groupDiff   p_dir     z_dir   cohenD"]
    for _, row in out_df.iterrows():
        pe = f"<1e-4" if row["paired_permutation_p"] <= 1 / float(row["n_permutations"]) else f"{row['paired_permutation_p']:.5f}"
        pd_ = f"<1e-4" if row["group_diff_permutation_p"] <= 1 / float(row["n_permutations"]) else f"{row['group_diff_permutation_p']:.5f}"
        lines.append(
            f"  {row['tf']:<7} {row['n_cells']:>6} {row['obs_mean_shift']:>7.3f} {row['randomized_mean_shift']:>7.3f} "
            f"{row['obs_gt_rand_fraction']:>9.3f} {row['paired_mean_diff_observed']:>10.3f} "
            f"{pe:>10} {row['paired_z']:>9.2f} {row['group_diff_observed']:>11.3f} "
            f"{pd_:>7} {row['group_diff_z']:>7.2f} {row['cohen_d']:>7.2f}")
    lines.append("")
    lines.append(f"Fixed seed 20240817, N_perm={N_PERM} per TF. 'obs>rand' = fraction of cells where observed "
                 f"shift exceeds its paired randomized control. p reported as '<1e-4' when no permutation "
                 f"reached the observed statistic (lower resolvable bound = 1/N_perm).")
    (OUT_DIR / "celloracle_permutation_test_summary.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()