# -*- coding: utf-8 -*-
"""
Extended Figure S3 (no overlap fix): re-layouts the standard Figure S3 panels
(A-F) and adds two permutation-null-distribution panels (G, H) in a single
unified 4-row gridspec so that every original panel remains fully visible.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(r"d:\docker_run_pyscenic\scripts\manuscript_figures")
SRC = SCRIPT_DIR / "build_figS3_celloracle_perturbation_metrics_revised_v2_minorfix.py"

spec = importlib.util.spec_from_file_location("figS3_base", SRC)
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

KO_DIR = Path(r"d:\docker_run_pyscenic\celloracle_run\ko_round1")
OUT_DIR = Path(r"d:\docker_run_pyscenic\manuscript_output\figures_supplementary_revised_v2")
OUT_DIR.mkdir(parents=True, exist_ok=True)
STEM = "FigS3_celloracle_perturbation_metrics_revised_v2_minorfix_with_permutation"

TFS = ["NFE2L2", "THRB", "BHLHE40", "SOX2"]
N_PERM = 10000
RNG = np.random.default_rng(20240817)
NULL_COLOR = "#9AA5B1"
OBS_COLOR = "#9A1F2D"


def load_scores(tf):
    return pd.read_csv(KO_DIR / tf / f"{tf}_state_shift_scores.csv")


def paired_effect(d_obs, d_rand, n_perm, rng):
    obs = np.linalg.norm(d_obs, axis=1)
    rand = np.linalg.norm(d_rand, axis=1)
    t_obs = float((obs - rand).mean())
    t_all = np.concatenate([obs[:, None], rand[:, None]], axis=1)
    flips = rng.integers(0, 2, size=(len(obs), n_perm))
    pick = np.where(flips == 0, t_all[:, 0:1], t_all[:, 1:2]).mean(axis=0)
    pick_rand = np.where(flips == 0, t_all[:, 1:2], t_all[:, 0:1]).mean(axis=0)
    null = pick - pick_rand
    p = float(max((null >= t_obs).sum() / n_perm, 1.0 / n_perm))
    z = (t_obs - null.mean()) / (null.std() + 1e-12)
    return t_obs, null, p, z


def group_direction(delta, group, n_perm, rng):
    g1, g2 = np.unique(group)[0], np.unique(group)[1]
    obs_diff = delta[group == g1].mean(axis=0) - delta[group == g2].mean(axis=0)
    axis = obs_diff / (np.linalg.norm(obs_diff) + 1e-12)
    t_obs = float(np.linalg.norm(obs_diff))
    null = np.empty(n_perm)
    for i in range(n_perm):
        perm = rng.permutation(len(group))
        d1 = delta[perm][group == g1].mean(axis=0)
        d2 = delta[perm][group == g2].mean(axis=0)
        null[i] = float((d1 - d2) @ axis)
    p = float(max((np.abs(null) >= np.abs(t_obs)).sum() / n_perm, 1.0 / n_perm))
    z = (np.abs(t_obs) - np.abs(null).mean()) / (np.abs(null).std() + 1e-12)
    return t_obs, null, p, z


def draw_null(ax, null, obs, p, z, title, xlabel):
    ax.hist(null, bins=40, color=NULL_COLOR, edgecolor="white", alpha=0.9)
    ax.axvline(obs, color=OBS_COLOR, lw=2.2)
    edge = OBS_COLOR if p < 0.05 else "#6B7280"
    pstr = "p < 1e-4" if p <= 1 / len(null) else f"p = {p:.5f}"
    ax.text(0.98, 0.97, f"{pstr}\nz = {z:.2f}", transform=ax.transAxes, ha="right",
            va="top", fontsize=9.5, bbox=dict(boxstyle="round,pad=0.3", fc="#F8FAFC", ec=edge, lw=1.1))
    ax.set_title(title, fontsize=8.2, pad=4)
    ax.set_xlabel(xlabel, fontsize=7.0)
    ax.set_ylabel("Permutation count", fontsize=7.0)
    ax.tick_params(labelsize=6.6)


def add_panel_label(ax, label, x=-0.08, y=1.06):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=13, fontweight="bold",
            ha="right", va="bottom", color="#111111")


def main():
    base.setup_style()
    metrics = base.load_metrics()
    vectors = base.load_vectors()
    xlim, ylim = base.compute_limits(vectors)

    fig = plt.figure(figsize=(13.3, 15.0), constrained_layout=False)
    gs = fig.add_gridspec(
        4, 2,
        height_ratios=[0.95, 1.55, 1.55, 1.45],
        hspace=0.42,
        wspace=0.28,
        left=0.065,
        right=0.985,
        top=0.95,
        bottom=0.05,
    )

    # row 0: A, B
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    base.build_panel_a(ax_a, metrics, fig)
    base.build_panel_b(ax_b, metrics)

    # row 1: C (NFE2L2), D (THRB)
    base.build_tf_block(fig, gs[1, 0], "NFE2L2", "C", vectors, xlim, ylim, seed=17)
    base.build_tf_block(fig, gs[1, 1], "THRB", "D", vectors, xlim, ylim, seed=23)
    # row 2: E (BHLHE40), F (SOX2)
    base.build_tf_block(fig, gs[2, 0], "BHLHE40", "E", vectors, xlim, ylim, seed=31)
    base.build_tf_block(fig, gs[2, 1], "SOX2", "F", vectors, xlim, ylim, seed=43)

    # row 3: G (effect), H (direction)
    ax_g = fig.add_subplot(gs[3, 0])
    ax_h = fig.add_subplot(gs[3, 1])

    all_eff = {}
    all_dir = {}
    for tf in TFS:
        df = load_scores(tf)
        d_obs = df[["delta_x", "delta_y"]].to_numpy(float)
        d_rand = df[["delta_random_x", "delta_random_y"]].to_numpy(float)
        group = df["group"].to_numpy()
        all_eff[tf] = paired_effect(d_obs, d_rand, N_PERM, RNG)
        all_dir[tf] = group_direction(d_obs, group, N_PERM, RNG)

    draw_null(ax_g, all_eff[TFS[0]][1], all_eff[TFS[0]][0], all_eff[TFS[0]][2], all_eff[TFS[0]][3],
              "G: paired observed-vs-randomized effect (all TFs)", "mean paired diff (obs − randomized)")
    for tf in TFS[1:]:
        ax_g.axvline(all_eff[tf][0], color=base.TF_COLORS[tf], lw=1.6, ls="--")
    add_panel_label(ax_g, "G")

    draw_null(ax_h, all_dir[TFS[0]][1], all_dir[TFS[0]][0], all_dir[TFS[0]][2], all_dir[TFS[0]][3],
              "H: group specificity (internal_control vs lesion)", "signed group mean-vector projection")
    for tf in TFS[1:]:
        ax_h.axvline(all_dir[tf][0], color=base.TF_COLORS[tf], lw=1.6, ls="--")
    add_panel_label(ax_h, "H")

    fig.suptitle("Complete CellOracle in silico perturbation metrics and randomized-control assessment",
                 x=0.5, y=0.975, fontsize=12.5, fontweight="bold")

    png = OUT_DIR / f"{STEM}.png"
    pdf = OUT_DIR / f"{STEM}.pdf"
    tiff = OUT_DIR / f"{STEM}.tiff"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    try:
        fig.savefig(tiff, dpi=450, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    except TypeError:
        fig.savefig(tiff, dpi=450, bbox_inches="tight")
    plt.close(fig)

    print("saved:", png, pdf, tiff)
    for tf in TFS:
        print(f"{tf}: eff z={all_eff[tf][3]:.2f} | dir z={all_dir[tf][3]:.2f}")


if __name__ == "__main__":
    main()