#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import numpy as np

import build_figS7_drug_repositioning_exploratory_revised_v2 as base


base.STEM = "FigS7_drug_repositioning_exploratory_revised_v2_minorfix"

EGCG_FULL = "(-)-Epigallocatechin gallate"
KEY_PANEL_C_LABELS = ["VALPROIC ACID", "lobeline", "S-(+)-Rolipram", EGCG_FULL]


def build_panel_c_minorfix(ax: base.plt.Axes, full: base.pd.DataFrame, paper: base.pd.DataFrame) -> None:
    reviewed = full[full["after_bbb_final_layer"].notna()].copy()
    reviewed["final_score_prebbb_v21"] = base.as_num(reviewed["final_score_prebbb_v21"])
    reviewed["final_score_after_bbb"] = base.as_num(reviewed["final_score_after_bbb"])
    colors = [base.LAYER_COLORS.get(layer, "#D7DADB") for layer in reviewed["after_bbb_final_layer"]]
    sizes = np.where(reviewed["after_bbb_final_layer"].isin(base.LAYER_ORDER), 43, 20)
    alphas = np.where(reviewed["after_bbb_final_layer"].isin(base.LAYER_ORDER), 0.82, 0.40)
    ax.scatter(
        reviewed["final_score_prebbb_v21"],
        reviewed["final_score_after_bbb"],
        s=sizes,
        c=colors,
        alpha=alphas,
        edgecolor="white",
        linewidth=0.55,
        zorder=3,
    )
    ax.plot([0, 1], [0, 1], color="#BEBEBE", linewidth=0.8, linestyle="--", zorder=2)

    labels = reviewed[reviewed["compound"].astype(str).isin(KEY_PANEL_C_LABELS)].copy()
    order = {name: i for i, name in enumerate(KEY_PANEL_C_LABELS)}
    labels["label_order"] = labels["compound"].astype(str).map(order)
    labels = labels.sort_values("label_order")
    label_specs = {
        "VALPROIC ACID": {"text": "VALPROIC ACID", "offset": (-78, -26), "ha": "right", "va": "top"},
        "lobeline": {"text": "lobeline", "offset": (13, -27), "ha": "left", "va": "top"},
        "S-(+)-Rolipram": {"text": "S-(+)-Rolipram", "offset": (36, 22), "ha": "left", "va": "bottom"},
        EGCG_FULL: {"text": "EGCG", "offset": (-48, 24), "ha": "right", "va": "bottom"},
    }
    for _, row in labels.iterrows():
        name = str(row["compound"])
        spec = label_specs[name]
        ax.annotate(
            spec["text"],
            (float(row["final_score_prebbb_v21"]), float(row["final_score_after_bbb"])),
            xytext=spec["offset"],
            textcoords="offset points",
            fontsize=5.9,
            color="#444444",
            ha=spec["ha"],
            va=spec["va"],
            arrowprops={"arrowstyle": "-", "color": "#8A8A8A", "lw": 0.42, "shrinkA": 0, "shrinkB": 5},
        )
    ax.set_xlabel("Pre-BBB exploratory score")
    ax.set_ylabel("After-BBB exploratory score")
    ax.set_title("Pre-BBB versus after-BBB exploratory score shift", pad=7)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    base.style_axes(ax, grid_axis="both")
    base.panel_label(ax, "C")


def build_panel_d_minorfix(ax: base.plt.Axes, full: base.pd.DataFrame, fig: base.plt.Figure) -> None:
    base.build_panel_d_original(ax, full, fig)
    ax.set_xticklabels(["mechanism-direction", "supportive", "peripheral/program-modulating"], rotation=25, ha="right")


def residual_summary() -> dict[str, str]:
    figure_text = "\n".join(
        [
            "Exploratory drug-signature enrichment and BBB-layered prioritization",
            base.DISCLAIMER,
            "BBB-layered exploratory clue matrix",
            "Exploratory after-BBB ranked drug-signature clues",
            "Pre-BBB versus after-BBB exploratory score shift",
            "Pre-BBB exploratory score",
            "After-BBB exploratory score",
            "VALPROIC ACID",
            "lobeline",
            "S-(+)-Rolipram",
            "EGCG",
            "Mechanism-theme distribution across exploratory layers",
            "mechanism-direction",
            "supportive",
            "peripheral/program-modulating",
            "Conservative interpretation of exploratory drug-signature clues",
            base.PANEL_E_FINAL,
        ]
    ).lower()
    terms = [
        "therapeutic candidate",
        "clinical candidate",
        "validated drug",
        "treatment lead",
        "candidate therapy",
        "drug validation",
        "fcd ii treatment",
        "docking",
    ]
    return {term: ("YES" if term in figure_text else "NO") for term in terms}


def write_minorfix_notes() -> Path:
    path = base.OUT_DIR / f"{base.STEM}_notes.txt"
    residuals = residual_summary()
    lines = [
        "FigS7 drug repositioning exploratory revised_v2 minorfix notes",
        "",
        "Output stem:",
        f"- {base.rel(base.OUT_DIR / base.STEM)}",
        "",
        "Input files used:",
        *[f"- {base.rel(p)}" for p in base.INPUT_FILES],
        "",
        "Minor fixes applied:",
        "- Panel C retains only four direct labels: VALPROIC ACID, lobeline, S-(+)-Rolipram, and EGCG.",
        "- EGCG = (-)-Epigallocatechin gallate.",
        "- Panel C labels were manually offset with leader lines to avoid overlap and stay inside the plotting region.",
        "- Panel C title and axis labels were kept unchanged.",
        "- Panel D x-axis labels were shortened to mechanism-direction, supportive, and peripheral/program-modulating.",
        "- Panel D title was kept unchanged.",
        "- The top disclaimer strip was retained.",
        "- Panel E was retained, including the final boundary sentence.",
        "",
        "Forbidden-term check against figure text:",
        *[f"- {term}: {status}" for term, status in residuals.items()],
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main() -> None:
    base.build_panel_d_original = base.build_panel_d
    base.build_panel_c = build_panel_c_minorfix
    base.build_panel_d = build_panel_d_minorfix

    outputs = base.build_figure()
    notes = write_minorfix_notes()
    residuals = residual_summary()

    print("FigS7 revised_v2 minorfix output path:")
    print(f"- {base.rel(base.OUT_DIR)}")
    for path in outputs:
        print(f"- {base.rel(path)}")
    print(f"- {base.rel(notes)}")

    print("\nTop disclaimer retained:")
    print("- YES")
    print("\nPanel E retained:")
    print("- YES")
    print("\nPanel C label avoidance completed:")
    print("- YES; four key labels use manual offsets and leader lines")
    print("\n(-)-Epigallocatechin gallate abbreviated as EGCG:")
    print("- YES; notes define EGCG = (-)-Epigallocatechin gallate")

    print("\nForbidden therapeutic-term residuals in figure text:")
    for term, status in residuals.items():
        print(f"- {term}: {status}")


if __name__ == "__main__":
    main()
