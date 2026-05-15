#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

from matplotlib.lines import Line2D

import build_figS5_supportive_external_evidence_revised_v2 as base


base.STEM = "FigS5_supportive_external_evidence_revised_v2_minorfix"

base.FIGURE_TEXT_ITEMS = [
    "Supportive external evidence from GSE140393 and GSE190452",
    "A. GSE140393 single-group expression support",
    "Positive fraction",
    "Sample ID",
    "Candidate TF",
    "B. GSE140393 ranked single-group support score",
    "Support score",
    "C. GSE190452 cross-syndrome expression and target-set support",
    "Expression support score",
    "Target-set support score",
    "D. GSE190452 cluster-level expression-support pattern",
    "Directionally signed expression-support score",
    "Cluster",
    "E. Conservative interpretation of supportive evidence",
    *base.TF_ORDER,
    *base.INTERPRETATION.values(),
]


def build_panel_c_minorfix(ax: base.plt.Axes, round2: base.pd.DataFrame) -> None:
    required = ["tf", "expr_support_score", "regulon_support_score"]
    missing = [col for col in required if col not in round2.columns]
    if missing:
        raise ValueError(f"Missing required GSE190452 columns for Panel C: {missing}")
    plot_df = round2.set_index("tf").reindex(base.TF_ORDER).reset_index()
    plot_df["expr_support_score"] = base.as_num(plot_df["expr_support_score"])
    plot_df["target_set_support_score"] = base.as_num(plot_df["regulon_support_score"])

    ax.scatter(
        plot_df["expr_support_score"],
        plot_df["target_set_support_score"],
        s=96,
        c=[base.TF_COLORS[tf] for tf in plot_df["tf"]],
        edgecolor="white",
        linewidth=0.85,
        alpha=0.92,
        zorder=3,
    )
    offsets = {
        "NFE2L2": (-62, 13),
        "THRB": (8, 9),
        "BHLHE40": (8, -14),
        "SOX2": (8, 12),
    }
    for _, row in plot_df.iterrows():
        tf = str(row["tf"])
        ax.annotate(
            tf,
            (float(row["expr_support_score"]), float(row["target_set_support_score"])),
            xytext=offsets.get(tf, (5, 5)),
            textcoords="offset points",
            fontsize=6.6,
            color=base.TF_COLORS[tf],
            arrowprops={"arrowstyle": "-", "color": "#8A8A8A", "lw": 0.45, "shrinkA": 0, "shrinkB": 5},
        )
    ax.set_xlim(-0.05, 1.08)
    ax.set_ylim(-0.05, 1.08)
    ax.set_xlabel("Expression support score")
    ax.set_ylabel("Target-set support score")
    ax.set_title("GSE190452 cross-syndrome expression and target-set support", pad=7)
    base.style_axes(ax, grid_axis="both")
    base.panel_label(ax, "C")


def build_panel_e_minorfix(ax: base.plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.0, 0.86, "Conservative interpretation of supportive evidence", fontsize=8.5, fontweight="bold", ha="left", va="top")
    for i, tf in enumerate(base.TF_ORDER):
        x0 = 0.02 + i * 0.245
        ax.plot([x0, x0 + 0.055], [0.55, 0.55], color=base.TF_COLORS[tf], linewidth=3.2, solid_capstyle="round")
        ax.text(x0 + 0.065, 0.64, tf, fontsize=7.2, fontweight="bold", ha="left", va="center", color="#222222")
        ax.text(x0 + 0.065, 0.39, base.INTERPRETATION[tf], fontsize=7.0, ha="left", va="center", color="#333333")
    base.panel_label(ax, "E", x=-0.035, y=0.98)


def write_minorfix_notes() -> Path:
    path = base.OUT_DIR / f"{base.STEM}_notes.txt"
    residuals = residual_summary()
    lines = [
        "FigS5 supportive external evidence revised_v2 minorfix notes",
        "",
        "Output stem:",
        f"- {base.rel(base.OUT_DIR / base.STEM)}",
        "",
        "Input files used:",
        *[f"- {base.rel(p)}" for p in base.INPUT_FILES],
        "",
        "Minor fixes applied:",
        "- Panel C legend was removed because all points have direct TF labels.",
        "- Panel C xlim was set to -0.05 to 1.08 and ylim was set to -0.05 to 1.08.",
        "- Panel C NFE2L2 label was moved to the upper-left side of its point to avoid boundary crowding.",
        "- Panel C THRB, BHLHE40, and SOX2 labels were manually offset to avoid overlap.",
        "- Panel C title and axis labels were kept unchanged.",
        "- Panel D title, data source, and colorbar label were kept unchanged.",
        '- Panel E title was changed to "Conservative interpretation of supportive evidence".',
        "- Panel E candidate interpretation text was kept unchanged.",
        "",
        "Forbidden-term check against figure text:",
        *[f"- {term}: {status}" for term, status in residuals.items()],
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def residual_summary() -> dict[str, str]:
    figure_text = "\n".join(base.FIGURE_TEXT_ITEMS).lower()
    terms = [
        "validation",
        "formal validation",
        "independent validation",
        "external validation",
        "validated",
        "confirmed",
        "regulon replication",
        "projected regulon support",
        "internal_control",
    ]
    return {term: ("YES" if term in figure_text else "NO") for term in terms}


def required_term_summary() -> dict[str, str]:
    figure_text = "\n".join(base.FIGURE_TEXT_ITEMS).lower()
    terms = [
        "supportive external evidence",
        "single-group expression support",
        "cross-syndrome expression and target-set support",
        "cluster-level expression-support pattern",
        "internal-control",
    ]
    return {term: ("YES" if term in figure_text else "NO") for term in terms}


def main() -> None:
    base.build_panel_c = build_panel_c_minorfix
    base.build_panel_e = build_panel_e_minorfix
    outputs = base.build_figure()
    notes = write_minorfix_notes()
    residuals = residual_summary()
    required_terms = required_term_summary()

    print("FigS5 revised_v2 minorfix output path:")
    print(f"- {base.rel(base.OUT_DIR)}")
    for path in outputs:
        print(f"- {base.rel(path)}")
    print(f"- {base.rel(notes)}")

    print("\nPanel C legend deleted:")
    print("- YES")
    print("\nPanel C xlim/ylim adjusted:")
    print("- YES; xlim = -0.05 to 1.08; ylim = -0.05 to 1.08")
    print("\nPanel C NFE2L2 label position fixed:")
    print("- YES; label moved to the upper-left/left side of the point")
    print("\nPanel E title modified:")
    print('- YES; "Conservative interpretation of supportive evidence"')

    print("\nForbidden-term residuals in figure text:")
    for term, status in residuals.items():
        print(f"- {term}: {status}")

    print("\nRequired retained terms in figure text:")
    for term, status in required_terms.items():
        print(f"- {term}: {status}")


if __name__ == "__main__":
    main()
