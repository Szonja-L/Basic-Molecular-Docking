#!/usr/bin/env python3
"""
generate_figures.py
===================
Produces five publication-quality figures for the EGFR molecular docking
portfolio project.  Figures are saved to ./figures/ as 300-dpi PNGs.

Run after molecular_docking_pipeline.py (or standalone – data is embedded).
"""

import os
import math
import warnings
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import matplotlib.gridspec as gridspec
import seaborn as sns

warnings.filterwarnings("ignore")
matplotlib.rcParams.update({
    "font.family":       "DejaVu Sans",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.grid":         True,
    "grid.alpha":        0.25,
    "grid.linewidth":    0.6,
})

OUT_DIR = os.path.join(os.path.dirname(__file__), "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# SHARED DATA
# ─────────────────────────────────────────────────────────────────────────────
RT = 0.001987 * 298.15   # ≈ 0.5921 kcal/mol

compounds   = ["Aspirin", "Curcumin", "Quercetin", "Erlotinib",
               "Gefitinib", "Afatinib", "Lapatinib", "Osimertinib"]
energies    = [-5.2, -7.9, -8.4, -9.8, -10.2, -10.7, -11.3, -11.8]
ki_nm       = [math.exp(dg / RT) * 1e9 for dg in energies]
comp_types  = [
    "Negative Control",
    "Natural Compound", "Natural Compound",
    "FDA-Approved (1st Gen)", "FDA-Approved (1st Gen)",
    "FDA-Approved (2nd Gen)", "FDA-Approved (Dual)",
    "FDA-Approved (3rd Gen)",
]

# colour map: consistent palette throughout
TYPE_COLORS = {
    "Negative Control":          "#B0B0B0",
    "Natural Compound":          "#52B788",
    "FDA-Approved (1st Gen)":    "#4895EF",
    "FDA-Approved (2nd Gen)":    "#4361EE",
    "FDA-Approved (Dual)":       "#7B2D8B",
    "FDA-Approved (3rd Gen)":    "#E63946",
}
colors = [TYPE_COLORS[t] for t in comp_types]

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1 – Binding Affinity Comparison Bar Chart
# ─────────────────────────────────────────────────────────────────────────────
def fig1_binding_affinity():
    fig, ax = plt.subplots(figsize=(10, 7))

    bars = ax.barh(compounds, energies, color=colors, edgecolor="white",
                   linewidth=0.8, height=0.62, zorder=3)

    # Ki annotation on each bar
    for i, (bar, ki) in enumerate(zip(bars, ki_nm)):
        ki_str = f"{ki:.0f} nM" if ki >= 1 else f"{ki*1000:.0f} pM"
        ax.text(bar.get_width() - 0.08, bar.get_y() + bar.get_height() / 2,
                ki_str, va="center", ha="right", fontsize=8.5,
                color="white", fontweight="bold")

    # Threshold lines
    ax.axvline(-10.0, color="#E63946", ls="--", lw=1.4, alpha=0.7, zorder=2,
               label="High-affinity threshold (−10 kcal/mol)")
    ax.axvline(-8.0,  color="#F4A261", ls=":",  lw=1.4, alpha=0.7, zorder=2,
               label="Moderate-affinity threshold (−8 kcal/mol)")

    # Legend for compound types
    legend_patches = [mpatches.Patch(color=c, label=l)
                      for l, c in TYPE_COLORS.items()]
    legend_patches += [
        plt.Line2D([0], [0], color="#E63946", ls="--", lw=1.4,
                   label="High-affinity threshold (−10 kcal/mol)"),
        plt.Line2D([0], [0], color="#F4A261", ls=":",  lw=1.4,
                   label="Moderate-affinity (−8 kcal/mol)"),
    ]
    ax.legend(handles=legend_patches, loc="lower right",
              fontsize=8, framealpha=0.9, edgecolor="#cccccc")

    ax.set_xlabel("Binding Free Energy  ΔG (kcal/mol)", fontsize=12)
    ax.set_title("EGFR Kinase Binding Affinity — AutoDock Vina Scores\n"
                 "Target: EGFR Tyrosine Kinase (PDB: 1IEP)", fontsize=13, pad=14)
    ax.set_xlim(-13.5, 0)
    ax.tick_params(axis="both", labelsize=10)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "01_binding_affinity_comparison.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓  Figure 1 saved → {path}")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2 – Drug-Likeness Radar Chart
# ─────────────────────────────────────────────────────────────────────────────
def fig2_drug_likeness_radar():
    """
    Radar chart of Lipinski / Veber properties normalised to their upper
    limits.  A score of 1.0 = exactly at the threshold; >1.0 = violation.
    """
    radar_compounds = {
        "Erlotinib":   dict(MW=393.44/500, LogP=2.70/5, HBD=1/5, HBA=7/10, TPSA=74.9/140, RotB=8/10),
        "Gefitinib":   dict(MW=446.90/500, LogP=3.74/5, HBD=1/5, HBA=8/10, TPSA=68.7/140, RotB=7/10),
        "Afatinib":    dict(MW=485.94/500, LogP=4.89/5, HBD=2/5, HBA=9/10, TPSA=99.2/140, RotB=9/10),
        "Lapatinib":   dict(MW=581.06/500, LogP=5.11/5, HBD=2/5, HBA=8/10, TPSA=100.4/140, RotB=8/10),
        "Osimertinib": dict(MW=499.61/500, LogP=4.94/5, HBD=2/5, HBA=7/10, TPSA=91.6/140, RotB=10/10),
        "Quercetin":   dict(MW=302.24/500, LogP=1.54/5, HBD=5/5, HBA=7/10, TPSA=131.4/140, RotB=1/10),
        "Curcumin":    dict(MW=368.38/500, LogP=3.29/5, HBD=2/5, HBA=6/10, TPSA=93.1/140, RotB=8/10),
    }
    labels = ["MW\n(/500)", "LogP\n(/5)", "HBD\n(/5)", "HBA\n(/10)", "TPSA\n(/140)", "RotB\n(/10)"]
    N      = len(labels)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]

    # Distinct colours – limit to 7 FDA / natural (skip aspirin in radar)
    radar_colors = ["#4895EF", "#3A86FF", "#4361EE", "#7B2D8B",
                    "#E63946", "#52B788", "#2D9E6B"]

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))

    for (name, props), col in zip(radar_compounds.items(), radar_colors):
        vals  = list(props.values()) + [list(props.values())[0]]
        ax.plot(angles, vals, "o-", lw=1.8, color=col, label=name, alpha=0.9)
        ax.fill(angles, vals, alpha=0.06, color=col)

    # Threshold ring at 1.0
    theta_ring = np.linspace(0, 2 * np.pi, 200)
    ax.plot(theta_ring, [1.0] * 200, "--", color="#E63946", lw=1.6,
            alpha=0.8, label="Lipinski threshold = 1.0")

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, fontsize=10.5)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0, 1.25])
    ax.set_yticklabels(["0.25", "0.50", "0.75", "1.00 ⚠", "1.25"], fontsize=8)
    ax.set_ylim(0, 1.35)
    ax.set_title("Drug-Likeness Radar\nNormalised Lipinski / Veber Properties",
                 fontsize=13, pad=22)
    ax.legend(loc="upper right", bbox_to_anchor=(1.32, 1.12),
              fontsize=9, framealpha=0.9)

    path = os.path.join(OUT_DIR, "02_drug_likeness_radar.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓  Figure 2 saved → {path}")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 3 – Protein–Ligand Interaction Heatmap
# ─────────────────────────────────────────────────────────────────────────────
def fig3_interaction_heatmap():
    residues = ["Met793", "Thr790", "Lys745", "Glu762",
                "Cys797", "Val726", "Ala743", "Phe699", "Asp855"]
    # ordered best-to-worst binder (top of heatmap = strongest binder)
    ordered_cmpds = ["Osimertinib", "Lapatinib", "Afatinib",
                     "Gefitinib", "Erlotinib", "Quercetin", "Curcumin", "Aspirin"]
    matrix = {
        "Erlotinib":   [2, 1, 0, 1, 0, 1, 1, 1, 0],
        "Gefitinib":   [2, 1, 0, 1, 0, 1, 1, 1, 0],
        "Afatinib":    [2, 1, 1, 1, 4, 1, 1, 0, 0],
        "Lapatinib":   [2, 0, 1, 1, 0, 1, 1, 3, 1],
        "Osimertinib": [2, 0, 1, 1, 4, 1, 1, 1, 1],
        "Quercetin":   [1, 2, 0, 0, 0, 1, 0, 3, 0],
        "Curcumin":    [1, 0, 0, 1, 0, 1, 1, 0, 0],
        "Aspirin":     [0, 0, 0, 0, 0, 1, 1, 0, 1],
    }
    data = np.array([matrix[c] for c in ordered_cmpds])

    fig, ax = plt.subplots(figsize=(11, 6))
    # Custom discrete colormap
    cmap = matplotlib.colors.ListedColormap(
        ["#F5F5F5", "#A8D5BA", "#4895EF", "#1A5EB8", "#E63946"])
    norm = matplotlib.colors.BoundaryNorm([0, 1, 2, 3, 4, 5], cmap.N)

    im = ax.imshow(data, cmap=cmap, norm=norm, aspect="auto")

    ax.set_xticks(range(len(residues)))
    ax.set_xticklabels(residues, fontsize=10.5, rotation=35, ha="right")
    ax.set_yticks(range(len(ordered_cmpds)))
    ax.set_yticklabels(ordered_cmpds, fontsize=10.5)

    # Value annotations
    interaction_labels = {0: "", 1: "Hydrophobic", 2: "H-bond",
                          3: "π–π", 4: "Covalent"}
    symbol = {0: "", 1: "HΦ", 2: "H", 3: "π", 4: "COV"}
    for i in range(len(ordered_cmpds)):
        for j in range(len(residues)):
            val = data[i, j]
            if val > 0:
                txt_color = "white" if val >= 3 else "#1a1a2e"
                ax.text(j, i, symbol[val], ha="center", va="center",
                        fontsize=9.5, fontweight="bold", color=txt_color)

    # Colorbar legend
    cbar = plt.colorbar(im, ax=ax, ticks=[0.5, 1.5, 2.5, 3.5, 4.5],
                        fraction=0.025, pad=0.04)
    cbar.ax.set_yticklabels(["None", "Hydrophobic", "H-bond", "π–π stack", "Covalent"],
                             fontsize=8.5)

    # Hinge region highlight
    ax.add_patch(FancyBboxPatch((-0.5, -0.5), 1, len(ordered_cmpds),
                                boxstyle="round,pad=0.05",
                                linewidth=2.2, edgecolor="#E63946",
                                facecolor="none", clip_on=False,
                                zorder=5))
    ax.text(0, -1.1, "Hinge", ha="center", fontsize=8.5,
            color="#E63946", fontweight="bold")

    ax.set_title("Protein–Ligand Interaction Fingerprint\nEGFR Kinase Active Site Residues (PDB: 1IEP)",
                 fontsize=13, pad=14)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "03_interaction_heatmap.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓  Figure 3 saved → {path}")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 4 – ΔG vs. Ki Scatter (Theoretical Curve + Compound Points)
# ─────────────────────────────────────────────────────────────────────────────
def fig4_dg_ki_scatter():
    fig, ax = plt.subplots(figsize=(9, 6))

    # Theoretical curve
    dg_range = np.linspace(-13, -3, 300)
    ki_range = np.array([math.exp(dg / RT) * 1e9 for dg in dg_range])
    ax.semilogy(dg_range, ki_range, "-", color="#AAAAAA", lw=2.2,
                zorder=1, label="Theoretical: Ki = exp(ΔG/RT)")

    # Data points
    scatter = ax.scatter(energies, ki_nm, c=colors, s=140,
                         zorder=5, edgecolors="white", linewidths=0.9)

    # Labels
    offsets = {
        "Osimertinib": (-0.05,  1.4),
        "Lapatinib":   (-0.05,  1.4),
        "Afatinib":    (-0.05,  1.4),
        "Gefitinib":   ( 0.2,   0.5),
        "Erlotinib":   ( 0.2,   0.5),
        "Quercetin":   ( 0.2,   0.5),
        "Curcumin":    ( 0.2,   0.5),
        "Aspirin":     ( 0.2,   0.5),
    }
    for name, e, ki in zip(compounds, energies, ki_nm):
        dx, dy_factor = offsets[name]
        ax.annotate(name, xy=(e, ki),
                    xytext=(e + dx, ki * dy_factor if dy_factor > 1 else ki / 1.8),
                    fontsize=8.5, color="#1a1a2e",
                    arrowprops=dict(arrowstyle="-", color="#AAAAAA", lw=0.8))

    # Potency zones
    ax.axhline(10,   color="#E63946", ls="--", lw=1.2, alpha=0.8,
               label="10 nM  (high potency)")
    ax.axhline(1000, color="#F4A261", ls=":",  lw=1.2, alpha=0.8,
               label="1 μM  (moderate potency)")

    ax.set_xlabel("Binding Free Energy  ΔG (kcal/mol)", fontsize=12)
    ax.set_ylabel("Inhibition Constant  Ki (nM)  [log scale]", fontsize=12)
    ax.set_title("ΔG vs. Inhibition Constant\nRelationship: Ki = exp(ΔG / RT)", fontsize=13, pad=14)
    ax.set_xlim(-13, -3.5)
    ax.set_ylim(0.5, 5e5)
    ax.invert_xaxis()

    # Custom legend for compound types
    legend_patches = [mpatches.Patch(color=c, label=l) for l, c in TYPE_COLORS.items()]
    legend_patches += [
        plt.Line2D([0], [0], color="#AAAAAA", lw=2,   label="Theoretical curve"),
        plt.Line2D([0], [0], color="#E63946", ls="--", lw=1.4, label="10 nM threshold"),
        plt.Line2D([0], [0], color="#F4A261", ls=":",  lw=1.4, label="1 μM threshold"),
    ]
    ax.legend(handles=legend_patches, fontsize=8, loc="upper right",
              framealpha=0.9, edgecolor="#cccccc")
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "04_ki_scatter_plot.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓  Figure 4 saved → {path}")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 5 – ADMET / Physicochemical Property Bubble Chart
# ─────────────────────────────────────────────────────────────────────────────
def fig5_admet_bubble():
    data = {
        "Compound":   ["Erlotinib", "Gefitinib", "Afatinib", "Lapatinib",
                       "Osimertinib", "Quercetin", "Curcumin", "Aspirin"],
        "MW":         [393.44, 446.90, 485.94, 581.06, 499.61, 302.24, 368.38, 180.16],
        "LogP":       [2.70,   3.74,   4.89,   5.11,   4.94,   1.54,   3.29,   1.19],
        "TPSA":       [74.9,   68.7,   99.2,   100.4,  91.6,   131.4,  93.1,   63.6],
        "BindE":      [-9.8,  -10.2,  -10.7,  -11.3,  -11.8,  -8.4,   -7.9,   -5.2],
        "Type":       comp_types,
        "RO5_pass":   [True,   True,   True,   False,  True,   True,   True,   True],
    }
    df = pd.DataFrame(data)
    # Bubble size proportional to |ΔG|
    sizes = (np.abs(df["BindE"]) ** 2) * 12

    fig, ax = plt.subplots(figsize=(10, 7))
    for _, row in df.iterrows():
        edge = "#333333" if row["RO5_pass"] else "#E63946"
        lw   = 0.8       if row["RO5_pass"] else 2.2
        ax.scatter(row["LogP"], row["MW"],
                   s=sizes[_], color=TYPE_COLORS[row["Type"]],
                   edgecolors=edge, linewidths=lw,
                   alpha=0.85, zorder=5)
        ax.annotate(row["Compound"], xy=(row["LogP"], row["MW"]),
                    xytext=(row["LogP"] + 0.10, row["MW"] + 10),
                    fontsize=8.5, color="#1a1a2e")

    # Lipinski limit lines
    ax.axvline(5,   color="#E63946", ls="--", lw=1.4, alpha=0.8,
               label="LogP = 5  (Lipinski limit)")
    ax.axhline(500, color="#E63946", ls=":",  lw=1.4, alpha=0.8,
               label="MW = 500  (Lipinski limit)")

    ax.set_xlabel("LogP  (lipophilicity)", fontsize=12)
    ax.set_ylabel("Molecular Weight  (g/mol)", fontsize=12)
    ax.set_title("Physicochemical Space — Lipinski Rule of Five\n"
                 "Bubble size ∝ |binding energy|  |  Red border = RO5 violation",
                 fontsize=12, pad=14)
    ax.set_xlim(0, 6.5)
    ax.set_ylim(100, 680)

    # Legend
    legend_patches = [mpatches.Patch(color=c, label=l) for l, c in TYPE_COLORS.items()]
    legend_patches += [
        plt.Line2D([0], [0], color="#E63946", ls="--", lw=1.4, label="LogP = 5"),
        plt.Line2D([0], [0], color="#E63946", ls=":",  lw=1.4, label="MW = 500"),
        mpatches.Patch(facecolor="white", edgecolor="#E63946", linewidth=2.2, label="RO5 violation"),
    ]
    ax.legend(handles=legend_patches, fontsize=8, loc="upper left",
              framealpha=0.9, edgecolor="#cccccc")
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "05_admet_bubble_chart.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓  Figure 5 saved → {path}")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Generating figures...")
    fig1_binding_affinity()
    fig2_drug_likeness_radar()
    fig3_interaction_heatmap()
    fig4_dg_ki_scatter()
    fig5_admet_bubble()
    print(f"\nAll figures saved to {OUT_DIR}")
