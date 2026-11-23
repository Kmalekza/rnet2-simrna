#!/usr/bin/env python3
"""Figures and PDF report for RibonanzaNet2 + SimRNA on 6TFF, 6WJR, 6Q57.  (needs matplotlib, pandas, reportlab)

  python analysis/make_report.py  ->  results/figures/*.png, results/report.pdf
Inputs (all in this repository): results/validation/{scores.csv, frames_rmsd.csv, pairs_by_arm.json, snapshots.csv.gz},
results/pipeline/*/report.json. RNAnneal numbers are hard-coded below from its preprint (see README).
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
RES = REPO / "results"; OUT = RES; FIG = RES / "figures"; FIG.mkdir(exist_ok=True)
TGT = ["6TFF", "6WJR", "6Q57"]
ARMS = ["nocon", "eterna", "rnet2", "xtal"]
LAB = {"nocon": "no pairs", "eterna": "EternaFold", "rnet2": "RibonanzaNet2", "xtal": "crystal pairs"}
COL = {"nocon": "#9aa5a1", "eterna": "#7d8f3a", "rnet2": "#2f6f9f", "xtal": "#c2693a", "rnanneal": "#7a5c99"}
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 200})

sc = pd.read_csv(RES / "validation" / "scores.csv")
top = sc[sc.cluster == 1].set_index(["target", "arm"])
fr = pd.read_csv(RES / "validation" / "frames_rmsd.csv").set_index(["target", "arm"])
pairs = json.load(open(RES / "validation" / "pairs_by_arm.json"))
pipe = {t: json.load(open(RES / "pipeline" / t / "report.json")) for t in TGT}
sn = pd.read_csv(RES / "validation" / "snapshots.csv.gz")
snap = {(t, a): (g.energy.to_numpy(), g.rmsd_c4.to_numpy()) for (t, a), g in sn.groupby(["target", "arm"])}

# ------------------------------------------------------------------ figure 1: final models per arm
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))
x = np.arange(len(TGT)); w = 0.2
for k, (m, lab) in enumerate([("TM", "TM-score (largest cluster)"), ("RMSD", "RMSD, whole chain (Å)"), ("INF", "INF (base-pair agreement)")]):
    for j, a in enumerate(ARMS):
        ax[k].bar(x + (j - 1.5) * w, [top.loc[(t, a), m] for t in TGT], w, color=COL[a], label=LAB[a])
    ax[k].set_xticks(x, TGT); ax[k].set_title(lab, fontsize=9)
ax[0].axhline(0.45, ls="--", lw=0.8, color="k"); ax[0].text(2.45, 0.46, "same fold", ha="right", fontsize=7)
ax[0].set_ylim(0, 0.6); ax[2].set_ylim(0, 1.05)
ax[0].legend(fontsize=7, frameon=False, loc="upper left", ncol=2)
fig.tight_layout(); fig.savefig(FIG / "fig1_arms.png"); plt.close(fig)

# ------------------------------------------------------------------ figure 2: all snapshots, RMSD distributions
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6), sharey=False)
for k, t in enumerate(TGT):
    data = [snap[(t, a)][1] for a in ARMS]
    vp = ax[k].violinplot(data, positions=range(4), showextrema=False, widths=0.85)
    for b, a in zip(vp["bodies"], ARMS):
        b.set_facecolor(COL[a]); b.set_alpha(0.55); b.set_edgecolor("none")
    ax[k].scatter(range(4), [d.min() for d in data], marker="v", color="k", s=14, zorder=3, label="closest snapshot")
    ax[k].scatter(range(4), [top.loc[(t, a), "RMSD"] for a in ARMS], marker="o", facecolor="white",
                  edgecolor="k", s=18, zorder=3, label="chosen model (largest cluster)")
    ax[k].axhline(5, ls="--", lw=0.8, color="k")
    ax[k].set_xticks(range(4), ["none", "Eterna", "RNet2", "crystal"], fontsize=7.5)
    ax[k].set_title(f"{t}: all 32,080 snapshots", fontsize=9); ax[k].set_ylim(0, None)
ax[0].set_ylabel("RMSD to crystal (Å)"); ax[0].legend(fontsize=6.5, frameon=False, loc="upper right")
fig.tight_layout(); fig.savefig(FIG / "fig2_snapshots.png"); plt.close(fig)

# ------------------------------------------------------------------ figure 3: energy vs RMSD, RibonanzaNet2 arm
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))
for k, t in enumerate(TGT):
    E, r = snap[(t, "rnet2")]; low = E <= np.quantile(E, 0.01)
    ax[k].scatter(r[~low], E[~low], s=1, color="#b8c4cf", rasterized=True, label="all snapshots")
    ax[k].scatter(r[low], E[low], s=3, color=COL["rnet2"], rasterized=True, label="lowest-energy 1%")
    rho = pd.Series(E).corr(pd.Series(r), method="spearman")
    ax[k].set_title(f"{t}  (rank corr. {rho:.2f})", fontsize=9); ax[k].set_xlabel("RMSD to crystal (Å)")
    ax[k].axvline(5, ls="--", lw=0.8, color="k")
ax[0].set_ylabel("SimRNA energy"); ax[0].legend(fontsize=6.5, frameon=False, markerscale=3, loc="upper right")
fig.tight_layout(); fig.savefig(FIG / "fig3_energy.png"); plt.close(fig)

# ------------------------------------------------------------------ figure 4: comparison with RNAnneal
rn_top = {"6TFF": 14.0, "6WJR": np.nan, "6Q57": 10.0}        # 6Q57: "> 10 A" (lower bound); 6WJR not resolvable
rn_best = {"6TFF": 2.29, "6WJR": 15.0, "6Q57": 7.5}
ours_top = {t: top.loc[(t, "rnet2"), "RMSD"] for t in TGT}
ours_best = {t: fr.loc[(t, "rnet2"), "best_sampled"] for t in TGT}
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.6))
w = 0.36
for k, (title, rn, ou) in enumerate([("Chosen answer", rn_top, ours_top), ("Closest shape generated", rn_best, ours_best)]):
    b1 = ax[k].bar(x - w / 2, [rn[t] for t in TGT], w, color=COL["rnanneal"], label="RNAnneal (preprint)")
    b2 = ax[k].bar(x + w / 2, [ou[t] for t in TGT], w, color=COL["rnet2"], label="RibonanzaNet2 + SimRNA")
    for xi, t in zip(x, TGT):
        if np.isnan(rn[t]):
            ax[k].text(xi - w / 2, 0.5, "unclear", rotation=90, ha="center", va="bottom", fontsize=7)
    if k == 0:
        ax[k].text(2 - w / 2, rn_top["6Q57"] + 0.4, "> 10", ha="center", fontsize=7)
    ax[k].axhline(5, ls="--", lw=0.8, color="k"); ax[k].set_xticks(x, TGT); ax[k].set_title(title, fontsize=9)
ax[0].set_ylabel("RMSD to crystal (Å)"); ax[1].legend(fontsize=7, frameon=False, loc="upper left")
fig.tight_layout(); fig.savefig(FIG / "fig4_rnanneal.png"); plt.close(fig)

# ------------------------------------------------------------------ PDF
from reportlab.lib import colors  # noqa: E402
from reportlab.lib.pagesizes import letter  # noqa: E402
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet  # noqa: E402
from reportlab.lib.units import inch  # noqa: E402
from reportlab.platypus import Image, KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle  # noqa: E402

ss = getSampleStyleSheet()
H1 = ParagraphStyle("h1", parent=ss["Heading1"], fontSize=15, spaceAfter=6)
H2 = ParagraphStyle("h2", parent=ss["Heading2"], fontSize=11.5, spaceBefore=10, spaceAfter=4, textColor=colors.HexColor("#17201d"))
B = ParagraphStyle("b", parent=ss["BodyText"], fontSize=9, leading=12.2, spaceAfter=5)
S = ParagraphStyle("s", parent=B, fontSize=7.8, leading=10, textColor=colors.HexColor("#5b6a65"))
M = ParagraphStyle("m", parent=B, fontName="Courier", fontSize=7, leading=8.6)
W = 7.0 * inch


def tbl(rows, widths=None, hl=None):
    t = Table([[Paragraph(str(c), S if i else ParagraphStyle("th", parent=S, fontName="Helvetica-Bold")) for c in r]
               for i, r in enumerate(rows)], colWidths=widths, repeatRows=1)
    st = [("LINEBELOW", (0, 0), (-1, 0), 0.8, colors.HexColor("#17201d")), ("LINEBELOW", (0, 1), (-1, -1), 0.25, colors.HexColor("#c9d3cf")),
          ("VALIGN", (0, 0), (-1, -1), "TOP"), ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]
    for (c, r) in (hl or []):
        st.append(("BACKGROUND", (c, r), (c, r), colors.HexColor("#e3eef6")))
    t.setStyle(TableStyle(st)); return t


def fig(name, cap, h):
    return KeepTogether([Image(str(FIG / name), width=W, height=h), Paragraph(cap, S), Spacer(1, 6)])


f2 = lambda v: f"{v:.2f}"; f1 = lambda v: f"{v:.1f}"
story = [Paragraph("RibonanzaNet2 + SimRNA: sequence-only RNA 3D prediction", H1),
         Paragraph("Riboswitches 6TFF, 6WJR and 6Q57, compared with RNAnneal. All numbers are from results/validation and "
                   "results/pipeline in this repository; RNAnneal numbers are from its preprint.", S),
         Spacer(1, 6),
         Paragraph("Summary", H2),
         Paragraph("The protocol takes only an RNA sequence. RibonanzaNet2, a frozen RNA language model, predicts the base pairs, "
                   "including pseudoknots; SimRNA, a coarse-grained Monte Carlo folding program, builds 3D structures under those pairs; "
                   "the centre of the largest cluster of low-energy shapes is the prediction. The crystal structure is used only to score "
                   "the result.", B),
         Paragraph("<b>Findings.</b> (1) RibonanzaNet2 pairs reach the 3D models: higher base-pair agreement with the crystal than with EternaFold "
                   "(INF 0.72 vs 0.67 on average), and the only predicted source that builds 6Q57's pseudoknot (5 of 5 pairs). "
                   "(2) No model has the right overall fold (all TM-scores below 0.45); even the crystal's own pairs leave the models "
                   "about 12 Å off on average. The limit is the 3D arrangement of the helices, not the pairs. "
                   "(3) Across all 32,080 snapshots per run, a shape within 5 Å of the crystal appears only for 6TFF, and SimRNA's "
                   "energy barely ranks shapes by correctness. (4) Against RNAnneal on the same three RNAs, RibonanzaNet2 + SimRNA reaches a little "
                   "closer on 6WJR and 6Q57 and matches it on 6TFF, while RNAnneal's learned score chooses far better. Combining the two "
                   "addresses both weak points.", B),
         Paragraph("Protocol", H2),
         tbl([["Step", "What it does", "Settings"],
              ["1  Base pairs", "RibonanzaNet2 (48 layers, frozen); attention head layer 25 head 1, chosen once on 612 training RNAs "
               "(SPOT-RNA PDB TR1 + 500 bpRNA TR0; these include 4LVW, a 99% copy of 6Q57); decoded greedily", "canonical pairs only, |i-j| >= 4, one partner each, "
               "round(0.3 L) pairs, crossing pairs kept as pseudoknots"],
              ["2  3D folding", "SimRNA 3.20 replica-exchange Monte Carlo; pairs given as secondary-structure restraints, one "
               "dot-bracket line per pseudoknot layer", "8 runs x 10 replicas x 1.6 M steps (10% of the SimRNA benchmark length); "
               "32,080 snapshots"],
              ["3  Choose", "pool snapshots, keep the lowest-energy 1%, cluster, rebuild all atoms", "cluster at 0.1 L Å; "
               "model_1 = largest cluster, model_2..5 = next"],
              ["4  Score (optional)", "against an experimental structure of the same RNA", "TM-score and core RMSD (US-align, C1'), "
               "whole-chain RMSD (residue by residue), INF = sqrt(precision x recall) of canonical pairs (rnapolis)"]],
             [1.05 * inch, 3.2 * inch, 2.75 * inch]),
         Spacer(1, 4),
         Paragraph("Command: <font name='Courier'>bash pipeline/run.sh NAME SEQUENCE [REFERENCE.cif CHAIN]</font>. Three Slurm "
                   "jobs (prep 1 min, 8 folding runs of 10-40 min in parallel, post 1-5 min); outputs model_1..5.pdb, pairs.txt, "
                   "report.md/json.", B),
         Paragraph("Pipeline results (sequence only)", H2)]
rows = [["RNA", "Pairs predicted", "model_1 cluster", "TM", "RMSD Å", "core RMSD Å", "INF", "best of model_2..5"]]
desc = {"6TFF": "16, all correct; no pseudoknot", "6WJR": "34; 8 as pseudoknot, 2 real", "6Q57": "27; full 5-pair pseudoknot"}
for t in TGT:
    m = pipe[t]["models"]; m1 = m[0]; alt = min(m[1:], key=lambda d: d["RMSD"])
    rows.append([t, desc[t], f"{m1['cluster_pct']:.0f}%", f2(m1["TM"]), f1(m1["RMSD"]), f1(m1["RMSD_core"]), f2(m1["INF"]),
                 f"{alt['model'][:-4]}: TM {f2(alt['TM'])}, {f1(alt['RMSD'])} Å"])
story += [tbl(rows, [0.5 * inch, 1.7 * inch, 0.75 * inch, 0.45 * inch, 0.6 * inch, 0.75 * inch, 0.45 * inch, 1.8 * inch]),
          Spacer(1, 4),
          Paragraph("Predicted pairs (pairs.txt):", S),
          Paragraph("6TFF " + pipe["6TFF"]["pairs"].replace("<", "&lt;").replace(">", "&gt;"), M),
          Paragraph("6Q57 " + pipe["6Q57"]["pairs"], M),
          Paragraph("6WJR " + pipe["6WJR"]["pairs"], M),
          Paragraph("TM-score above 0.45 means the same overall fold. A large first cluster means the runs agree, not that the model "
                    "is right: on 6Q57, 96% of the low-energy shapes agree on a structure 13.8 Å from the crystal. With fixed random seeds the pipeline reproduces the snapshots and clusters of the validation runs below; the final all-atom rebuild adds a small random variation (6WJR model_1: 21.9 Å here, 21.7 Å in the validation tables).", S),
          PageBreak(),
          Paragraph("Validation: is RibonanzaNet2 the right source of pairs?", H2),
          Paragraph("The same SimRNA run was repeated with four sources of base pairs. Only the pairs differ, so differences come from "
                    "the pairs. EternaFold represents the pairing step RNAnneal uses (no pseudoknots); crystal pairs are a ceiling "
                    "control and the only place the crystal was ever given to SimRNA.", B)]
rows = [["Arm", "6TFF (20 pairs, 2 pk)", "6WJR (32, 4 pk)", "6Q57 (30, 5 pk)"]]
for a in ARMS:
    rows.append([LAB[a]] + ["-" if a == "nocon" else f"{pairs[t][a]['correct']}/{pairs[t][a]['n']} correct, pk {pairs[t][a]['pk_found']}" for t in TGT])
story += [tbl(rows, [1.3 * inch, 1.9 * inch, 1.9 * inch, 1.9 * inch], hl=[(1, 3), (2, 3), (3, 3)]), Spacer(1, 6),
          fig("fig1_arms.png", "<b>Figure 1.</b> Largest-cluster model for each RNA and source of pairs. Left: TM-score (dashed line, "
              "0.45 = same fold). Middle: whole-chain C1' RMSD. Right: INF, agreement between the model's and the crystal's canonical base pairs, sqrt(precision x recall).", 2.5 * inch)]
rows = [["Mean over 3 RNAs", "no pairs", "EternaFold", "RibonanzaNet2", "crystal pairs"]]
for m, lab, f in [("TM", "TM-score", f2), ("RMSD", "RMSD (Å)", f1), ("INF", "INF", f2)]:
    rows.append([lab] + [f(np.mean([top.loc[(t, a), m] for t in TGT])) for a in ARMS])
pk = lambda a: ", ".join(f"{t} {int(top.loc[(t, a), 'pk_formed'])}/{int(top.loc[(t, a), 'pk_ref'])}" for t in TGT)
rows.append(["Pseudoknot pairs formed"] + [pk(a) for a in ARMS])
story += [tbl(rows, [1.4 * inch, 1.3 * inch, 1.3 * inch, 1.5 * inch, 1.5 * inch], hl=[(3, 1), (3, 3)]), Spacer(1, 4),
          Paragraph("RibonanzaNet2 pairs give the highest INF and TM of the predicted sources and the only predicted 6Q57 "
                    "pseudoknot, but not a lower RMSD than EternaFold (16.4 vs 15.2 Å). Every model is below TM 0.45.", B),
          Paragraph("Every snapshot: sampling or selection?", H2),
          fig("fig2_snapshots.png", "<b>Figure 2.</b> C4' RMSD to the crystal of all 32,080 snapshots per run (violins). Triangles: "
              "closest snapshot ever produced. Circles: the chosen model. Dashed line: 5 Å.", 2.6 * inch)]
rows = [["RNA", "no pairs", "EternaFold", "RibonanzaNet2", "crystal pairs", "% < 5 Å (RNet2)", "% < 15 Å (RNet2)", "RNet2 chosen model"]]
pct = lambda t, a, x: 100 * (snap[(t, a)][1] < x).mean()
for t in TGT:
    rows.append([t] + [f1(fr.loc[(t, a), "best_sampled"]) for a in ARMS] + [f"{pct(t, 'rnet2', 5):.1f}%", f"{pct(t, 'rnet2', 15):.2g}%",
                                                                              f1(top.loc[(t, 'rnet2'), 'RMSD'])])
rows2 = [["< 15 Å · median RMSD"] + [LAB[a] for a in ARMS]]
for t in TGT:
    rows2.append([t] + [f"{pct(t, a, 15):.2g}% · {np.median(snap[(t, a)][1]):.1f} Å" for a in ARMS])
story += [Paragraph("Closest snapshot ever produced (Å):", S),
          tbl(rows, [0.5 * inch, 0.75 * inch, 0.85 * inch, 1.0 * inch, 0.9 * inch, 0.95 * inch, 1.0 * inch, 1.05 * inch],
              hl=[(3, 1), (3, 2), (3, 3), (6, 3)]),
          Spacer(1, 4),
          Paragraph("Share of all snapshots within 15 Å of the crystal, and their median RMSD:", S),
          tbl(rows2, [1.4 * inch] * 5, hl=[(3, 2), (3, 3)]),
          Spacer(1, 6),
          fig("fig3_energy.png", "<b>Figure 3.</b> SimRNA energy against RMSD for every snapshot of the RibonanzaNet2 run. Blue: "
              "the lowest-energy 1%, from which the answer is chosen. The energy barely tracks correctness (rank correlation in the "
              "titles), so even when a near-crystal shape exists (6TFF, 3.1 Å) it is not chosen.", 2.5 * inch),
          Paragraph("Only 6TFF ever got within 5 Å; for 6WJR and 6Q57 the right shape was never generated. RibonanzaNet2 pairs give "
                    "the closest snapshot of any predicted source on all three RNAs, closer even than the crystal's own pairs on 6WJR "
                    "and 6Q57. On 6Q57, 40% of its snapshots are within 15 Å, against 5% with EternaFold: the pseudoknot "
                    "pulls the whole search toward the right region, almost as well as the crystal's pairs (49%).", B),
          PageBreak(),
          Paragraph("Comparison with RNAnneal", H2),
          Paragraph("RNAnneal (Herron, ..., Tiwary, bioRxiv 2026.02.01.703098) benchmarked 16 riboswitches including these three. "
                    "6TFF values are from the table in its Fig. 3e (exact); 6WJR and 6Q57 values were read off its Figs. 2e-h and 3a "
                    "(approximate, about ±1 Å and ±0.02).", B),
          tbl([["Step", "RNAnneal", "RibonanzaNet2 + SimRNA"],
               ["Base pairs", "EternaFold, LinearFold, LinearAliFold (no pseudoknots)", "RibonanzaNet2 head L25H1 (pseudoknots kept)"],
               ["3D building", "FARFAR2, then 100 ns MD on 200 models", "SimRNA replica-exchange Monte Carlo"],
               ["Candidates per RNA", "20,000", "32,080"],
               ["Choosing", "RNAnneal score (Thermodynamic Maps, learned per RNA); 10 ranked states", "largest cluster of the lowest-energy 1%; 5 models"],
               ["Compute per RNA", "about 20 µs of MD on GPUs (estimate from the methods)", "about 15-55 CPU-core hours, under an hour wall time"]],
              [1.3 * inch, 2.85 * inch, 2.85 * inch]),
          Spacer(1, 6),
          tbl([["RNA", "RNAnneal: top-ranked", "RNAnneal: closest", "RibonanzaNet2 + SimRNA: model_1", "RibonanzaNet2 + SimRNA: closest snapshot", "Pseudoknot built"],
               ["6TFF", "14.0 Å, INF 0.63", "2.29 Å, INF 0.90 (ranked 5th)", f"{f1(ours_top['6TFF'])} Å, INF {f2(top.loc[('6TFF','rnet2'),'INF'])}", f"{f1(ours_best['6TFF'])} Å", "neither"],
               ["6WJR", "unclear in the preprint", "~15 Å, INF <= 0.72", f"{f1(ours_top['6WJR'])} Å, INF {f2(top.loc[('6WJR','rnet2'),'INF'])}", f"{f1(ours_best['6WJR'])} Å", "neither"],
               ["6Q57", "> 10 Å, INF < 0.7", "~7.5 Å, INF <= 0.82", f"{f1(ours_top['6Q57'])} Å, INF {f2(top.loc[('6Q57','rnet2'),'INF'])}", f"{f1(ours_best['6Q57'])} Å", "RibonanzaNet2 + SimRNA, 5 of 5 pairs"]],
              [0.5 * inch, 1.25 * inch, 1.45 * inch, 1.25 * inch, 1.25 * inch, 1.3 * inch], hl=[(4, 1), (4, 2), (4, 3), (5, 3)]),
          Spacer(1, 6),
          fig("fig4_rnanneal.png", "<b>Figure 4.</b> RMSD to the crystal of the chosen answer (left) and of the closest shape generated "
              "(right). RNAnneal values from the preprint; its 6Q57 answer is only reported as above 10 Å (bar at 10), and its 6WJR "
              "answer is unclear (text and figures disagree). Dashed line: 5 Å.", 2.6 * inch),
          tbl([["Choosing the right shape", "6TFF", "6WJR", "6Q57"],
               ["RNAnneal score, BEDROC", "~0.80", "~0.71", "~0.90"],
               ["Rosetta score, BEDROC", "~0.16", "~0.14", "~0.57"],
               ["DES-AMBER force field, BEDROC", "~0.20", "~0.54", "~0.55"],
               ["SimRNA energy vs RMSD, rank correlation (RibonanzaNet2 + SimRNA)"] + [f2(fr.loc[(t, 'rnet2'), 'spearman_E_rmsd']) for t in TGT]],
              [3.0 * inch, 1.3 * inch, 1.3 * inch, 1.3 * inch], hl=[(1, 1), (2, 1), (3, 1)]),
          Spacer(1, 4),
          Paragraph("<b>Reading the comparison.</b> On 6TFF the two pipelines behave the same way: a wrong shape first (14.0 vs 13.8 Å) "
                    "and a near-crystal shape only as a minor alternative (2.3 vs 3.1 Å). On 6WJR and 6Q57 RibonanzaNet2 + SimRNA reaches a little "
                    "closer (13.6 vs ~15 Å, 6.5 vs ~7.5 Å), and only RibonanzaNet2 + SimRNA forms 6Q57's pseudoknot: that is the RibonanzaNet2 pairs at "
                    "work, at a fraction of the compute. RNAnneal's score ranks near-crystal shapes far better than Rosetta, the force "
                    "field or SimRNA's energy. The natural next step is RibonanzaNet2 pairs to generate candidates and RNAnneal's "
                    "score to choose among them.", B),
          Paragraph("Caveats", H2),
          Paragraph("SimRNA runs were 10% of its published benchmark length. TM-score is harsh for RNAs of 50-110 nt. RNAnneal's RMSD "
                    "is probably over all atoms, RibonanzaNet2 + SimRNA is scored over one backbone atom per nucleotide (C1' for models, C4' for snapshots), so "
                    "differences of a few Å are not meaningful. RNAnneal's 6WJR and 6Q57 numbers are read off plots, and its text "
                    "says 6WJR had a representative within 5 Å while its Fig. 2e shows no candidate closer than about 15 Å. "
                    "A near-identical copy of 6Q57 (4LVW, 99%) was among the 612 RNAs used to choose RibonanzaNet2's attention head; "
                    "it affected only which of 576 heads was read, not the model itself.", B),
          Paragraph("Files", H2),
          Paragraph("Pipeline: scripts/run_slurm.sh, scripts/run_local.sh, rnet2_simrna/. Validation runs: analysis/run_validation.sh, "
                    "analysis/collect_validation.py. Results: results/validation/{summary.csv, scores.csv, frames_rmsd.csv, snapshots.csv.gz}, "
                    "results/pipeline/{6TFF,6WJR,6Q57}/. This report: analysis/make_report.py.", S)]


def footer(c, d):
    c.saveState(); c.setFont("Helvetica", 7); c.setFillColor(colors.HexColor("#5b6a65"))
    c.drawString(0.75 * inch, 0.5 * inch, "RibonanzaNet2 + SimRNA protocol report")
    c.drawRightString(7.75 * inch, 0.5 * inch, f"page {d.page}"); c.restoreState()


SimpleDocTemplate(str(OUT / "report.pdf"), pagesize=letter, leftMargin=0.75 * inch, rightMargin=0.75 * inch,
                  topMargin=0.7 * inch, bottomMargin=0.75 * inch, title="RibonanzaNet2 + SimRNA report").build(
    story, onFirstPage=footer, onLaterPages=footer)
print("wrote", OUT / "report.pdf")
