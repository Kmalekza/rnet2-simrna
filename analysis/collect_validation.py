#!/usr/bin/env python3
"""Collect the validation runs (analysis/run_validation.sh) into the tables used by analysis/make_report.py.

  python analysis/collect_validation.py [RUNS_DIR]     (default results/validation_runs)
Writes results/validation/{scores.csv, summary.csv, frames_rmsd.csv, snapshots.csv.gz}:
  scores.csv       per target x arm x cluster: TM, whole-chain RMSD, core RMSD, INF, pseudoknot pairs formed
  summary.csv      per target x arm: the largest-cluster model and the best of the five
  snapshots.csv.gz every SimRNA snapshot: energy and C4' RMSD to the crystal (residue matched, Kabsch superposition)
  frames_rmsd.csv  per target x arm: closest snapshot, share within 5 / 8 A, best in the lowest-energy 1%, energy-RMSD rank correlation
"""
import json
import sys
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from rnet2_simrna.score import wc_pairs  # noqa: E402

RUNS = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "results" / "validation_runs"
OUT = REPO / "results" / "validation"; OUT.mkdir(parents=True, exist_ok=True)
TG = pd.read_csv(REPO / "data" / "targets.csv").set_index("pdb_id")
REF = json.load(open(REPO / "data" / "reference_pairs.json"))
ARMS = ["nocon", "eterna", "rnet2", "xtal"]


def crystal_c4(t):
    st = gemmi.read_structure(str(REPO / TG.loc[t].structure)); st.setup_entities(); out = {}
    for res in st[0][TG.loc[t].chain].get_polymer():
        a = res.find_atom("C4'", "*")
        if a is not None:
            out[res.label_seq - 1] = a.pos.tolist()
    return out


def rmsd_many(X, Y):
    X = X - X.mean(1, keepdims=True); Y = Y - Y.mean(0)
    H = np.einsum("fni,nj->fij", X, Y); U, S, Vt = np.linalg.svd(H)
    S[:, -1] *= np.sign(np.linalg.det(np.einsum("fij,fjk->fik", U, Vt)))
    return np.sqrt(np.maximum((X ** 2).sum((1, 2)) + (Y ** 2).sum() - 2 * S.sum(1), 0) / X.shape[1])


scores, snaps, frames = [], [], []
for t in TG.index:
    seq = REF[t]["seq"]; L = len(seq); pk = {tuple(p) for p in REF[t]["pk"]}
    ref = crystal_c4(t); idx = sorted(i for i in ref if i < L); Y = np.array([ref[i] for i in idx])
    for a in ARMS:
        d = RUNS / f"{t}_{a}"; rep = json.load(open(d / "report.json"))
        for k, m in enumerate(rep["models"], 1):
            pm = wc_pairs(d / m["model"], seq, "A")
            scores.append(dict(target=t, arm=a, cluster=k, cluster_pct=m["cluster_pct"], TM=m["TM"], RMSD=m["RMSD"],
                               RMSD_core=m["RMSD_core"], INF=m["INF"], pk_formed=len(pm & pk), pk_ref=len(pk), pairs_model=len(pm)))
        E, Et, C = [], [], []
        with open(d / "runs" / "all.trafl") as fh:
            for h in fh:
                c = fh.readline(); p = h.split(); Et.append(float(p[2])); E.append(float(p[3]))
                C.append(np.array(c.split(), dtype=float).reshape(L, 5, 3)[idx, 1])
        r = rmsd_many(np.array(C), Y); E = np.array(E); low = np.argsort(E)[: max(1, len(E) // 100)]
        snaps.append(pd.DataFrame(dict(target=t, arm=a, frame=np.arange(len(r)), energy=np.round(E, 4),
                                       energy_with_restraints=np.round(Et, 4), rmsd_c4=np.round(r, 6))))
        frames.append(dict(target=t, arm=a, n_frames=len(r), n_res=len(idx), best_sampled=r.min(), pct_below_5A=100 * (r < 5).mean(),
                           pct_below_8A=100 * (r < 8).mean(), best_in_lowE1pct=r[low].min(), median_lowE1pct=np.median(r[low]),
                           spearman_E_rmsd=pd.Series(E).corr(pd.Series(r), method="spearman")))
        print(t, a, "done", flush=True)

sc = pd.DataFrame(scores); sc.to_csv(OUT / "scores.csv", index=False, float_format="%.4f")
top = sc[sc.cluster == 1].set_index(["target", "arm"]); g = sc.groupby(["target", "arm"])
pd.DataFrame({"TM_top1": top.TM, "TM_best5": g.TM.max(), "RMSD_top1": top.RMSD, "RMSD_best5": g.RMSD.min(), "INF_top1": top.INF,
              "pk_top1": top.pk_formed.astype(str) + "/" + top.pk_ref.astype(str), "cluster1_pct": top.cluster_pct}
             ).loc[top.index].reset_index().to_csv(OUT / "summary.csv", index=False, float_format="%.3f")
pd.DataFrame(frames).to_csv(OUT / "frames_rmsd.csv", index=False, float_format="%.3f")
pd.concat(snaps).to_csv(OUT / "snapshots.csv.gz", index=False, compression="gzip")
print("wrote", OUT)
