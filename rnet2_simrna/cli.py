#!/usr/bin/env python3
"""RibonanzaNet2 + SimRNA, step by step. Normally driven by scripts/run_slurm.sh or scripts/run_local.sh.

  python -m rnet2_simrna.cli prep NAME SEQUENCE_OR_FASTA OUTDIR [--ref FILE] [--chain C] [--iters N] [--runs R] [--pairs JSON|none]
      checks the sequence, predicts base pairs with RibonanzaNet2 (or takes them from --pairs), writes the SimRNA inputs
      (OUTDIR/runs/{seq, ss, config.dat}) and OUTDIR/pairs.txt
  python -m rnet2_simrna.cli post OUTDIR
      pools the SimRNA snapshots, keeps the lowest-energy 1%, clusters at 0.1 L A, rebuilds all atoms for the 5 largest
      clusters -> OUTDIR/model_1.pdb (the prediction) .. model_5.pdb; scores them if --ref was given; writes report.md/json

Environment: RN2_DIR (RibonanzaNet2 checkpoint dir, default weights/ribonanzanet2), SIMRNA_DIR (SimRNA folder holding the
SimRNA binary and its data/ directory, default third_party/SimRNA_64bitIntel_Linux), USALIGN, ANNOTATOR (see score.py).
--pairs is for the validation runs only: a JSON list of 0-based [i, j] pairs used instead of RibonanzaNet2, or 'none'.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RN2_DIR = Path(os.environ.get("RN2_DIR", REPO / "weights" / "ribonanzanet2"))
SIMRNA_DIR = Path(os.environ.get("SIMRNA_DIR", REPO / "third_party" / "SimRNA_64bitIntel_Linux"))
BR = ["()", "[]", "{}", "<>"]


def read_sequence(arg: str) -> str:
    if Path(arg).is_file():
        arg = "".join(l.strip() for l in open(arg) if not l.startswith(">"))
    return re.sub(r"\s", "", arg).upper().replace("T", "U")


def prep(a):
    from .lm import layers
    seq = read_sequence(a.sequence); bad = set(seq) - set("ACGU")
    if bad or not 10 <= len(seq) <= 500:
        sys.exit(f"sequence must be 10-500 nt of A/C/G/U (T is read as U); got {len(seq)} nt, other letters: {sorted(bad)}")
    out = Path(a.outdir); (out / "runs").mkdir(parents=True, exist_ok=True); L = len(seq)
    if a.pairs is None:
        import torch
        from .lm import load_model, predict_pairs
        dev = "cuda" if torch.cuda.is_available() else "cpu"
        pairs = predict_pairs(load_model(RN2_DIR, dev), seq, dev); source = "RibonanzaNet2 L25H1"
    elif a.pairs == "none":
        pairs, source = [], "no pairs"
    else:
        pairs = sorted(tuple(p) for p in json.load(open(a.pairs))); source = f"pairs from {Path(a.pairs).name}"
    lay = layers(pairs, L); ss_lines, db = [], ["."] * L
    for k, ly in enumerate(lay):
        line = ["."] * L
        for i, j in ly:
            line[i], line[j] = "(", ")"; db[i], db[j] = BR[min(k, 3)]
        ss_lines.append("".join(line))
    (out / "runs" / "seq").write_text(seq + "\n")
    if ss_lines:
        (out / "runs" / "ss").write_text("\n".join(ss_lines) + "\n")
    cfg = (SIMRNA_DIR / "config.dat").read_text()
    cfg = re.sub(r"NUMBER_OF_ITERATIONS \d+", f"NUMBER_OF_ITERATIONS {a.iters}", cfg)
    cfg = re.sub(r"TRA_WRITE_IN_EVERY_N_ITERATIONS \d+", f"TRA_WRITE_IN_EVERY_N_ITERATIONS {max(1, a.iters // 400)}", cfg)
    (out / "runs" / "config.dat").write_text(cfg)
    n_pk = sum(len(x) for x in lay[1:])
    (out / "pairs.txt").write_text(f">{a.name} {source}, {len(pairs)} pairs, {n_pk} in pseudoknots\n{seq}\n{''.join(db)}\n")
    json.dump(dict(name=a.name, seq=seq, L=L, source=source, pairs=[list(p) for p in pairs], n_pk=n_pk,
                   ref=str(Path(a.ref).resolve()) if a.ref else None, chain=a.chain, iters=a.iters, runs=a.runs),
              open(out / "meta.json", "w"), indent=1)
    print(open(out / "pairs.txt").read())


def sh(cmd, cwd):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)


def post(a):
    out = Path(a.outdir).resolve(); run = out / "runs"; meta = json.load(open(out / "meta.json"))
    L = meta["L"]; thr = f"{0.1 * L:.1f}"
    trafl = sorted(run.glob("r??_??.trafl"))
    if not trafl:
        sys.exit("no SimRNA trajectories found in " + str(run))
    with open(run / "all.trafl", "w") as fo:
        for f in trafl:
            fo.write(f.read_text())
    for f in run.glob("all_thrs*"):
        f.unlink()
    if not (run / "data").exists():
        (run / "data").symlink_to(SIMRNA_DIR / "data")
    r = sh([str(SIMRNA_DIR / "clustering"), "all.trafl", "0.01", thr], run); log = r.stdout + r.stderr
    (run / "clustering.log").write_text(log)
    size = {int(k): float(p) for k, p in re.findall(r"^cluster (\d+) contains \d+ \(([\d.]+)%\)", log, re.M)}
    models = []
    for k in range(1, 6):
        c = sorted(run.glob(f"all_thrs*A_clust{k:02d}.trafl"))
        if not c:
            continue
        sh([str(SIMRNA_DIR / "SimRNA_trafl2pdbs"), "r01_01-000001.pdb", c[0].name, "1", "AA"], run)
        aa = run / (c[0].stem + "-000001_AA.pdb")
        (out / f"model_{k}.pdb").write_text(aa.read_text()); aa.unlink()
        models.append(dict(model=f"model_{k}.pdb", cluster_pct=size.get(k)))
    rep = dict(name=meta["name"], length=L, sequence=meta["seq"], pairs_source=meta.get("source"),
               pairs=open(out / "pairs.txt").read().split("\n")[2], n_pairs=len(meta["pairs"]), n_pk_pairs=meta["n_pk"],
               snapshots=sum(1 for _ in open(run / "all.trafl")) // 2, models=models, reference=meta["ref"])
    note = ""
    if meta["ref"]:
        from .score import ref_sequence, score
        rs, rch = ref_sequence(meta["ref"], meta["chain"])
        if rs != meta["seq"]:
            note = (f"Reference chain {rch} sequence ({len(rs)} nt) differs from the input ({L} nt); TM and RMSD use "
                    "US-align's alignment, INF uses residue numbering and may be off.")
        for m in models:
            m.update(score(out / m["model"], meta["ref"], meta["chain"], meta["seq"]))
    json.dump(rep, open(out / "report.json", "w"), indent=1)
    lines = [f"# {meta['name']}: RibonanzaNet2 + SimRNA", "", f"Sequence ({L} nt): `{meta['seq']}`", "",
             f"Base pairs ({meta.get('source')}; [ ] = pseudoknot):", "", "```", meta["seq"], rep["pairs"], "```", "",
             f"{rep['n_pairs']} pairs, {rep['n_pk_pairs']} in pseudoknots. SimRNA: {meta['runs']} runs x 10 replicas x "
             f"{meta['iters']:,} steps = {rep['snapshots']:,} snapshots; lowest-energy 1% clustered at {thr} A.", "",
             "**Prediction: `model_1.pdb`** (centre of the largest cluster). `model_2..5.pdb` are the next clusters.", ""]
    if meta["ref"]:
        lines += [f"Scored against `{Path(meta['ref']).name}`" + (f" chain {meta['chain']}" if meta["chain"] else "") + ".", "",
                  "| model | cluster size | TM-score | RMSD (A) | core RMSD (A) | INF |", "|---|---|---|---|---|---|"]
        lines += [f"| {m['model']} | {m['cluster_pct']:.0f}% | {m['TM']:.3f} | {m['RMSD']:.2f} | {m['RMSD_core']:.2f} | {m['INF']:.3f} |"
                  for m in models]
        lines += ["", "TM-score: 0-1, above 0.45 means the same overall fold. RMSD: whole chain, C1', residue by residue. "
                  "Core RMSD: only the residues US-align could superpose. INF: agreement of canonical base pairs with the "
                  "reference, sqrt(precision x recall)."] + (["", f"Note: {note}"] if note else [])
    else:
        lines += ["| model | cluster size |", "|---|---|"] + [f"| {m['model']} | {m['cluster_pct']:.0f}% |" for m in models]
        lines += ["", "No reference structure given, so TM-score, RMSD and INF were not computed. A large first cluster "
                  "means the runs agree; it is not a guarantee of accuracy."]
    (out / "report.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("prep"); p.add_argument("name"); p.add_argument("sequence"); p.add_argument("outdir")
    p.add_argument("--ref"); p.add_argument("--chain"); p.add_argument("--iters", type=int, default=1_600_000)
    p.add_argument("--runs", type=int, default=8); p.add_argument("--pairs")
    q = sp.add_parser("post"); q.add_argument("outdir")
    a = ap.parse_args()
    {"prep": prep, "post": post}[a.cmd](a)


if __name__ == "__main__":
    main()
