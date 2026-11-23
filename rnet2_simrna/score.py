# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Karim Malekzadeh
"""Scoring a model against a reference structure: TM-score and RMSD (US-align, C1'), INF (rnapolis canonical pairs).

External programs, found on PATH or through environment variables:
  USALIGN    US-align binary            (default: USalign)
  ANNOTATOR  rnapolis 'annotator' CLI   (default: annotator; pip install rnapolis)
"""
from __future__ import annotations

import csv
import os
import subprocess
import tempfile
from pathlib import Path

import gemmi
import numpy as np

USALIGN = os.environ.get("USALIGN", "USalign")
ANNOTATOR = os.environ.get("ANNOTATOR", "annotator")
CANON = {"AU", "UA", "GC", "CG", "GU", "UG"}
NONSTD = {"GTP": "G", "GDP": "G", "ATP": "A", "CTP": "C", "UTP": "U"}


def _chain(st, chain):
    return st[0][chain] if chain and st[0].find_chain(chain) else st[0][0]


def ref_sequence(path, chain=None) -> tuple[str, str]:
    """Full deposited sequence (mmCIF entity) when available, otherwise the modelled residues."""
    st = gemmi.read_structure(str(path)); st.setup_entities(); ch = _chain(st, chain)
    ent = st.get_entity_of(ch.get_polymer())
    names = [n.split(",")[0] for n in ent.full_sequence] if ent is not None and ent.full_sequence else [r.name for r in ch.get_polymer()]
    return "".join(n if n in "ACGU" else NONSTD.get(n, "X") for n in names), ch.name


def res_index(path, chain=None) -> dict:
    """(chain, residue number) -> 0-based sequence index; mmCIF label_seq when present, else residue order."""
    st = gemmi.read_structure(str(path)); st.setup_entities(); ch = _chain(st, chain); out = {}
    for k, r in enumerate(ch.get_polymer()):
        out[(ch.name, str(r.seqid.num) + (r.seqid.icode.strip() or ""))] = r.label_seq - 1 if r.label_seq is not None else k
    return out


def clean_pdb(path, out, chain=None, strip_o2=False) -> Path:
    """Polymer chain only, no ligands, waters or hydrogens, written as PDB."""
    st = gemmi.read_structure(str(path)); st.setup_entities(); st.remove_ligands_and_waters(); st.remove_hydrogens()
    if chain and st[0].find_chain(chain):
        for ch in list(st[0]):
            if ch.name != chain:
                st[0].remove_chain(ch.name)
    if strip_o2:
        for ch in st[0]:
            for r in ch:
                for k in reversed(range(len(r))):
                    if r[k].name == "O2'":
                        del r[k]
    st.write_pdb(str(out)); return Path(out)


def _key(nt: str):
    ch, rest = nt.split(".", 1)
    num = "".join(c for c in rest.split("^")[0] if c.isdigit() or c == "-")
    return ch, num + (rest.split("^")[1] if "^" in rest else "")


def wc_pairs(path, seq: str, chain=None) -> set[tuple[int, int]]:
    """Canonical cWW pairs (AU/GC/GU) annotated by rnapolis, as 0-based sequence indices."""
    idx = res_index(path, chain)
    with tempfile.TemporaryDirectory() as d:
        pdb = clean_pdb(path, Path(d) / "s.pdb", chain, strip_o2=True); out = Path(d) / "bp.csv"
        subprocess.run([ANNOTATOR, "-f", "-c", str(out), str(pdb)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        res = set()
        for row in csv.DictReader(open(out)):
            if row["type"] != "base pair" or row["classification-1"] != "cWW":
                continue
            a, b = _key(row["nt1"]), _key(row["nt2"])
            if a in idx and b in idx:
                i, j = sorted((idx[a], idx[b]))
                if seq[i] + seq[j] in CANON:
                    res.add((i, j))
    return res


def _usalign(model, ref, extra):
    o = subprocess.run([USALIGN, str(model), str(ref), "-atom", " C1'", "-mol", "RNA", "-outfmt", "2"] + extra,
                       capture_output=True, text=True, check=True).stdout
    return [l for l in o.splitlines() if not l.startswith("#")][0].split("\t")


def score(model, ref, chain, seq) -> dict:
    """TM (normalised by the reference), RMSD over the whole chain residue by residue, core RMSD, INF."""
    with tempfile.TemporaryDirectory() as d:
        rp = clean_pdb(ref, Path(d) / "ref.pdb", chain)
        f = _usalign(model, rp, []); g = _usalign(model, rp, ["-TMscore", "7"])
    pm = wc_pairs(model, seq, "A"); pr = wc_pairs(ref, seq, chain); tp = len(pm & pr)
    inf = float(np.sqrt(tp / len(pm) * tp / len(pr))) if tp else 0.0
    return dict(TM=float(f[3]), RMSD=float(g[4]), RMSD_core=float(f[4]), INF=inf,
                pairs_model=len(pm), pairs_ref=len(pr), pairs_shared=tp)
