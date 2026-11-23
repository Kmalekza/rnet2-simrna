# 6TFF: RibonanzaNet2 + SimRNA

Sequence (52 nt): `GGCUUCAACAACCCCGUAGGUUGGGCCGAAAGGCAGCGAAUCUACUGGAGCC`

Predicted base pairs (RibonanzaNet2; [ ] = pseudoknot):

```
GGCUUCAACAACCCCGUAGGUUGGGCCGAAAGGCAGCGAAUCUACUGGAGCC
((((((.........(((((((..(((....)))....))))))).))))))
```

16 pairs, 0 in pseudoknots. SimRNA: 8 runs x 10 replicas x 1,600,000 steps = 32,080 snapshots; lowest-energy 1% clustered at 5.2 A.

**Prediction: `model_1.pdb`** (centre of the largest cluster). `model_2..5.pdb` are the next clusters.

Scored against `data/structures/6TFF.cif` chain A.

| model | cluster size | TM-score | RMSD (A) | core RMSD (A) | INF |
|---|---|---|---|---|---|
| model_1.pdb | 63% | 0.256 | 13.80 | 3.77 | 0.922 |
| model_2.pdb | 16% | 0.342 | 9.82 | 2.96 | 0.896 |
| model_3.pdb | 7% | 0.269 | 13.90 | 3.17 | 0.949 |
| model_4.pdb | 2% | 0.273 | 8.39 | 3.77 | 0.949 |
| model_5.pdb | 2% | 0.376 | 6.41 | 3.34 | 0.922 |

TM-score: 0-1, above 0.45 means the same overall fold. RMSD: whole chain, C1', residue by residue. Core RMSD: only the residues US-align could superpose. INF: share of base pairs in common with the reference (0-1).
