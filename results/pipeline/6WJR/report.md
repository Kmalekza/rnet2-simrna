# 6WJR: RibonanzaNet2 + SimRNA

Sequence (112 nt): `GGAUCUUCGGGGCAGGGUGAAAUUCCCGACCGGUGGUAUAGUCCACGAAAGUAUUUGCUUUGAUUUGGUGAAAUUCCAAAACCGACAGUAGAGUCUGGAUGAGAGAAGAUUC`

Predicted base pairs (RibonanzaNet2; [ ] = pseudoknot):

```
GGAUCUUCGGGGCAGGGUGAAAUUCCCGACCGGUGGUAUAGUCCACGAAAGUAUUUGCUUUGAUUUGGUGAAAUUCCAAAACCGACAGUAGAGUCUGGAUGAGAGAAGAUUC
(((((((((.....(((.......))).[[).((((.....{)))).(((...[[[[[)))(.(((((.......)))))..)(...]]]]]]])...}.....))))))))
```

34 pairs, 8 in pseudoknots. SimRNA: 8 runs x 10 replicas x 1,600,000 steps = 32,080 snapshots; lowest-energy 1% clustered at 11.2 A.

**Prediction: `model_1.pdb`** (centre of the largest cluster). `model_2..5.pdb` are the next clusters.

Scored against `data/structures/6WJR.cif` chain X.

| model | cluster size | TM-score | RMSD (A) | core RMSD (A) | INF |
|---|---|---|---|---|---|
| model_1.pdb | 58% | 0.189 | 21.89 | 5.18 | 0.574 |
| model_2.pdb | 25% | 0.218 | 19.05 | 3.97 | 0.568 |
| model_3.pdb | 6% | 0.204 | 24.25 | 4.40 | 0.625 |
| model_4.pdb | 6% | 0.198 | 22.93 | 3.51 | 0.364 |
| model_5.pdb | 2% | 0.219 | 20.57 | 4.93 | 0.552 |

TM-score: 0-1, above 0.45 means the same overall fold. RMSD: whole chain, C1', residue by residue. Core RMSD: only the residues US-align could superpose. INF: share of base pairs in common with the reference (0-1).
