# 6Q57: RibonanzaNet2 + SimRNA

Sequence (89 nt): `GGAGAGUAGAUGAUUCGCGUUAAGUGUGUGUGAAUGGGAUGUCGUCACACAACGAAGCGAGAGCGCGGUGAAUCAUUGCAUCCGCUCCA`

Predicted base pairs (RibonanzaNet2; [ ] = pseudoknot):

```
GGAGAGUAGAUGAUUCGCGUUAAGUGUGUGUGAAUGGGAUGUCGUCACACAACGAAGCGAGAGCGCGGUGAAUCAUUGCAUCCGCUCCA
((((.....((((((((((....(...((((((.(.[[[[[..)))))))..)..........)))...)))))))..]]]]].)))).
```

27 pairs, 5 in pseudoknots. SimRNA: 8 runs x 10 replicas x 1,600,000 steps = 32,080 snapshots; lowest-energy 1% clustered at 8.9 A.

**Prediction: `model_1.pdb`** (centre of the largest cluster). `model_2..5.pdb` are the next clusters.

Scored against `data/structures/6Q57.cif` chain A.

| model | cluster size | TM-score | RMSD (A) | core RMSD (A) | INF |
|---|---|---|---|---|---|
| model_1.pdb | 96% | 0.319 | 13.75 | 4.30 | 0.658 |
| model_2.pdb | 2% | 0.246 | 14.75 | 4.60 | 0.590 |
| model_3.pdb | 1% | 0.344 | 11.61 | 4.66 | 0.604 |
| model_4.pdb | 0% | 0.288 | 11.90 | 4.62 | 0.721 |

TM-score: 0-1, above 0.45 means the same overall fold. RMSD: whole chain, C1', residue by residue. Core RMSD: only the residues US-align could superpose. INF: share of base pairs in common with the reference (0-1).
