# RibonanzaNet2 + SimRNA

Sequence-only RNA 3D structure prediction that combines an RNA language model with a physics-style folding program:

1. **Base pairs.** [RibonanzaNet2](https://www.kaggle.com/models/shujun717/ribonanzanet2), used frozen, reads the sequence. One attention head (layer 25, head 1) is read as a base-pair map and decoded into pairs, including pseudoknots.
2. **3D folding.** [SimRNA 3.20](https://genesilico.pl/software/stand-alone/simrna) folds a coarse-grained chain by replica-exchange Monte Carlo, with those pairs as restraints.
3. **Choosing.** The lowest-energy 1% of snapshots are clustered, and the centre of the largest cluster is the prediction. All atoms are then rebuilt.

The only input is the sequence. An experimental structure is optional and is used only to score the result.

The protocol was tested on three riboswitches from the RNAnneal benchmark: **6TFF** (NAD⁺), **6WJR** (FMN, apo) and **6Q57** (THF). They are compared with RNAnneal (Herron, …, Tiwary, bioRxiv 2026.02.01.703098). The full write-up is in [`results/report.pdf`](results/report.pdf).

## Results in brief

| | 6TFF | 6WJR | 6Q57 |
|---|---|---|---|
| RibonanzaNet2 pairs: correct / predicted, pseudoknot pairs found | 16/16, 0 of 2 | 23/34, 2 of 4 | 22/27, **5 of 5** |
| `model_1.pdb`: TM-score · whole-chain RMSD · INF | 0.26 · 13.8 Å · 0.92 | 0.19 · 21.9 Å · 0.57 | 0.32 · 13.8 Å · 0.66 |
| Closest of 32,080 snapshots | **3.1 Å** | 13.6 Å | 6.5 Å |
| RNAnneal: top-ranked / closest candidate (from the preprint) | 14.0 Å / 2.3 Å | unclear / ~15 Å | >10 Å / ~7.5 Å |

* The language model's pairs reach the 3D models. They give higher base-pair agreement than EternaFold pairs, and they are the only predicted source that builds 6Q57's pseudoknot. On 6Q57, 40% of snapshots are within 15 Å of the crystal, against 5% with EternaFold pairs.
* No model has the right overall fold (TM < 0.45). Even the crystal's own pairs leave SimRNA about 12 Å off. The limit is the 3D arrangement of the helices, not the pairs.
* SimRNA's energy barely ranks shapes by correctness, so a near-crystal snapshot on 6TFF (3.1 Å) is not chosen. RNAnneal's learned score is much better at this step. RibonanzaNet2 pairs to generate candidates plus RNAnneal's score to choose among them is a natural combination.

## Install

Requirements: Linux, Python ≥ 3.10, and:

| Package | Tested version | Purpose |
|---|---|---|
| torch | 2.7.1 | RibonanzaNet2 (CPU is enough; a GPU is used if present) |
| numpy | 1.26.4 | |
| pandas | 2.3.1 | |
| einops | 0.8.2 | RibonanzaNet2 |
| pyyaml | 6.0.2 | RibonanzaNet2 config |
| gemmi | 0.6.7 | structure files |
| rnapolis | 0.17.0 | base-pair annotation; provides the `annotator` command |
| matplotlib | 3.10.5 | figures (analysis only) |
| reportlab | 4.x | PDF report (analysis only) |

```bash
pip install -r requirements.txt
```

Two components can't be redistributed here, so you download them yourself:

* **RibonanzaNet2 weights** (public Kaggle model, about 400 MB): `bash scripts/get_weights.sh` puts them in `weights/ribonanzanet2/`. Set `RN2_DIR` to use another location.
* **SimRNA 3.20** (free for academic use, under its own license): download `SimRNA_64bitIntel_Linux.tgz` from <https://genesilico.pl/software/simrna/> and unpack it into `third_party/`. Set `SIMRNA_DIR` to use another location.
* **US-align** (for scoring): build it from <https://zhanggroup.org/US-align/> and put `USalign` on your `PATH`, or set `USALIGN`.

## Run

```bash
# on a Slurm cluster: 3 jobs (prep ~1 min, 8 folding runs of 10-40 min in parallel, post ~1-5 min)
SLURM_ACCOUNT=myaccount SLURM_PARTITION=batch bash scripts/run_slurm.sh myrna GGAGAGUAGAUGAUUCGCG...

# on one machine (8 runs x 10 threads by default; lower RUNS/THREADS, or ITERS=160000 for a quick test)
bash scripts/run_local.sh myrna myrna.fasta

# optional: score against an experimental structure of the same RNA (path and chain)
bash scripts/run_slurm.sh 6Q57 GGAGAGUAGAUGAUUCG... data/structures/6Q57.cif A
```

Output in `results/pipeline/NAME/`:

| File | Contents |
|---|---|
| `model_1.pdb` | the prediction: all-atom 3D structure, centre of the largest cluster |
| `model_2..5.pdb` | centres of the next four clusters |
| `pairs.txt` | RibonanzaNet2 base pairs in dot-bracket; `[ ]`, `{ }` are pseudoknot layers |
| `report.md`, `report.json` | cluster sizes and, with a reference, TM-score, whole-chain RMSD, core RMSD and INF per model |

Settings (environment variables): `RUNS=8`, `ITERS=1600000` (10% of the SimRNA benchmark length), `OUT=results/pipeline`, `PYTHON=python`.

**Scores**, only when a reference is given:

* **TM-score** is US-align on C1′ atoms, normalised by the reference. Above 0.45 means the same overall fold.
* **RMSD** is C1′, residue by residue, over the whole chain.
* **Core RMSD** covers only the residues US-align could superpose.
* **INF** is √(precision × recall) of canonical base pairs, annotated with rnapolis.

## How the base pairs are read

RibonanzaNet2 has 48 layers × 12 attention heads.

* **Choosing the head.** Head 288 (layer 25, head 1) was chosen once, as the head whose attention best matched known base pairs on 612 training RNAs: SPOT-RNA PDB TR1 plus 500 bpRNA TR0.
* **Benchmark accuracy.** On 185 test RNAs (SPOT-RNA PDB TS1/TS2/TS3, five project targets including these three, and 60 structures released after April 2025), its pair accuracy is within 0.02 F1 of IPknot and EternaFold. On the post-2025 structures it finds about 10× more pseudoknot pairs.
* **Decoding.** The map is made symmetric and corrected for background (average product correction). Pairs are then taken greedily: canonical only (AU, GC, GU), |i − j| ≥ 4, one partner per nucleotide, round(0.3 L) pairs. Crossing pairs are kept as pseudoknots.
* **Caveat.** The training pool includes 4LVW, a 99% copy of 6Q57. It influenced only which of the 576 heads is read, not the model.

## Reproducing the validation and the report

```bash
bash analysis/run_validation.sh             # 12 Slurm runs: 3 RNAs x 4 sources of pairs (none, EternaFold, RibonanzaNet2, crystal)
python analysis/collect_validation.py       # -> results/validation/*.csv, snapshots.csv.gz
python analysis/make_report.py              # -> results/figures/*.png, results/report.pdf
```

`make_report.py` runs straight from the tables shipped in `results/validation/`, so the figures and PDF can be rebuilt without re-running SimRNA. `snapshots.csv.gz` stands in for the 13 GB of trajectories: it holds the energy and C4′ RMSD to the crystal of every snapshot.

The crystal-pairs arm is a ceiling control. It shows the best SimRNA can do with perfect pairs, and it is the only place a crystal structure is given to SimRNA. The protocol itself uses only RibonanzaNet2's pairs.

## Repository layout

```
rnet2_simrna/      lm.py (pairs from RibonanzaNet2), cli.py (prep/post), score.py, ribonanzanet.py (model code, Apache-2.0)
scripts/           run_slurm.sh, run_local.sh, get_weights.sh
analysis/          run_validation.sh, collect_validation.py, make_report.py
data/              targets.csv, structures/*.cif (PDB), reference_pairs.json, validation_pairs_*.json
results/           pipeline/{6TFF,6WJR,6Q57}/, validation/*.csv, figures/, report.pdf
licenses/          Apache-2.0 license of the RibonanzaNet2 model code
```

## Caveats

* SimRNA was run at 10% of its published benchmark length.
* TM-score is harsh for RNAs of 50–110 nt.
* The full-atom rebuild adds a small random variation (6WJR `model_1`: 21.9 Å here, 21.7 Å in the validation tables).
* RNAnneal's numbers for 6WJR and 6Q57 were read off the preprint's figures (about ±1 Å). Its text and Fig. 2e disagree on 6WJR.
* RNAnneal's RMSD is probably over all atoms; here it is one backbone atom per nucleotide.

## Credits

* `rnet2_simrna/ribonanzanet.py` is copied unchanged from NVIDIA's RNAPro (Apache-2.0, © NVIDIA and the AlQuraishi Laboratory; see `licenses/`).
* The RibonanzaNet2 weights are by Shujun He et al. (Kaggle `shujun717/ribonanzanet2`).
* SimRNA: Boniecki et al., *Nucleic Acids Res.* 2016.
* US-align: Zhang et al., *Nature Methods* 2022.
* rnapolis: Zok et al.
* EternaFold: Wayment-Steele et al., *Nature Methods* 2022.
* RNAnneal: Herron et al., bioRxiv 2026.
* Crystal structures are from the PDB (6TFF, 6WJR, 6Q57).

No license has been chosen yet for the original code in this repository; all rights reserved.
