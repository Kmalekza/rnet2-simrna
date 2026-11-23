# Tools and prompts

This is a record of how RibonanzaNet2 + SimRNA was built, as the exercise asks ("report every tool you used including as
many of the prompts as possible"). The work was done in one Claude Code session on 22–23 November 2026 on a Slurm
cluster.

* **Prompts** are quoted exactly as typed, typos included; email addresses are replaced by `[email]`.
* **Time stamps** are UTC.
* **Under each prompt** is a short note of what was done and which files it produced.
* **Numbering** is the order of the prompts in the session; three messages sent while Claude was still working
  are marked "mid-task".

This session continued an earlier project on the same exercise: RNA-FM / structRFM language models with RNAPro 3D
prediction. That work has its own prompt log, which is not part of this repository. Prompts from this session that were
only about that earlier work (RNAPro, RhoFold+, Boltz-2) are listed briefly in section 3, so nothing is hidden.

## 1. Tools

| Tool | Version | Used for |
|---|---|---|
| **Claude Code** (Anthropic), model Claude Opus 5.5 | CLI, Nov 2026 | coding agent for the whole exercise: wrote and ran all scripts, submitted the cluster jobs, analysed results, drafted the pages, the report and this log |
| Claude artifacts (claude.ai) | – | the interactive protocol page during development; its final version is `README.html` / `docs/index.html` |
| **RibonanzaNet2** (Kaggle `shujun717/ribonanzanet2`) | alpha v1 checkpoint | the language model: base pairs from one attention head |
| RibonanzaNet model code (`rnet2_simrna/ribonanzanet.py`, from NVIDIA RNAPro, Apache-2.0) | RNAPro commit 83a0b90 | loading and running RibonanzaNet2 |
| PyTorch | 2.7.1 | running RibonanzaNet2 (GPU for benchmarks, CPU in the pipeline) |
| **SimRNA** (Boniecki et al. 2016), academic license | 3.20 | the tertiary structure model: replica-exchange Monte Carlo folding, clustering, all-atom rebuild |
| EternaFold | source build, commit 87b9aac | comparison source of base pairs (validation only) |
| ViennaRNA RNAfold, LinearFold-V, IPknot | 2.8.0, –, 1.1.0 | classical baselines in the base-pair benchmark used to choose the head |
| US-align | 20241108 | TM-score and RMSD |
| rnapolis | 0.17.0 | base-pair annotation of crystal structures and models (INF) |
| gemmi | 0.6.7 | reading PDB/mmCIF files |
| numpy, pandas, scikit-learn, scipy | 1.26, 2.3, –, – | analysis, logistic-regression probes, Wilcoxon tests |
| matplotlib, reportlab | 3.10, 4.x | figures and `results/report.pdf` |
| pypdf, PyMuPDF | – | reading the RNAnneal preprint PDF, including its figures |
| SPOT-RNA datasets (bpRNA TR0/TS0, PDB TR1/TS1/TS2/TS3), ArchiveII (MultiMolecule mirror), RCSB search API | – | the base-pair benchmark; the post-April-2025 PDB test set |
| REDIAL (`tiwarylab/redial`) | commit ff1648b | mutation-perturbation contact maps, applied to RibonanzaNet2 |
| Slurm (batch and gpu partitions) | – | all heavy runs |
| git, GitHub, GitHub Pages | – | this repository and its web page |

Papers read for the exercise: REDIAL (bioRxiv 2026.05.11.724344) and RNAnneal (bioRxiv 2026.02.01.703098). The REDIAL
paper was read through Claude Code's web-fetch tool. RNAnneal's per-target numbers were read from the PDF, supplied by
the user, after bioRxiv refused automated downloads.

## 2. Prompts that built this protocol

### Choosing the language model

**4 · 2026-11-22 22:19** — "what  RibonanzaNet2 do?"
Explained in plain language: a 48-layer RNA model trained on chemical-probing data from about 2 million RNAs, with a
running per-pair representation.

**5 · 22:47** — "which of the RNA language model are better, RNA-FM, structRFM or RibonanzaNet2?"
Compared the three from earlier results, and proposed testing RibonanzaNet2 the same way.

**6 · 22:52** — "yes, test RibonanzaNet2 the same way"
Ran RibonanzaNet2 through the same base-pair readouts used for the other two models:
* attention averages and the pair representation;
* REDIAL perturbation at the last layer;
* a probe of each attention head, held out one target at a time.

**7 · 23:03** — "do the test for more RNA( and standard test)"
Built a base-pair benchmark.
* **Training pool:** 612 RNAs (SPOT-RNA PDB TR1 plus 500 bpRNA TR0).
* **Test sets:** bpRNA TS0, PDB TS1/TS2/TS3, ArchiveII, and 60 new PDB RNAs released after 2025-04-01 and below 80%
  sequence identity to older entries.
* **Compared:** RNA-FM, structRFM, RibonanzaNet2, RNAfold, LinearFold, EternaFold and IPknot.
* **Head choice:** RibonanzaNet2's best head on the training pool was **layer 25, head 1 (index 288)**. This is the head
  the pipeline reads.
* **Result:** on the new PDB set it finds about 30% of pseudoknot pairs, against at most 4% for the classical tools.

**8 · 23:25** — "you can submit job in the GPU or CPU nodes. they are avaiable"
From then on, jobs were spread across both partitions.

**17 · 2026-11-23 00:41** — "did we used REDIAL in our main protocols or is it possible to use REDIAL to improve RibonanzaNet2?"
Explained where REDIAL had been used so far, and proposed running it on all 48 RibonanzaNet2 layers, including its pair
representation.

**18 · 00:45** — "yes, run REDIAL on RibonanzaNet2"
Ran REDIAL on every layer of RibonanzaNet2, on both its sequence and pair representations.
* **Pair representation, averaged over all layers:** with no training, it tied the trained head and the classical tools
  on 185 test RNAs (Wilcoxon p > 0.05).
* **Sequence representation (the official REDIAL recipe):** much weaker.

**19 · 15:00** — "does it mean that REDIAL increase the accuracy of the RibonanzaNet2?"
No. REDIAL is another way of reading the same model: it matches the trained head without training, but doesn't
exceed it.

### Designing the combination

**20 · 17:04** — "how prof towary used RNA foundation model(Language Model), and tertiary structure prediction model in his two papers?"
Read REDIAL with the web-fetch tool, and RNAnneal from earlier notes.
* **REDIAL:** the language models are the object of study, and RhoFold+ serves as a memorization example.
* **RNAnneal:** no language model; EternaFold-type pairs, then FARFAR2, then MD, then a learned score.

**21 · 17:14** — "i just want to compine two model( one Language Model, and another tertiary structure prediction model) what do you suggest?"
Recommended RibonanzaNet2 pairs passed to a physics-based 3D builder. The reasoning:
* pairs given as restraints must be obeyed, unlike in AlphaFold-style models;
* it fills RNAnneal's pseudoknot gap;
* it can't memorize PDB structures.

FARFAR2 needs a Rosetta license, so SimRNA was offered as the immediately available alternative.

**22 · 17:17** — "start with SimRNA"
Set up SimRNA and the four-way validation.
* **Install:** downloaded SimRNA 3.20 (academic license) and read its manual.
* **Inputs:** wrote `src/simrna_fold.py`, with pseudoknots split into separate dot-bracket layers.
* **Timing test:** 64,000 steps took 76 s on 112 nt, so 1.6 M steps (10% of the benchmark length) take about 30 min.
* **Runs:** submitted 5 RNAs × 4 sources of pairs (none, EternaFold, RibonanzaNet2, crystal) × 8 runs × 10 replicas.

**23 · 17:23** — "how about run them on GPU node. they should be faster?"
SimRNA is CPU-only, but the GPU nodes had 1,272 idle CPU cores. The waiting tasks were moved to `batch,gpu`.
* **Problem:** the account's billing-minutes limit then held them back.
* **Fix:** the time limit was cut from 24 h to 2 h, and 138 of 160 runs started.

**24 · 17:24** — "do we have FARFAR2 here?"
No. Rosetta is not installed, and it needs the user's academic license. The options were explained.

**27 · 17:54** — "check the simulations result?"
Tested clustering and all-atom rebuild on 6TFF, then scored 6TFF and 6Q57.
* **Correction:** US-align's default RMSD covers only the aligned residues, which made it look like 3 Å. The whole-chain
  RMSD (12–27 Å) was added to the scoring.

**28 · 18:02** — "does RNAneal compare it with crystal strcuture and the scoring is base on that?"
No. RNAnneal's score never sees the crystal structure. The crystal is used only to evaluate:
* the near-native label (INF ≥ 0.85, or RMSD < 5 Å);
* AUAC and BEDROC ranking;
* whether a near-native shape was sampled at all.

**29 · 18:09** — "check the rest of the simulation results"
Wrote `src/simrna_frames.py`: the RMSD of all 32,080 snapshots per run.
* **Near-crystal shapes:** sampled only for 6TFF.
* **RibonanzaNet2 pairs:** give the closest snapshot of any predicted source.
* **SimRNA's energy:** barely tracks correctness.

**30 · 18:14** — "is it possible to add MSA and template  steps to this protocol to improve the effiency?"
Possible, but it only helps RNAs with relatives. For 6WJR and 6Q57 the templates are identical-sequence structures, so
the gain would be copying. Not run.

**32 · 18:32** — "what is Eternafold?" · **33 · 18:35** — "why did we use it at all, whne we have RibonanzaNet2?"
EternaFold is the comparison: it is RNAnneal's current pairing step, and it can't predict pseudoknots. It is not part of
the protocol.

### Packaging the protocol (sequence in, structure out)

**35 · 18:38** — "i want ot change this protocol in a way that only get the RNA sequence(like RNAneal), and report the  3d Structure(.pdb) and RM,RMSD, TNF" (sent twice)
Wrote the sequence-only pipeline (the origin of `scripts/run_slurm.sh` and `rnet2_simrna/cli.py`).
* **Steps:** three chained Slurm jobs (prep, fold, post).
* **Output:** `model_1..5.pdb`, `pairs.txt` and `report.md/json`.
* **Scores:** TM-score, RMSD and INF are computed only when an experimental structure is given.
* **Test:** the first test run failed in the report step (clustering output went to stderr); fixed.

**37 · 18:48** — "run the pipeline on 6TFF, 6WJR and 6Q57"
Ran it. Fixed random seeds reproduced the validation clusters. Fixed a misleading sequence-mismatch note for 6WJR, whose
crystal is missing one modelled residue. Results are in `results/pipeline/`.

**38 · 18:49** — "still it shows the crystal structure as inout od SimRNA"
The code was already correct: SimRNA only receives `seq` and `ss`. The page's chart still showed the validation
experiment, so it was redrawn as the sequence-only pipeline, and the crystal-pairs arm was relabelled as a control.

**39 · 18:51** — "what is Training RNAs SPOT-RNA PDB TR1 + 500 bpRNA?"
Explained the head-selection set.
* **Overlap found:** 4LVW, a 99% copy of 6Q57, is in the training pool. This is now stated as a caveat everywhere.

**40 · 19:03** — "did RNAneal predict the structure of these three 6TFF, 6WJR, and 6Q57?"
Yes, all three are in its benchmark. Per-target details were missing from the notes.

**41 · 19:05** — "yes, pull the per-target numbers from the preprint"
The web fetch returned only medians; bioRxiv then blocked automated downloads (Cloudflare 1015).

**42 · 19:12** — "i add a pdf here as 2026.02.01.703098v1.full.pdf  is the pdf the paper"
Read the PDF.
* **Text:** extracted with pypdf.
* **Figures:** Fig. 2 and Fig. 3 rendered with PyMuPDF and read off visually.
* **6TFF:** exact values from the Fig. 3e table.
* **6WJR and 6Q57:** approximate values from the plots.
* **Inconsistency recorded:** the text says 6WJR has a representative within 5 Å; Fig. 2e shows none closer than
  about 15 Å.

### The page, the report and the repository

**25, 26 · 17:32, 17:34** — "create something like https://claude.ai/artifact/VhbkoJKuDuWgsjgpzCj6r4 for what we create with SimRNA and RibonanzaNet2" / "create the chart to show how the protocol od combination of the SimRNA and RibonanzaNet2 wors. likw what we create for the main protocol in the https://claude.ai/artifact/VhbkoJKuDuWgsjgpzCj6r4"
Built the protocol page: data flow, lanes, comparison arms and a SimRNA diagram.

**36 · 18:41** — "update https://claude.ai/artifact/BMLKaMXxh2xDf1qZCWeiBr too"
Added the "Run it on any sequence" section.

**Mid-task · 18:43** — "in the https://claude.ai/artifact/BMLKaMXxh2xDf1qZCWeiBr, only keep the result for the 6TFF, 6WJR, and 6Q57"
Removed 9HRF and 12CJ (two post-2025 aptamers that had also been run) from the page, and recomputed the means over the
three riboswitches.

**44 · 19:41** — "yes, add it to the page, and a comparsion of our result with RNAneal to the https://claude.ai/artifact/BMLKaMXxh2xDf1qZCWeiBr,  also remove all discussion about RNApro in it"
Added the pipeline results and section 7 (comparison with RNAnneal), and removed every mention of RNAPro.
* **Recorded change:** see the mid-task message under prompt 36: keeping only the three exercise targets changed the
  averages, and RibonanzaNet2 pairs then no longer gave a lower RMSD than EternaFold (16.4 vs 15.2 Å). The page says so.

**45 · 20:11** — "create report.pdf too which add the plot andresult overthere"
Wrote `src/simrna_report.py` (now `analysis/make_report.py`): four figures and a 6-page PDF.

**Mid-task · 20:13** — "why do we have The crystal-pairs arm is a ceiling control: it shows the best SimRNA can do with perfect pairs, and is the only place the crystal was ever fed to SimRNA. The protocol itself uses only the rnet2 arm in the "validation experiment, not part of the protocol""
Explained that the crystal-pairs arm tests whether the 3D errors come from the pairs or from SimRNA. Perfect pairs still
left the models about 12 Å off, so the limit is the 3D step.

**47 · 20:18** — "in part 6 "Every snapshot, not just the answer" check the snapshots < 15 Å (RNet2)"
On 6Q57, 40% of RibonanzaNet2 snapshots are within 15 Å of the crystal, against 5% with EternaFold pairs. Added to the
page and the PDF.

**48, 49 · 20:24, 20:26** — "what no pairs mean here?" / "what do you mean give the pair with crystal?"
Explained the two control arms: SimRNA with no restraints, and SimRNA given the crystal's base pairs (not its
coordinates) as a ceiling.

**50 · 20:30** — "in the part 7 "Compared with RNAnneal"  instead of using the langauge of our method use "RibonanzaNet2 + SimRNA" i want to send this protocol as a report to pro. tiwary"
Replaced "ours" with "RibonanzaNet2 + SimRNA" on the page and in the PDF.

**51 · 20:33** — "check teh whole protocol to see if we miss anything or any contradiction. forxample in the chat you forget to add the report.pdf"
Fixed 11 problems:
* inconsistent 6WJR values (21.9 vs 21.7 Å, from the random all-atom rebuild);
* 2.1 vs 2.2 Å for one value;
* two different INF definitions;
* the missing 4LVW caveat;
* run-time and accuracy wording;
* no link to the PDF;
* target-count differences from RNAnneal's Table I;
* the protocol tree implying the crystal is an input;
* a reference to a page the reader hasn't seen;
* the US-align version, found later.

**52 · 20:38** — "i want to put RibonanzaNet2 + SimRNA( with its code and result) on git hub."
Built this repository as a standalone package.
* **Checks:** the repackaged code predicts identical pairs; a local end-to-end run passed; the rebuilt PDF matches.
* **Size:** the trajectories (13 GB) are replaced by `snapshots.csv.gz`.
* **Not included:** SimRNA, the weights and the preprint PDF, for license and size reasons.
* **Correction:** one snapshot value at exactly 5 Å changed with 3-decimal rounding, so the full precision was stored.

**53 · 20:51** — "i want to add licence to my code"
MIT license, with an exception note for the Apache-2.0 RibonanzaNet file; SPDX headers in each file.

**54 · 20:55** — "i want to create README.html  for the github resporitoty that is like the protocole that we create https://claude.ai/artifact/BMLKaMXxh2xDf1qZCWeiBr so it is easier to understand the protocol"
Made a standalone copy of the page with this repository's paths.

**55–58, 63–66 · 21:04–21:26** — "my email for my git hub is [email] what els do you need to push it to a resporitoty", "my github username is Kmalekza", "wher sould i do the step 2", "i done it only put the visibility public", "i conform it", "set up the GitHub Pages for README.html", "i save the setting"
Pushed the repository and published the page.
* **Authorship:** commits re-authored to the GitHub email.
* **Key:** Claude Code's attempt to create an SSH key was blocked by its own safety check, so the user created the key
  and added it to GitHub.
* **Published:** pushed over SSH; GitHub Pages serves `docs/index.html`.

**67 · 21:30** — "i just see a plain text that we i clic on it doe not show the flow chat. is shows so cold line"
The Mermaid flowchart didn't render on GitHub Pages. It was replaced by a hand-drawn inline SVG with no external library.

**68–70 · 21:34–21:49** — the reply to Prof. Tiwary's exercise email; "where is TOOLS_AND_PROMPTS.md?"; "yes, create it and push to github"
Drafted the reply email and wrote this file.

## 3. Other prompts in the same session (not part of this protocol)

These are about the earlier RNAPro-based work, or are plain questions. They are listed for completeness:

* 1 "where are we and what was our last discussion"
* 2 "how RNApro works? …"
* 3 "what model weights means?"
* 9 "is it possible to use two The advisor models( RibonanzaNet2, and structRFM) or  more than one structRFM in RNAPro?"
* 10 "run the RhoFold+ template test on 9HRF and 12CJ"
* 11 "how its works for 6TFF, 6WJR, and 6Q57?"
* 12 "what did we do in the extentions?"
* 13 "in the RNAPro model, is it possibel instead of the chnaging the advisor(RibonanzaNet2) change the 3d model prediction?"
* 14 "yes, run Boltz-2 with RibonanzaNet2 pairs and crystal pairs"
* 15 "how about for 6TFF, 6WJR, and 6Q57?"
* 16 "run the tertiary contact test on 9HRF and 12CJ"
* mid-task, 00:33 "how RibonanzaNet2 create constrain? does it change any thin in the Boltz-2 protocaol and steps?"
* 31 "this work is great."

One lesson from these runs shaped this protocol. In Boltz-2 (prompts 14–16), base-pair and tertiary-contact
constraints, even from the crystal, barely changed the 3D result. That is why a builder that obeys restraints (SimRNA)
was chosen for the combination.
