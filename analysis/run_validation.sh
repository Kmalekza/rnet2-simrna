#!/bin/bash
# Validation runs: the same SimRNA protocol with four sources of base pairs, on 6TFF, 6WJR and 6Q57 (12 runs of 8 x 10 replicas).
#   nocon   no pairs (SimRNA alone)
#   eterna  EternaFold MEA pairs          (data/validation_pairs_<target>_eterna.json)
#   rnet2   RibonanzaNet2 L25H1 pairs     (predicted, as in the pipeline)
#   xtal    crystal pairs, ceiling control (data/validation_pairs_<target>_xtal.json)
# Only the pairs differ between arms. The crystal structure is also passed in, for scoring only.
#   bash analysis/run_validation.sh            (Slurm; then: python analysis/collect_validation.py)
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
export OUT=${OUT:-$REPO/results/validation_runs}
tail -n +2 $REPO/data/targets.csv | while IFS=, read -r t structure chain length seq desc; do
  for a in nocon eterna rnet2 xtal; do
    case $a in
      nocon) P=none ;;
      rnet2) P= ;;
      *)     P=$REPO/data/validation_pairs_${t}_${a}.json ;;
    esac
    PAIRS=$P bash $REPO/scripts/run_slurm.sh ${t}_${a} $seq $REPO/$structure $chain
  done
done
