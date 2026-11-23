#!/bin/bash
# RNA sequence -> 3D structure with RibonanzaNet2 + SimRNA, on one machine (no Slurm).
#
#   bash scripts/run_local.sh NAME SEQUENCE_OR_FASTA [REFERENCE.cif|pdb [CHAIN]]
#
# Runs RUNS SimRNA runs side by side with THREADS OpenMP threads each (default 8 runs x 10 threads = 80 cores;
# lower RUNS or THREADS on a smaller machine, or set ITERS lower for a quick test, e.g. ITERS=160000).
# Same settings and outputs as run_slurm.sh.
set -euo pipefail
[ $# -lt 2 ] && { sed -n 2,9p "$0"; exit 1; }
NAME=$1; SEQ=$2; REF=${3:-}; CHAIN=${4:-}
REPO=$(cd "$(dirname "$0")/.." && pwd)
PY=${PYTHON:-python}
export SIMRNA_DIR=$(realpath -m ${SIMRNA_DIR:-$REPO/third_party/SimRNA_64bitIntel_Linux})
export RN2_DIR=$(realpath -m ${RN2_DIR:-$REPO/weights/ribonanzanet2})
RUNS=${RUNS:-8}; ITERS=${ITERS:-1600000}; THREADS=${THREADS:-10}
OUT=$(realpath -m ${OUT:-$REPO/results/pipeline})/$NAME
[ -f "$SEQ" ] && SEQ=$(realpath "$SEQ")
EXTRA=""; [ -n "$REF" ] && EXTRA="--ref $(realpath $REF)"; [ -n "$CHAIN" ] && EXTRA="$EXTRA --chain $CHAIN"
[ -n "${PAIRS:-}" ] && EXTRA="$EXTRA --pairs $([ "$PAIRS" = none ] && echo none || realpath "$PAIRS")"
SSARG=$([ "${PAIRS:-}" = none ] && echo "" || echo "-S ss")

(cd $REPO && $PY -m rnet2_simrna.cli prep $NAME "$SEQ" $OUT --iters $ITERS --runs $RUNS $EXTRA)
cd $OUT/runs && ln -sfn $SIMRNA_DIR/data data
for i in $(seq 1 $RUNS); do
  r=$(printf r%02d $i)
  OMP_NUM_THREADS=$THREADS $SIMRNA_DIR/SimRNA -c config.dat -E 10 -R $((1000 + i)) -s seq $SSARG -o $r > $r.log 2>&1 &
done
wait
cd $REPO && $PY -m rnet2_simrna.cli post $OUT
