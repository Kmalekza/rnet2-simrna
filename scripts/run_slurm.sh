#!/bin/bash
# RNA sequence -> 3D structure with RibonanzaNet2 + SimRNA, on a Slurm cluster.
#
#   bash scripts/run_slurm.sh NAME SEQUENCE_OR_FASTA [REFERENCE.cif|pdb [CHAIN]]
#
# Submits three jobs that run one after another:
#   1 prep  RibonanzaNet2 base pairs + SimRNA inputs          (4 CPUs, ~1 min; uses a GPU if the job has one)
#   2 fold  RUNS SimRNA replica-exchange runs, 10 cores each   (~10-40 min each, in parallel)
#   3 post  cluster, rebuild all atoms, score, report           (1 CPU, ~1-5 min)
# Output in OUT/NAME/: model_1.pdb (the prediction), model_2..5.pdb, pairs.txt, report.md, report.json
#
# Settings (environment variables):
#   RUNS=8  ITERS=1600000  OUT=results/pipeline  PYTHON=python
#   SLURM_ACCOUNT=<account>  SLURM_PARTITION=<partition[,partition]>   (passed to sbatch when set)
#   RN2_DIR, SIMRNA_DIR, USALIGN, ANNOTATOR                            (see README)
#   PAIRS=<json>|none   validation only: use these base pairs instead of RibonanzaNet2
set -euo pipefail
[ $# -lt 2 ] && { sed -n 2,18p "$0"; exit 1; }
NAME=$1; SEQ=$2; REF=${3:-}; CHAIN=${4:-}
REPO=$(cd "$(dirname "$0")/.." && pwd)
PY=${PYTHON:-python}
export SIMRNA_DIR=$(realpath -m ${SIMRNA_DIR:-$REPO/third_party/SimRNA_64bitIntel_Linux})
export RN2_DIR=$(realpath -m ${RN2_DIR:-$REPO/weights/ribonanzanet2})
RUNS=${RUNS:-8}; ITERS=${ITERS:-1600000}
OUT=$(realpath -m ${OUT:-$REPO/results/pipeline})/$NAME
[ -f "$SEQ" ] && SEQ=$(realpath "$SEQ")
[ -n "$REF" ] && REF=$(realpath "$REF")
mkdir -p $OUT/logs
SB="sbatch --parsable ${SLURM_ACCOUNT:+-A $SLURM_ACCOUNT} ${SLURM_PARTITION:+-p $SLURM_PARTITION} --export=ALL"
EXTRA=""; [ -n "$REF" ] && EXTRA="--ref $REF"; [ -n "$CHAIN" ] && EXTRA="$EXTRA --chain $CHAIN"
[ -n "${PAIRS:-}" ] && EXTRA="$EXTRA --pairs $([ "$PAIRS" = none ] && echo none || realpath "$PAIRS")"
SSARG=$([ "${PAIRS:-}" = none ] && echo "" || echo "-S ss")

j1=$($SB -J ${NAME}_prep -c 4 --mem=16G -t 00:30:00 -o $OUT/logs/prep.out \
     --wrap "cd $REPO && $PY -m rnet2_simrna.cli prep $NAME '$SEQ' $OUT --iters $ITERS --runs $RUNS $EXTRA")

j2=$($SB -J ${NAME}_fold -c 10 --mem=8G -t 03:00:00 -a 1-$RUNS -d afterok:$j1 -o $OUT/logs/fold_%a.out \
     --wrap "cd $OUT/runs && ln -sfn $SIMRNA_DIR/data data && export OMP_NUM_THREADS=10 && \
             r=\$(printf r%02d \$SLURM_ARRAY_TASK_ID) && \
             $SIMRNA_DIR/SimRNA -c config.dat -E 10 -R \$((1000 + SLURM_ARRAY_TASK_ID)) -s seq $SSARG -o \$r > \$r.log 2>&1")

j3=$($SB -J ${NAME}_post -c 1 --mem=16G -t 01:00:00 -d afterok:$j2 -o $OUT/logs/post.out \
     --wrap "cd $REPO && $PY -m rnet2_simrna.cli post $OUT")

echo "submitted: prep $j1 -> fold $j2 ($RUNS runs) -> post $j3"
echo "results will be in $OUT  (report.md, model_1.pdb)"
