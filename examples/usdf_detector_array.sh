#!/bin/bash
# Whole-detector run at USDF as a Slurm array: N tasks share the tiles, then one stitch job.
#
#   VISIT=2025062100433 DETECTOR=144 bash examples/usdf_detector_array.sh
#
# Check partition / account for your allocation (`sacctmgr show assoc user=$USER`).
set -euo pipefail
VISIT=${VISIT:-2025062100433}
DETECTOR=${DETECTOR:-144}
NTASKS=${NTASKS:-25}                  # ~100 tiles of 400 px per LSSTCam detector -> ~4 per task
OUT=${OUT:-$HOME/meas_crowd_results/v${VISIT}_d${DETECTOR}}
PKG=${PKG:-$HOME/meas_crowd}
mkdir -p "$OUT/logs"

SETUP="source /sdf/group/rubin/sw/w_latest/loadLSST.bash && setup lsst_distrib && setup -k -r $PKG && export OMP_NUM_THREADS=1"

FIT=$(sbatch --parsable --job-name=crowd_fit --partition=milano --account=rubin:developers \
  --array=0-$((NTASKS - 1)) --cpus-per-task=1 --mem=8G --time=04:00:00 \
  --output="$OUT/logs/fit_%a.log" \
  --wrap="$SETUP && python $PKG/bin.src/crowdDetector.py fit $VISIT $DETECTOR --out $OUT")

sbatch --dependency=afterany:$FIT --job-name=crowd_stitch --partition=milano \
  --account=rubin:developers --cpus-per-task=1 --mem=16G --time=00:30:00 \
  --output="$OUT/logs/stitch.log" \
  --wrap="$SETUP && python $PKG/bin.src/crowdDetector.py stitch --out $OUT && python $PKG/bin.src/crowdDetector.py plot --out $OUT"

echo "submitted fit array $FIT ($NTASKS tasks) + stitch; results in $OUT"
