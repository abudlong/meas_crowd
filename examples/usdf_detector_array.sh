#!/bin/bash
# This file is part of meas_crowd.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Whole detector as a Slurm array at USDF:  VISIT=... DETECTOR=... bash examples/usdf_detector_array.sh
set -euo pipefail
VISIT=${VISIT:-2025062100433}
DETECTOR=${DETECTOR:-144}
NTASKS=${NTASKS:-25}                  # ~100 tiles of 400 px per LSSTCam detector -> ~4 per task
OUT=${OUT:-$HOME/meas_crowd_results/v${VISIT}_d${DETECTOR}}
PKG=${PKG:-$HOME/meas_crowd}
mkdir -p "$OUT/logs"

SETUP="source /sdf/group/rubin/sw/w_latest/loadLSST.bash && setup lsst_distrib && setup -k -r $PKG && export OMP_NUM_THREADS=1"

# partition / account: see sacctmgr show assoc user=$USER
FIT=$(sbatch --parsable --job-name=crowd_fit --partition=milano --account=rubin:developers \
  --array=0-$((NTASKS - 1)) --cpus-per-task=1 --mem=8G --time=04:00:00 \
  --output="$OUT/logs/fit_%a.log" \
  --wrap="$SETUP && python $PKG/bin.src/crowdDetector.py fit $VISIT $DETECTOR --out $OUT")

sbatch --dependency=afterany:$FIT --job-name=crowd_stitch --partition=milano \
  --account=rubin:developers --cpus-per-task=1 --mem=16G --time=00:30:00 \
  --output="$OUT/logs/stitch.log" \
  --wrap="$SETUP && python $PKG/bin.src/crowdDetector.py stitch --out $OUT && python $PKG/bin.src/crowdDetector.py plot --out $OUT"

echo "submitted fit array $FIT ($NTASKS tasks) + stitch; results in $OUT"
