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

# Compare two settings on tiles across one detector at USDF:  OUT=... bash examples/usdf_validate.sh
set -euo pipefail
OUT=${OUT:-/sdf/home/a/abudlong/project/repo-main-logs/crowdedFields/october2026/crowded_fields/7a-newTest}
PKG=${PKG:-$HOME/meas_crowd}
mkdir -p "$OUT/logs"
NTASKS=$(python3 "$PKG/examples/validate_settings.py" list | wc -l)   # settings x tiles (12)

SETUP="source /sdf/group/rubin/sw/w_latest/loadLSST.bash && setup lsst_distrib && setup -k -r $PKG && export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1"

# partition / account: see sacctmgr show assoc user=$USER
FIT=$(sbatch --parsable --job-name=crowd_validate --partition=milano --account=rubin:developers \
  --array=0-$((NTASKS - 1)) --cpus-per-task=1 --mem=8G --time=02:00:00 \
  --output="$OUT/logs/task_%a.log" \
  --wrap="$SETUP && python $PKG/examples/validate_settings.py fit --out $OUT")

sbatch --dependency=afterany:$FIT --job-name=crowd_compare --partition=milano \
  --account=rubin:developers --cpus-per-task=1 --mem=8G --time=00:30:00 \
  --output="$OUT/logs/compare.log" \
  --wrap="$SETUP && python $PKG/examples/validate_settings.py compare --out $OUT"

echo "submitted $NTASKS fit tasks ($FIT) + compare; results in $OUT/validation_summary.txt"
