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

# Fit cutouts of one visit as a Slurm job at USDF:  sbatch examples/usdf_job.sh
#SBATCH --job-name=meas_crowd
# partition / account: see sacctmgr show assoc user=$USER
#SBATCH --partition=milano
#SBATCH --account=rubin:developers
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=meas_crowd_%j.log

source /sdf/group/rubin/sw/w_latest/loadLSST.bash
setup lsst_distrib
setup -k -r "$HOME/meas_crowd"          # where you unpacked and built the package
export OMP_NUM_THREADS=1

python "$HOME/meas_crowd/examples/run_visit_batch.py" \
    --visit 2025062100433 --detector 144 \
    --cutouts 1000,1000,300 2500,2000,300 \
    --n-jobs 2 --out "$HOME/meas_crowd_results/v2025062100433_d144"
