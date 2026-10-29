#!/bin/bash
# Slurm batch job for USDF.  Submit from a rubin-devl node:   sbatch examples/usdf_job.sh
# Check the partition / account names that apply to you (`sacctmgr show assoc user=$USER`).
#SBATCH --job-name=meas_crowd
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
