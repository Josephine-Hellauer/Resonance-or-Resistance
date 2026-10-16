#!/bin/bash
#SBATCH -J phase_ii_s_sweep
#SBATCH -N 1
#SBATCH --tasks-per-node=1
#SBATCH -t 20:00:00
#SBATCH -p small-short

cd "$SLURM_SUBMIT_DIR"
source ~/.bashrc
conda activate soc-env

python -u src/phase_ii_s_dual_contagion_sweep.py
