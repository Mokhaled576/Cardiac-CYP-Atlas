#!/bin/bash
#SBATCH --account=def-aelkadi
#SBATCH --time=00:20:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=2
set -euo pipefail
cd "$HOME/scratch/HeartMap_CYP_Fresh"
module load python/3.11.5
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
"$HOME/HeartMap/venv_states/bin/python" -u fit_donor_cyp_links_v2.py results/Donor_CYP_Links_2969410
