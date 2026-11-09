#!/bin/bash
#SBATCH --account=def-aelkadi
#SBATCH --time=00:25:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2
set -euo pipefail
cd "$HOME/scratch/HeartMap_CYP_Fresh"
module load StdEnv/2023 gcc/12.3 r-bundle-bioconductor/3.20
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
Rscript --vanilla donor_cyp_link_scores_v1.R "results/CYP_Disease_2962292/Cohort_${SLURM_ARRAY_TASK_ID}" "results/Donor_CYP_Links_${SLURM_ARRAY_JOB_ID}/Cohort_${SLURM_ARRAY_TASK_ID}"
