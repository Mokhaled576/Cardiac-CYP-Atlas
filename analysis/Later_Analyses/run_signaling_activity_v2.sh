#!/bin/bash
#SBATCH --account=def-aelkadi
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=18G
set -euo pipefail
cd "$HOME/scratch/HeartMap_CYP_Fresh"
module load StdEnv/2023 gcc/12.3 r-bundle-bioconductor/3.20
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
Rscript --vanilla fit_signaling_activity_v2.R "results/CYP_Disease_2962292/Cohort_${SLURM_ARRAY_TASK_ID}" "results/Signaling_Activity_${SLURM_ARRAY_JOB_ID}/Cohort_${SLURM_ARRAY_TASK_ID}"
