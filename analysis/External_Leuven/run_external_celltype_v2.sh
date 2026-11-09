#!/bin/bash
#SBATCH --account=def-aelkadi
#SBATCH --time=00:30:00
#SBATCH --mem=56G
#SBATCH --cpus-per-task=2
set -euo pipefail
cd "$HOME/scratch/HeartMap_CYP_Fresh/external_validation/GSE298023"
module load StdEnv/2023 gcc/12.3 r-bundle-bioconductor/3.20
if [ -d "$HOME/R_libs_4.4" ]; then export R_LIBS_USER="$HOME/R_libs_4.4"; fi
test -s PUBLISHED_OBJECT_INSPECTED.txt
Rscript --vanilla fit_external_celltype_v2.R GSE298023_HH_dataset_seurat_object_2026.rds sample1 seurat_clusters "Validation_$SLURM_JOB_ID"
