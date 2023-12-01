#!/bin/bash

CPUS=$1

sbatch --job-name=sidis-isospin-errors-ncteq --cpus-per-task=$CPUS sidis_isospin_errors.job nCTEQ15HQ_FullNuc_1_1 0 5
sbatch --job-name=sidis-isospin-errors-ncteq --cpus-per-task=$CPUS sidis_isospin_errors.job nCTEQ15HQ_FullNuc_1_1 6 10
sbatch --job-name=sidis-isospin-errors-ncteq --cpus-per-task=$CPUS sidis_isospin_errors.job nCTEQ15HQ_FullNuc_1_1 11 15
sbatch --job-name=sidis-isospin-errors-ncteq --cpus-per-task=$CPUS sidis_isospin_errors.job nCTEQ15HQ_FullNuc_1_1 16 20
sbatch --job-name=sidis-isospin-errors-ncteq --cpus-per-task=$CPUS sidis_isospin_errors.job nCTEQ15HQ_FullNuc_1_1 21 25
sbatch --job-name=sidis-isospin-errors-ncteq --cpus-per-task=$CPUS sidis_isospin_errors.job nCTEQ15HQ_FullNuc_1_1 26 30
sbatch --job-name=sidis-isospin-errors-ncteq --cpus-per-task=$CPUS sidis_isospin_errors.job nCTEQ15HQ_FullNuc_1_1 31 38