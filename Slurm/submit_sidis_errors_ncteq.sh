#!/bin/bash

CPUS=$1

sbatch --job-name=sidis-errors-ncteq --cpus-per-task=$CPUS sidis_errors.job nCTEQ15HQ_FullNuc_56_26 0 5
sbatch --job-name=sidis-errors-ncteq --cpus-per-task=$CPUS sidis_errors.job nCTEQ15HQ_FullNuc_56_26 6 10
sbatch --job-name=sidis-errors-ncteq --cpus-per-task=$CPUS sidis_errors.job nCTEQ15HQ_FullNuc_56_26 11 15
sbatch --job-name=sidis-errors-ncteq --cpus-per-task=$CPUS sidis_errors.job nCTEQ15HQ_FullNuc_56_26 16 20
sbatch --job-name=sidis-errors-ncteq --cpus-per-task=$CPUS sidis_errors.job nCTEQ15HQ_FullNuc_56_26 21 25
sbatch --job-name=sidis-errors-ncteq --cpus-per-task=$CPUS sidis_errors.job nCTEQ15HQ_FullNuc_56_26 26 30
sbatch --job-name=sidis-errors-ncteq --cpus-per-task=$CPUS sidis_errors.job nCTEQ15HQ_FullNuc_56_26 31 38
