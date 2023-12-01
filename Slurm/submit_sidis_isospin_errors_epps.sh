#!/bin/bash

CPUS=$1

sbatch --job-name=sidis-isospin-errors-epps --cpus-per-task=$CPUS sidis_isospin_errors.job CT18ANLO 0 10
sbatch --job-name=sidis-isospin-errors-epps --cpus-per-task=$CPUS sidis_isospin_errors.job CT18ANLO 11 20
sbatch --job-name=sidis-isospin-errors-epps --cpus-per-task=$CPUS sidis_isospin_errors.job CT18ANLO 21 30
sbatch --job-name=sidis-isospin-errors-epps --cpus-per-task=$CPUS sidis_isospin_errors.job CT18ANLO 31 40
sbatch --job-name=sidis-isospin-errors-epps --cpus-per-task=$CPUS sidis_isospin_errors.job CT18ANLO 41 50
sbatch --job-name=sidis-isospin-errors-epps --cpus-per-task=$CPUS sidis_isospin_errors.job CT18ANLO 51 58