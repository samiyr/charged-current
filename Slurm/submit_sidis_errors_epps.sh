#!/bin/bash

CPUS=$1

sbatch --job-name=sidis-errors-epps --cpus-per-task=$CPUS sidis_errors.job EPPS21nlo_CT18Anlo_Fe56 0 10
sbatch --job-name=sidis-errors-epps --cpus-per-task=$CPUS sidis_errors.job EPPS21nlo_CT18Anlo_Fe56 11 20
sbatch --job-name=sidis-errors-epps --cpus-per-task=$CPUS sidis_errors.job EPPS21nlo_CT18Anlo_Fe56 21 30
sbatch --job-name=sidis-errors-epps --cpus-per-task=$CPUS sidis_errors.job EPPS21nlo_CT18Anlo_Fe56 31 40
sbatch --job-name=sidis-errors-epps --cpus-per-task=$CPUS sidis_errors.job EPPS21nlo_CT18Anlo_Fe56 41 50
sbatch --job-name=sidis-errors-epps --cpus-per-task=$CPUS sidis_errors.job EPPS21nlo_CT18Anlo_Fe56 51 60
sbatch --job-name=sidis-errors-epps --cpus-per-task=$CPUS sidis_errors.job EPPS21nlo_CT18Anlo_Fe56 61 70
sbatch --job-name=sidis-errors-epps --cpus-per-task=$CPUS sidis_errors.job EPPS21nlo_CT18Anlo_Fe56 71 80
sbatch --job-name=sidis-errors-epps --cpus-per-task=$CPUS sidis_errors.job EPPS21nlo_CT18Anlo_Fe56 81 90
sbatch --job-name=sidis-errors-epps --cpus-per-task=$CPUS sidis_errors.job EPPS21nlo_CT18Anlo_Fe56 91 100
sbatch --job-name=sidis-errors-epps --cpus-per-task=$CPUS sidis_errors.job EPPS21nlo_CT18Anlo_Fe56 101 106
