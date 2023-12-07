#!/bin/bash

CPUS=$1

sbatch --job-name sidis-scales-epps --cpus-per-task=$CPUS --array=0-16 generic.job sidis.differential.scales EPPS21nlo_CT18Anlo_Fe56
sbatch --job-name sidis-scales-epps --cpus-per-task=$CPUS --array=0-6 generic.job sidis.integrated.scales.min3 EPPS21nlo_CT18Anlo_Fe56
sbatch --job-name sidis-scales-epps --cpus-per-task=$CPUS --array=0-6 generic.job sidis.integrated.scales.min0 EPPS21nlo_CT18Anlo_Fe56