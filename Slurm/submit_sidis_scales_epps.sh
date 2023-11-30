#!/bin/bash

CPUS=$1

sbatch --job-name=sidis-scales-epps --cpus-per-task=$CPUS sidis_scales.job EPPS21nlo_CT18Anlo_Fe56
