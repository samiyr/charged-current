#!/bin/bash

CPUS=$1

sbatch --job-name=sidis-others-epps --cpus-per-task=$CPUS sidis_others.job EPPS21nlo_CT18Anlo_Fe56
