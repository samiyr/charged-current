#!/bin/bash

CPUS=$1

declare -a ARGS=(sidis.differential.errors sidis.integrated.errors.min3 sidis.integrated.errors.min0)

for ARG in "${ARGS[@]}"
do
	sbatch --job-name sidis-errors-epps --cpus-per-task=$CPUS --array=0-106 generic.job $ARG EPPS21nlo_CT18Anlo_Fe56
done