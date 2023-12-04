#!/bin/bash

CPUS=$1

declare -a ARGS=(sidis.differential.isospin.errors)

for ARG in "${ARGS[@]}"
do
	sbatch --job-name sidis-isospin-errors-epps --cpus-per-task=$CPUS --array=0-58 generic.job $ARG CT18ANLO
done