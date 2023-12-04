#!/bin/bash

CPUS=$1

declare -a ARGS=(sidis.differential.isospin.errors)

for ARG in "${ARGS[@]}"
do
	sbatch --job-name sidis-isospin-errors-ncteq --cpus-per-task=$CPUS --array=0-38 generic.job $ARG nCTEQ15HQ_FullNuc_1_1
done