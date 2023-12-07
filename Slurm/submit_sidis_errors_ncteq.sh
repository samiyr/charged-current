#!/bin/bash

CPUS=$1

declare -a ARGS=(sidis.differential.errors sidis.integrated.errors.min3 sidis.integrated.errors.min0)

for ARG in "${ARGS[@]}"
do
	sbatch --job-name sidis-errors-ncteq --cpus-per-task=$CPUS --array=0-38 generic.job $ARG nCTEQ15HQ_FullNuc_56_26
done

sbatch --job-name sidis-isospin-errors-ncteq --cpus-per-task=$CPUS --array=0-38 generic.job sidis.differential.isospin.errors nCTEQ15HQ_FullNuc_1_1