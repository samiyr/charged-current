#!/bin/bash

CPUS=$1

declare -a ARGS=(dis.differential.charm.errors dis.differential.total.errors dis.integrated.charm.errors dis.integrated.total.errors)

for ARG in "${ARGS[@]}"
do
	sbatch --job-name dis-errors-ncteq --cpus-per-task=$CPUS --array=0-38 generic.job $ARG nCTEQ15HQ_FullNuc_56_26
done