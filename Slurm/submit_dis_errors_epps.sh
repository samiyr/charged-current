#!/bin/bash

CPUS=$1

declare -a ARGS=(dis.differential.charm.errors dis.differential.total.errors dis.integrated.charm.errors dis.integrated.total.errors)

for ARG in "${ARGS[@]}"
do
	sbatch --job-name dis-errors-epps --cpus-per-task=$CPUS --array=0-106 generic.job $ARG EPPS21nlo_CT18Anlo_Fe56
done