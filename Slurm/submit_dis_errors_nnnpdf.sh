#!/bin/bash

CPUS=$1

declare -a ARGS=(dis.differential.charm.errors dis.differential.total.errors dis.integrated.charm.errors dis.integrated.total.errors)

for ARG in "${ARGS[@]}"
do
	sbatch --job-name dis-errors-nnnpdf --cpus-per-task=$CPUS --array=0-200 generic.job $ARG nNNPDF30_nlo_as_0118_A56_Z26
done