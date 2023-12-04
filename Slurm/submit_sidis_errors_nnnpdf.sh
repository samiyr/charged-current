#!/bin/bash

CPUS=$1

declare -a ARGS=(sidis.differential.errors sidis.integrated.errors.min3 sidis.integrated.errors.min0)

for ARG in "${ARGS[@]}"
do
	sbatch --job-name sidis-errors-nnnpdf --cpus-per-task=$CPUS --array=0-200 generic.job $ARG nNNPDF30_nlo_as_0118_A56_Z26
done