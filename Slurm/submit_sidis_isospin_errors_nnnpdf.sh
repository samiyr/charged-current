#!/bin/bash

CPUS=$1

declare -a ARGS=(sidis.differential.isospin.errors)

for ARG in "${ARGS[@]}"
do
	sbatch --job-name sidis-isospin-errors-nnnpdf --cpus-per-task=$CPUS --array=0-200 generic.job $ARG nNNPDF30_nlo_as_0118_p
done