#!/bin/bash

CPUS=$1

declare -a ARGS=(dis.differential.charm.scales dis.differential.total.scales dis.integrated.charm.scales dis.integrated.total.scales)

for ARG in "${ARGS[@]}"
do
	sbatch --job-name dis-scales-epps --cpus-per-task=$CPUS --array=0-6 generic.job $ARG nCTEQ15HQ_FullNuc_56_26
done