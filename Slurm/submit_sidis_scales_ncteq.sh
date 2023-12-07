#!/bin/bash

CPUS=$1

sbatch --job-name sidis-scales-ncteq --cpus-per-task=$CPUS --array=0-16 generic.job sidis.differential.scales nCTEQ15HQ_FullNuc_56_26
sbatch --job-name sidis-scales-ncteq --cpus-per-task=$CPUS --array=0-6 generic.job sidis.integrated.scales.min3 nCTEQ15HQ_FullNuc_56_26
sbatch --job-name sidis-scales-ncteq --cpus-per-task=$CPUS --array=0-6 generic.job sidis.integrated.scales.min0 nCTEQ15HQ_FullNuc_56_26