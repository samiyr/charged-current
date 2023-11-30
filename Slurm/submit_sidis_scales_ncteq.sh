#!/bin/bash

CPUS=$1

sbatch --job-name=sidis-scales-ncteq --cpus-per-task=$CPUS sidis_scales.job nCTEQ15HQ_FullNuc_56_26
