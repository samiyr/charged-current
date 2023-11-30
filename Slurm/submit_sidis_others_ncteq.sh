#!/bin/bash

CPUS=$1

sbatch --job-name=sidis-others-ncteq --cpus-per-task=$CPUS sidis_others.job nCTEQ15HQ_FullNuc_56_26
