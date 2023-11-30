#!/bin/bash

CPUS=$1

sbatch --job-name=sidis-scales-nnnpdf --cpus-per-task=$CPUS sidis_scales.job nNNPDF30_nlo_as_0118_A56_Z26
