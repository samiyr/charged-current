#!/bin/bash

CPUS=$1

sbatch --job-name sidis-scales-nnnpdf --cpus-per-task=$CPUS --array=0-16 generic.job sidis.differential.scales nNNPDF30_nlo_as_0118_A56_Z26
sbatch --job-name sidis-scales-nnnpdf --cpus-per-task=$CPUS --array=0-6 generic.job sidis.integrated.scales.min3 nNNPDF30_nlo_as_0118_A56_Z26
sbatch --job-name sidis-scales-nnnpdf --cpus-per-task=$CPUS --array=0-6 generic.job sidis.integrated.scales.min0 nNNPDF30_nlo_as_0118_A56_Z26