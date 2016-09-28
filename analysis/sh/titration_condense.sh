#!/bin/bash

#SBATCH -J condense_titration
#SBATCH --array=1-4
#SBATCH -n 1
#SBATCH -p airoldi
#SBATCH -t 0
#SBATCH --mem=32768
#SBATCH --exclusive
#SBATCH -o logs/condense_%j.STDOUT
#SBATCH -e logs/condense_%j.STDERR
#SBATCH --mail-user=ablocker@gmail.com
#SBATCH --mail-type=ALL

EXPERIMENT=titration_control_${SLURM_ARRAY_TASK_ID}

Rscript R/extract_detections.R --outdir=results/ \
  config/${EXPERIMENT}.yml
Rscript R/extract_clusters.R --outdir=results/ \
  config/${EXPERIMENT}.yml

