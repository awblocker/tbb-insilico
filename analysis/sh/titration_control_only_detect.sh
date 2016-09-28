#!/bin/bash

#SBATCH -J detect_titration_control_only
#SBATCH --array=1-4
#SBATCH -n 1
#SBATCH -p airoldi
#SBATCH -t 0
#SBATCH --mem=32768
#SBATCH -o logs/summarize_%j.STDOUT
#SBATCH -e logs/summarize_%j.STDERR
#SBATCH --mail-user=ablocker@gmail.com
#SBATCH --mail-type=ALL

EXPERIMENT=titration_control_only_${SLURM_ARRAY_TASK_ID}

if test ! -f config/${EXPERIMENT}.yml
then
  echo config/${EXPERIMENT}.yml not found
  exit
fi

# Run detection for control samples.
cplate_detect_mcmc --all config/${EXPERIMENT}.yml

# Condense results.
Rscript R/extract_detections.R --outdir=results/ \
  config/${EXPERIMENT}.yml
Rscript R/extract_clusters.R --outdir=results/ \
  config/${EXPERIMENT}.yml
