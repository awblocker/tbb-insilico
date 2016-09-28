#!/bin/bash

#SBATCH -J summarize_titration
#SBATCH --array=0-7
#SBATCH -n 12
#SBATCH -p airoldi
#SBATCH -t 0
#SBATCH --mem=94208
#SBATCH --exclusive
#SBATCH -o logs/summarize_%j.STDOUT
#SBATCH -e logs/summarize_%j.STDERR
#SBATCH --mail-user=ablocker@gmail.com
#SBATCH --mail-type=ALL

ID=$((SLURM_ARRAY_TASK_ID / 2 + 1))
EXPERIMENT=titration_control_${ID}

if test ! -f config/${EXPERIMENT}.yml
then
  echo config/${EXPERIMENT}.yml not found
  exit
fi

if ((SLURM_ARRAY_TASK_ID % 2 == 0))
then
  NULL=""
else
  NULL="--null"
fi

# Run base pair level summaries
cplate_summarise_mcmc --all --mmap ${NULL} config/${EXPERIMENT}.yml
# Run hyperparameter summaries
cplate_summarise_params_mcmc --all ${NULL} config/${EXPERIMENT}.yml
# Run cluster level summaries
cplate_summarise_clusters_mcmc --all ${NULL} config/${EXPERIMENT}.yml
