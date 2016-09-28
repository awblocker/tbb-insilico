#!/bin/bash

#SBATCH -J mcmc_titration
#SBATCH --array=1-4
#SBATCH -n 72
#SBATCH -p stats
#SBATCH -t 0
#SBATCH --mem=64345
#SBATCH --exclusive
#SBATCH -o logs/%j.STDOUT
#SBATCH -e logs/%j.STDERR
#SBATCH --mail-user=ablocker@gmail.com
#SBATCH --mail-type=ALL

EXPERIMENT=titration_control_${SLURM_ARRAY_TASK_ID}
if test ! -f config/${EXPERIMENT}.yml
then
  echo config/${EXPERIMENT}.yml not found
  exit
fi

for (( CHROM = 1; CHROM <= 16; CHROM++ ))
do
  RESULT_ID=titration-${SLURM_ARRAY_TASK_ID}_chrom$(printf '%02d' $CHROM)
  DRAWS=results/mcmc_draws_${RESULT_ID}.tar
  if test ! -f ${DRAWS}
  then
    mpirun -np 72 cplate_deconvolve_mcmc \
      -c ${CHROM} config/${EXPERIMENT}.yml \
      1>logs/mcmc_${RESULT_ID}.STDOUT \
      2>logs/mcmc_${RESULT_ID}.STDERR
    wait
  fi
  DRAWS=results/mcmc_draws_control_${RESULT_ID}.tar
  if test ! -f ${DRAWS}
  then
    mpirun -np 72 cplate_deconvolve_mcmc \
      --null \
      -c ${CHROM} config/${EXPERIMENT}.yml \
      1>logs/mcmc_${RESULT_ID}.STDOUT \
      2>logs/mcmc_${RESULT_ID}.STDERR
    wait
  fi
done
