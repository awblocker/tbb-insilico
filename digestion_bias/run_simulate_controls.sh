#!/bin/bash

#SBATCH -J simulate_controls
#SBATCH --array=0-3
#SBATCH -p airoldi
#SBATCH -n 12
#SBATCH --exclusive
#SBATCH --hint=compute_bound
#SBATCH -t 129600
#SBATCH --mem=32768
#SBATCH --mail-user=ablocker@gmail.com
#SBATCH --mail-type=ALL

readonly REFERENCE_SEQUENCE="../Scerevisiae_chr.fna"
readonly READS=( $(ls ../XuData/titration/y_Sample*.txt) )
readonly REGIONS=( $(ls ../XuData/titration/regions_Sample*.txt) )
readonly LENGTHS=( $(ls ../XuData/titration/lengthDist_Sample_*.txt) )
readonly STATS=( $(ls ./joint_dinucleotide_statistics_Sample_NucDG_*.csv |
egrep -v unique) )
readonly OUTPUT_DIR=./control_$(basename ${READS[${SLURM_ARRAY_TASK_ID}]} |
sed 's/^y_//g' | sed 's/\.[a-zA-Z]*$//g')

CMD="python simulate_controls.py \
  --reference_sequence=$REFERENCE_SEQUENCE \
  --read_counts=${READS[${SLURM_ARRAY_TASK_ID}]} \
  --regions=${REGIONS[${SLURM_ARRAY_TASK_ID}]} \
  --lengths=${LENGTHS[${SLURM_ARRAY_TASK_ID}]} \
  --dinucleotide_statistics=${STATS[${SLURM_ARRAY_TASK_ID}]} \
  --output_dir=$OUTPUT_DIR"
echo $CMD
eval $CMD
