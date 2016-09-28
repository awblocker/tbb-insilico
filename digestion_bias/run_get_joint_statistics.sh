#!/bin/bash

#SBATCH -J joint_statistics
#SBATCH --array=0-3
#SBATCH -n 1
#SBATCH -t 0
#SBATCH --mem=4096
#SBATCH --mail-user=ablocker@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -p airoldi

WINDOW_SIZE=10
REFERENCE_GENOME="../Scerevisiae_chr.fna"
ALIGNED_READS=( $(ls ../XuData/titration/Sample*.sorted.bam) )

READS=${ALIGNED_READS[${SLURM_ARRAY_TASK_ID}]}

DESTINATION=$(basename $READS)
DESTINATION=joint_dinucleotide_statistics_${DESTINATION%%.*}.csv
python get_joint_dinucleotide_statistics.py \
    --aligned_reads_source=$READS \
    --reference_sequence_source=$REFERENCE_GENOME \
    --destination=$DESTINATION \
    --window_size=$WINDOW_SIZE
