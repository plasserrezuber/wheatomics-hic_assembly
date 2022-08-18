#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 10
#SBATCH --mem=40G #pour 10 job en parallel x 4G par job
#SBATCH --job-name=flagstt
#SBATCH --partition=smp

####################################################################################
#set up data directories
OUTPUT='/storage/scratch/palasser/juicer-1.6/RENANV2_pseudo/splits'

####################################################################################
module load gcc/8.1.0 samtools/1.9 parallel/20151222

####################################################################################
parallel -j 10 --bar "samtools flagstat {} > {}.flagstat.out 2> {}.flagstat.err" ::: $(find $OUTPUT -name "*.fastq.gz.bam")
