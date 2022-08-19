#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 384
#SBATCH --mem=1024G
#SBATCH --job-name=3DPipe
#SBATCH --partition=smp


####################################################################################
#set up data directories
OUTPUT='/home/palasser/3d-dna'
SCRATCHDIR='/storage/scratch/palasser/3d-dna'

####################################################################################
# copy input files in scratchdir
# mkdir -m 750 /storage/scratch/palasser
# chmod -R g+w $SCRATCHDIR
# cp -R $OUTPUT /storage/scratch/palasser &
# id=$!

# cp /home/palasser/juicer-1.6/references/Triticum_aestivum_RENAN_v2.fasta $SCRATCHDIR/Triticum_aestivum_RENAN_v2.fasta &
# id2=$!

# cp /home/palasser/results/juicer-1.6_results/RENANV2/aligned/merged_nodups.txt $SCRATCHDIR/merged_nodups.txt &
# id3=$!
# wait $id $id2 $id3

####################################################################################
# 3d-dna cmd
cd $SCRATCHDIR
source env.sh

./run-asm-pipeline.sh -r 4 Triticum_aestivum_RENAN_v2.fasta merged_nodups_SUB_on_RENANV2_scaff_awk_sort.txt 2> 3ddna_on_RENAN_V2_pseudo_sub.log


####################################################################################
#Move from scratch to output
#mv -u $SCRATCHDIR $OUTPUT && rm -Rf $SCRATCHDIR
