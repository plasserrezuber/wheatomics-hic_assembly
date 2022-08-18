#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --job-name=juicer_scaff_pseudo
#SBATCH --partition=smp
#SBATCH --dependency=afterok:36191520

####################################################################################
#Set up the temporary directory
SCRATCHDIR='/storage/scratch/palasser/juicer-1.6_scaff_sub'
#mkdir -p -m 770 $SCRATCHDIR

cd $SCRATCHDIR

####################################################################################
# #copy input files in scratchdir

# bwa index $SCRATCHDIR/references/RENANV2_SCAFFOLDS_PSEUDO_v1_order.fasta
# cut -f1,2 $SCRATCHDIR/references/RENANV2_SCAFFOLDS_PSEUDO_v1_order.fasta.fai > $SCRATCHDIR/references/RENANV2_SCAFFOLDS_PSEUDO_v1_order.fasta.sizes

####################################################################################
#generate restriction site file
#./scripts/generate_site_positions.py Arima RENANV2_SCAFFOLDS_PSEUDO_v1_order ./references/RENANV2_SCAFFOLDS_PSEUDO_v1_order.fasta
#wait $id

####################################################################################
# juicer cmd
./scripts/juicer.sh \
-g "RENANV2_SCAFFOLDS_PSEUDO" \
-D $SCRATCHDIR \
-d $SCRATCHDIR/RENANV2_scaff_sub \
-s "Arima" \
-S "chimeric" \
-p $SCRATCHDIR/references/RENANV2_SCAFFOLDS_PSEUDO_v1_order.fasta.sizes \
-y $SCRATCHDIR/references/RENANV2_SCAFFOLDS_PSEUDO_v1_order_Arima.txt \
-z $SCRATCHDIR/references/RENANV2_SCAFFOLDS_PSEUDO_v1_order.fasta

####################################################################################
#Move from scratch to output
#mv -u $SCRATCHDIR/RENANV2 $OUTPUT/ && rm -r $SCRATCHDIR
