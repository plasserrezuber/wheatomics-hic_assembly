#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --job-name=juicerV2
#SBATCH --partition=smp
## SBATCH --dependency=afterok:32274149

####################################################################################
## Set up the temporary directory
SCRATCHDIR='/storage/scratch/palasser/juicer-1.6'
INPUT='/home/palasser/results/juicer-1.6_results/references'
mkdir -p -m 770 $SCRATCHDIR

cd $SCRATCHDIR
# mkdir -p $SCRATCHDIR/Renan_v2_pseudo_v2/fastq
# mkdir $SCRATCHDIR/references

####################################################################################
## before (time: 5h)
#sbatch -p smp -c 1 --mem=32G --wrap="ml bwa; bwa index /storage/scratch/palasser/juicer-1.6/references/Renan_v13_v2.pseudo.v2.fa"

####################################################################################
## copy input files in scratchdir

# cp -R /home/palasser/juicer-1.6/scripts $SCRATCHDIR
# for f in "00" "01" "02"; do cp /home/palasser/results/juicer-1.6_results/RENANV2_fastq/fastq/${f}*.fastq.gz $SCRATCHDIR/Renan_v2_pseudo_v2/fastq/; done
# cp $INPUT/Renan_v13_v2.pseudo.v2.fa $SCRATCHDIR/references/
# cp $INPUT/Renan_v13_v2.pseudo.v2.fa.fai $SCRATCHDIR/references/
# cut -f1,2 $SCRATCHDIR/references/Renan_v13_v2.pseudo.v2.fa.fai > $SCRATCHDIR/references/Renan_v13_v2.pseudo.v2.fa.sizes

## generate restriction site file (time: 10h)
# $SCRATCHDIR/scripts/generate_site_positions.py Arima Renan_v13_v2.pseudo.v2 $SCRATCHDIR/references/Renan_v13_v2.pseudo.v2.fa
# mv $SCRATCHDIR/Renan_v13_v2.pseudo.v2_Arima.txt $SCRATCHDIR/references/

## copy BWA index
# cp $INPUT/Renan_v13_v2.pseudo.v2.fa.amb $SCRATCHDIR/references/
# cp $INPUT/Renan_v13_v2.pseudo.v2.fa.ann $SCRATCHDIR/references/
# cp $INPUT/Renan_v13_v2.pseudo.v2.fa.bwt $SCRATCHDIR/references/
# cp $INPUT/Renan_v13_v2.pseudo.v2.fa.pac $SCRATCHDIR/references/
# cp $INPUT/Renan_v13_v2.pseudo.v2.fa.sa $SCRATCHDIR/references/


####################################################################################
# juicer cmd
./scripts/juicer.sh \
-g "Renan_v13_v2.pseudo.v2" \
-D $SCRATCHDIR \
-d $SCRATCHDIR/Renan_v2_pseudo_v2 \
-s "Arima" \
-S "merge" \
-p $SCRATCHDIR/references/Renan_v13_v2.pseudo.v2.fa.sizes \
-y $SCRATCHDIR/references/Renan_v13_v2.pseudo.v2_Arima.txt \
-z $SCRATCHDIR/references/Renan_v13_v2.pseudo.v2.fa 

####################################################################################
#Move from scratch to output

