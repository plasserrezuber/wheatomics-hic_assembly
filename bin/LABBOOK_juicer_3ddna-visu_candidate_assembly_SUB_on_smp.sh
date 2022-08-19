#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --job-name=ju_SUBpseudo
#SBATCH --partition=smp
# #SBATCH --dependency=afterok:32274149

####################################################################################
#Set up the temporary directory
SCRATCHDIR='/storage/scratch/palasser/juicer-1.6_sub'
#mkdir -p -m 770 $SCRATCHDIR

cd $SCRATCHDIR
#mkdir references
#mkdir -p RENANV2_pseudo_sub/fastq

####################################################################################
# #copy input files in scratchdir
#cp -R /home/palasser/juicer-1.6/scripts $SCRATCHDIR
#cp /home/celmonat/juicer-1.6/generate_site_positions.py $SCRATCHDIR/scripts

#for f in "00" "01" "02"; do cp /home/palasser/results/juicer-1.6_results/RENANV2_fastq/splits/${f}*.fastq.gz $SCRATCHDIR/RENANV2_pseudo_sub/fastq/ ; done
#cp /home/palasser/results/juicer-1.6_results/references/RENANV2_SCAFFOLDS_PSEUDO_v1_order.fasta* $SCRATCHDIR/references/


#cut -f1,2 $SCRATCHDIR/references/RENANV2_SCAFFOLDS_PSEUDO_v1_order.fasta.fai > $SCRATCHDIR/references/RENANV2_SCAFFOLDS_PSEUDO_v1_order.fasta.sizes

####################################################################################
#generate restriction site file
#./scripts/generate_site_positions.py Arima RENANV2_SCAFFOLDS_PSEUDO_v1_order $SCRATCHDIR/references/RENANV2_SCAFFOLDS_PSEUDO_v1_order.fasta
#id=$!

#wait $id

####################################################################################
# juicer cmd
./scripts/juicer.sh \
-g "RENANV2_pseudo" \
-D $SCRATCHDIR \
-d $SCRATCHDIR/RENANV2_pseudo_sub \
-s "Arima" \
-S "dedup" \
-p $SCRATCHDIR/references/Renan_v13_v2.pseudo.v1.fa.sizes \
-y $SCRATCHDIR/references/RENANV2_pseudo_Arima.txt \
-z $SCRATCHDIR/references/Renan_v13_v2.pseudo.v1.fa #\
#-t 16


####################################################################################
#Move from scratch to output
#mv -u $SCRATCHDIR/RENANV2 $OUTPUT/ && rm -r $SCRATCHDIR


####################################################################################
## /!\/!\ ON NE PEUT PAS IMPORTER UN FICHIER ASSEMBLY FAIT A LA MAIN OU PAR 3d-dna/utils/generate-assembly-file-from-fasta.awk
## POUR ANNOTER ET CORRIGER L'ASSEMBLAGE D'UN FICHIER .hic OBTENU AVEC JUICER !!! soit 40G de fastq.gz au lieu de 285G
## -> Quand on a un assemblage candidat avec des pseudomol, il faut utiliser des scripts 3d-dna pour creer conjointement .assembly et .hic, cf ci-dessous.
## cf Cookbook_Juicer_3D-DNA.pdf page 5, et ce post: https://groups.google.com/g/3d-genomics/c/6OCqwTr9_ho

## TRAVAIL SUR LES PSEUDO
## LA SOLUTION est de travailler sur un sous set de donnees et cela donne le meme resultat que sur la totalite des donnees !!  
## alignement bwa mem (-t 16 -SP5M) d'un SUBSET de PE-reads HiC (24 fichiers fastq = 3 split * R1/R2 * 4 libraries) sur les pseudo de RENANV2 avec le pipe Juicer
## obtention d'un fichier merged_nodups.txt avec juicer et mis en forme awk et sort (separateur " " et non "\t")
awk -F' ' '{if ($2>$6 && $9>0 && $12>0) { print $5,$6,$7,$8,$1,$2,$3,$4,$9,$10,$11,$12,$13,$14,$15,$16 }; if ($2<=$6 && $9>0 && $12>0) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16 } }' /storage/scratch/palasser/juicer-1.6_sub/RENANV2_pseudo_sub/aligned/merged_nodups.txt > /storage/scratch/palasser/juicer-1.6_sub/RENANV2_pseudo_sub/aligned/merged_nodups_mapQ1_awk.txt
sort -t' ' -k2,2d -k6,6d /storage/scratch/palasser/juicer-1.6_sub/RENANV2_pseudo_sub/aligned/merged_nodups_mapQ1_awk.txt > /storage/scratch/palasser/juicer-1.6_sub/RENANV2_pseudo_sub/aligned/merged_nodups_RENANV2_pseudo_SUB_mapQ1_awk_sort.txt

## obtention d'un fichier .assembly avec /3d-dna/utils/generate-assembly-file-from-fasta.awk du pipe 3d-dna
awk -f /storage/scratch/palasser/3d-dna/utils/generate-assembly-file-from-fasta.awk /storage/scratch/palasser/juicer-1.6_sub/references/Renan_v13_v2.pseudo.v1.fa > /storage/scratch/palasser/juicer-1.6_sub/RENANV2_pseudo_sub/aligned/RENANV2_pseudov1_SUB_run-assembly-visualizer.assembly

## obtention d'un fichier .hic avec /3d-dna/visualize/run-assembly-visualizer.sh du pipe 3d-dna
sbatch -p smp -c 256 --mem=500G --wrap="ml java/oracle-1.8.0_45; export IBM_JAVA_OPTIONS='-Xmx500g -Xgcthreads256'; export _JAVA_OPTIONS='-Xmx500g -Xms500g'; /storage/scratch/palasser/3d-dna/visualize/run-assembly-visualizer.sh /storage/scratch/palasser/juicer-1.6_sub/RENANV2_pseudo_sub/aligned/RENANV2_pseudov1_SUB_run-assembly-visualizer.assembly /storage/scratch/palasser/juicer-1.6_sub/RENANV2_pseudo_sub/aligned/merged_nodups_RENANV2_pseudo_SUB_mapQ1_awk_sort.txt"

## ==> resultat dans juicebox JBAT: affichage des features des chromosomes, deja une base pour apporter des corrections et les sauvegarder
## mais pas d'affichage des scaffolds


## TRAVAIL SUR LES SCAFFOLDS POUR POUVOIR RETRAVAILLER CETTE ECHELLE D'ASSEMBLAGE DANS JUICEBOX ASSEMBLY TOOL (=JBAT)

## FASTA DES SCAFFOLDS DANS ORDRE PSEUDO
## mise en forme d'un fichier fasta des scaffolds dans l'ordre d'assemblage des pseudomolécules de RENANV2
## bed
for chr in "1A" "1B" "1D" "2A" "2B" "2D" "3A" "3B" "3D" "4A" "4B" "4D" "5A" "5B" "5D" "6A" "6B" "6D" "7A" "7B" "7D";
do
    gawk -v chrom=$chr '{ if ($3~chrom) { sum+=$2; print $3, sum, sum, $1, $4 } }' scaffAssignment.50000.csv.order.orientation.csv >> RENANV2_PSEUDO_scaff_v1.bed
done
## puis mise en forme du bed dans excel (premiere cellule deuxieme colonne = "0" a chaque debut de chrom, decalage vers le bas, suppression derniere cellule)

## fasta des scaffolds
ml bedtools
bedtools getfasta -name -fi /home/celmonat/Renan_v13_v2.pseudo.v1.fa \
-bed /home/palasser/data/RENAN_v2_pseudo/RENANV2_PSEUDO_v1_scaff.bed \
-fo /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_PSEUDO_v1_order.fasta

cp /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_PSEUDO_v1_order.fasta /storage/scratch/palasser/juicer-1.6_scaff_sub/references/
sbatch -p smp -c 1 --mem=500G --wrap="/storage/scratch/palasser/juicer-1.6_scaff_sub/scripts/generate_site_positions.py Arima RENANV2_SCAFFOLDS_PSEUDO_v1_order /storage/scratch/palasser/juicer-1.6_scaff_sub/references/RENANV2_SCAFFOLDS_PSEUDO_v1_order.fasta"
sbatch -p smp -c 1 --mem=256G --wrap="ml bwa; bwa index /storage/scratch/palasser/juicer-1.6_scaff_sub/references/RENANV2_SCAFFOLDS_PSEUDO_v1_order.fasta"
sbatch -p smp -c 1 --mem=128G --dependency=afterok:36191520 --wrap="cp /storage/scratch/palasser/juicer-1.6_scaff_sub/references/RENANV2_SCAFFOLDS_PSEUDO_v1_order* /home/palasser/results/juicer-1.6_results/references/"

## JUICER SUR SCAFFOLDS
## alignement bwa mem (-t 16 -SP5M) d'un SUBSET de PE-reads HiC (24 fichiers fastq = 3 split * R1/R2 * 4 libraries) sur les scaffolds dans ordre pseudo de RENANV2 avec Juicer

## /!\ Pas besoin d'attendre les inter.hic et inter_30.hic de juicer_tools pre, ETAPE QUI NE FONCTIONNE PAS sur le fichier merged_nodups.txt NON MIS EN FORME comme ci-dessous
## cree erreur: "Error: the chromosome combination 1_1 appears in multiple blocks"
## de plus, il faut obtenir la version "candidate assembly" avec 3d-dna pour pouvoir apporter des modifications a l'assemblage dans Juicebox

## obtention d'un fichier merged_nodups.txt avec juicer puis mise en forme awk et sort (separateur " " et non "\t")
awk -F' ' '{if ($2>$6 && $9>0 && $12>0) { print $5,$6,$7,$8,$1,$2,$3,$4,$9,$10,$11,$12,$13,$14,$15,$16 }; if ($2<=$6 && $9>0 && $12>0) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16 } }' /storage/scratch/palasser/juicer-1.6_scaff_sub/RENANV2_scaff_sub/aligned/merged_nodups.txt > /storage/scratch/palasser/juicer-1.6_scaff_sub/RENANV2_scaff_sub/aligned/merged_nodups_RENANV2_scaff_sub_awk.txt
sort -t' ' -k2,2d -k6,6d /storage/scratch/palasser/juicer-1.6_scaff_sub/RENANV2_scaff_sub/aligned/merged_nodups_RENANV2_scaff_sub_awk.txt > /storage/scratch/palasser/juicer-1.6_scaff_sub/RENANV2_scaff_sub/aligned/merged_nodups_RENANV2_scaff_sub_awk_sort.txt

## obtention d'un fichier .assembly avec /3d-dna/utils/generate-assembly-file-from-fasta.awk du pipe 3d-dna
awk -f /storage/scratch/palasser/3d-dna/utils/generate-assembly-file-from-fasta.awk /storage/scratch/palasser/juicer-1.6_scaff_sub/references/RENANV2_SCAFFOLDS_PSEUDO_v1_order.fasta > /storage/scratch/palasser/juicer-1.6_scaff_sub/references/RENANV2_SCAFFOLDS_PSEUDO_v1_order.assembly

## obtention d'un fichier .hic avec /3d-dna/visualize/run-assembly-visualizer.sh du pipe 3d-dna
sbatch -p smp -c 256 --mem=500G --wrap="ml java/oracle-1.8.0_45; export IBM_JAVA_OPTIONS='-Xmx500g -Xgcthreads256'; export _JAVA_OPTIONS='-Xmx500g -Xms500g'; /storage/scratch/palasser/3d-dna/visualize/run-assembly-visualizer.sh /storage/scratch/palasser/juicer-1.6_scaff_sub/references/RENANV2_SCAFFOLDS_PSEUDO_v1_order.assembly /storage/scratch/palasser/juicer-1.6_scaff_sub/RENANV2_scaff_sub/aligned/merged_nodups_RENANV2_scaff_sub_awk_sort.txt"

## JUICEBOX
# corrections apportees au fichier .assembly dans le logiciel interface Juicebox Assembly Tools

## CF. buildPseudomol_RENAN_V2_HiC_reviewed.sh
## =>> construction de pseudomol_v2 de RENANV2 tenant compte des corrections d'assemblage apportees grace aux donnes HiC

## ANALYSE JUICER SUR RENANV2_PSEUDOV2
## alignement bwa mem (-t 16 -SP5M) du meme SUBSET de PE-reads HiC (24 fichiers fastq = 3 split * R1/R2 * 4 libraries) sur les pseudoV2 de RENANV2 avec Juicer
## obtention d'un fichier merged_nodups.txt avec juicer puis mise en forme awk et sort (separateur " " et non "\t")
awk -F' ' '{if ($2>$6 && $9>0 && $12>0) { print $5,$6,$7,$8,$1,$2,$3,$4,$9,$10,$11,$12,$13,$14,$15,$16 }; if ($2<=$6 && $9>0 && $12>0) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16 } }' /storage/scratch/palasser/juicer-1.6/Renan_v2_pseudo_v2/aligned/merged_nodups.txt > /storage/scratch/palasser/juicer-1.6/Renan_v2_pseudo_v2/aligned/merged_nodups_RENANV2_pseudoV2_SUB_mapQ1_awk.txt
sort -t' ' -k2,2d -k6,6d /storage/scratch/palasser/juicer-1.6/Renan_v2_pseudo_v2/aligned/merged_nodups_RENANV2_pseudoV2_SUB_mapQ1_awk.txt > /storage/scratch/palasser/juicer-1.6/Renan_v2_pseudo_v2/aligned/merged_nodups_RENANV2_pseudoV2_SUB_mapQ1_awk_sort.txt
## obtention d'un fichier .hic avec /3d-dna/visualize/run-assembly-visualizer.sh du pipe 3d-dna
awk -f /storage/scratch/palasser/3d-dna/utils/generate-assembly-file-from-fasta.awk /storage/scratch/palasser/juicer-1.6/references/Renan_v13_v2.pseudo.v2.fa > /storage/scratch/palasser/juicer-1.6/Renan_v2_pseudo_v2/aligned/Renan_v13_v2.pseudo.v2.assembly
sbatch -p smp -c 256 --mem=500G --wrap="ml java/oracle-1.8.0_45; export IBM_JAVA_OPTIONS='-Xmx500g -Xgcthreads256'; export _JAVA_OPTIONS='-Xmx500g -Xms500g'; /storage/scratch/palasser/3d-dna/visualize/run-assembly-visualizer.sh /storage/scratch/palasser/juicer-1.6/Renan_v2_pseudo_v2/aligned/Renan_v13_v2.pseudo.v2.assembly /storage/scratch/palasser/juicer-1.6/Renan_v2_pseudo_v2/aligned/merged_nodups_RENANV2_pseudoV2_SUB_mapQ1_awk_sort.txt"

## visualisation pour confirmation des modifications apportees a l'assemblage

# analyses 3d-dna de novo: 3ddna avec 4 rounds sur scaffolds SUB pour un assemblage de novo:
#./run-asm-pipeline.sh -r 4 Triticum_aestivum_RENAN_v2.fasta merged_nodups_SUB_on_RENANV2_scaff_awk_sort.txt 2> 3ddna_on_RENAN_V2_pseudo_sub.log


###############################################################################
#juicer-1.6

"""Usage: juicer.sh [-g genomeID] [-d topDir] [-q queue] [-l long queue] [-s site]
                 [-a about] [-R end] [-S stage] [-p chrom.sizes path]
                 [-y restriction site file] [-z reference genome file]
                 [-C chunk size] [-D Juicer scripts directory]
                 [-Q queue time limit] [-L long queue time limit] [-e] [-h] [-x]
* [genomeID] must be defined in the script, e.g. "hg19" or "mm10" (default
  "hg19"); alternatively, it can be defined using the -z command
* [topDir] is the top level directory (default
  "/Users/nchernia/Downloads/neva-muck/UGER")
     [topDir]/fastq must contain the fastq files
     [topDir]/splits will be created to contain the temporary split files
     [topDir]/aligned will be created for the final alignment
* [queue] is the queue for running alignments (default "short")
* [long queue] is the queue for running longer jobs such as the hic file
  creation (default "long")
* [site] must be defined in the script, e.g.  "HindIII" or "MboI"
  (default "none")
* [about]: enter description of experiment, enclosed in single quotes
* [stage]: must be one of "chimeric", "merge", "dedup", "final", "postproc", or "early".
    -Use "chimeric" when alignments are done but chimeric handling has not finished
    -Use "merge" when alignment has finished but the merged_sort file has not
     yet been created.
    -Use "dedup" when the files have been merged into merged_sort but
     merged_nodups has not yet been created.
    -Use "final" when the reads have been deduped into merged_nodups but the
     final stats and hic files have not yet been created.
    -Use "postproc" when the hic files have been created and only
     postprocessing feature annotation remains to be completed.
    -Use "early" for an early exit, before the final creation of the stats and
     hic files
* [chrom.sizes path]: enter path for chrom.sizes file
* [restriction site file]: enter path for restriction site file (locations of
  restriction sites in genome; can be generated with the script
  (misc/generate_site_positions.py) )
* [reference genome file]: enter path for reference sequence file, BWA index
  files must be in same directory
* [chunk size]: number of lines in split files, must be multiple of 4
  (default 90000000, which equals 22.5 million reads)
* [Juicer scripts directory]: set the Juicer directory,
  which should have scripts/ references/ and restriction_sites/ underneath it
  (default /broad/aidenlab)
* [queue time limit]: time limit for queue, i.e. -W 12:00 is 12 hours
  (default 1200)
* [long queue time limit]: time limit for long queue, i.e. -W 168:00 is one week
  (default 3600)
* -f: include fragment-delimited maps from hic file creation
* -e: early exit
* -h: print this help and exit"""

#successives jobs in juicer.sh:
#head-
#split-

#-S "chimeric"
#count_ligation-
#align1-
#merge-
#aligncheck-

#-S "merge"
#fragmerge-

#-S "dedup"
#dedupguard-
#dedup-
#post_dedup-
#dupcheck-

# "mnd" merge_nodup.txt created between those steps
# 3D DNA has to be used for a draft genome
#Olga Dudchenko (2018):
#For draft genome, use 3D-DNA visualize/run-assembly-visualizer.sh to build .hic files, not juicer.

#-S "final"
#prestats-
#stats-
#stats30-
#fincln1-
#hic-
#hic30-
#hiccups_wrap-
#arrowhead_wrap-
#fincln-


#########################################################
# INFO sur l'alignement bwa mem avec option -SP5M:
#The -SP option is used to ensure the results are equivalent to that obtained by running bwa mem on each mate separately, while retaining the right formatting for paired-end reads.
#This option skips a step in bwa mem that forces alignment of a poorly aligned read given an alignment of its mate with the assumption that the two mates are part of a single genomic segment.

#The -5 option is used to report the 5' portion of chimeric alignments as the primary alignment. In Hi-C experiments, when a mate has chimeric alignments,
#typically, the 5' portion is the position of interest, while the 3' portion represents the same fragment as the mate. For chimeric alignments, bwa mem reports two alignments:
#one of them is annotated as primary and soft-clipped, retaining the full-length of the original sequence. The other end is annotated as hard-clipped and marked as either 'supplementary' or 'secondary'.
#The -5 option forces the 5'end to be always annotated as primary.

#The -M option is used to annotate the secondary/supplementary clipped reads as secondary rather than supplementary, for compatibility with some public software tools such as picard MarkDuplicates.


#########################################################
#NB: processed count files format:
#name str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 mapq1 mapq2
#Each line represents one read pair.
#str1/str2: strand of read end 1 and read end 2, forward or revere, 0 is forward
#chr1/chr2: chromosome of read end 1 and read end 2
#pos1/pos2: position of read end 1 and read end 2
#frag1/frag2: fragment number of read end 1 and read end 2 (depends on restriction site)
#mapq1/mapq2: mapping quality score of read end 1 and read end 2; this is ignored if the -q flag is not set but will one day be used for quality control metrics
