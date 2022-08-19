#!/bin/bash
#SBATCH --job-name=hicmano 
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --mem=90G
#SBATCH -p debug
#SBATCH --export=ALL

echo "Start on `hostname`"

INPUT='/home/palasser/results/juicer-1.6_results/RENANV2/aligned'

#filtre sur mapQ>1, filtre sur PE intra scaff, filtre sur PE intra fragment de restriction
cat $INPUT/merged_nodups.txt |tr ' ' '\t' |cut -f2,4,6,8,9,12 |gawk '{if ($5>1 && $6>1) {print $0} }' |gawk '{if ($1!=$3 && $2!=$4) {print $1,$3,$5,$6} }' > $INPUT/merged_nodups_mapQ2_scaff_paires.txt

#script perl pour compter le nb de paires de scaff + tri
cat $INPUT/merged_nodups_mapQ2_scaff_paires.txt |tr ' ' '\t'|gawk '{if ($3>=15 && $4>=15) {print $0} }' |cut -f1,2 |perl -ne 'chomp; $h{$_}++; END { foreach $k(sort keys %h){ print $k,"\t",$h{$k},"\n"; }}' |sort -k1,1 -k3,3nr > $INPUT/hic_RENAN_V2_mapQ15_scaff_paires.txt

cat $INPUT/merged_nodups_mapQ2_scaff_paires.txt |tr ' ' '\t'|gawk '{if ($3>=30 && $4>=30) {print $0} }' |cut -f1,2 |perl -ne 'chomp; $h{$_}++; END { foreach $k(sort keys %h){ print $k,"\t",$h{$k},"\n"; }}' |sort -k1,1 -k3,3nr > $INPUT/hic_RENAN_V2_mapQ30_scaff_paires.txt

#merge_nodups.txt format
# <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <cigar1> <sequence1> <mapq2> <cigar2> <sequence2> <readname1> <readname2>
# str = strand (0 for forward, anything else for reverse)
# chr = chromosome (must be a chromosome in the genome)
# pos = position
# frag = restriction site fragment
# mapq = mapping quality score
# cigar = cigar string as reported by aligner
# sequence = DNA sequence

