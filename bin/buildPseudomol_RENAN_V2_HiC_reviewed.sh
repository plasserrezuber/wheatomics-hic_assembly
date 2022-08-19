#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --job-name=buildPseudo
#SBATCH --partition=fast

#a mettre dans le .bashrc: alias awk='awk -v OFS="\t"'
# rappel: join -t$'\t'

#################################################################################################################################
##### mise en forme fichier .assembly format type BDD: correspondance debut fichier ('>'Noms des scaff, leurs num arbitraires)
##### et fin fichier (ordre scaff corrige grace aux data HiC):
# /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order.txt:
# colonne 1: num arbitraire donne a chaque scaffold par 3d-dna/run-assembly-visualizer.sh dans fichier .assembly + incrementation Juicebox qd scaff fragmente
# servira de cle de correspondance entre les deux tables debut/fin du fichier .assembly
# colonne 2: orientation strand 1/-1 modifie dans Juicebox (tous les scaff sont strand 1 avant modif, fasta ref = pseudomolecules dont scaff sont issus ici)
# colonne 3: ajout d'un numero avec la commande gawk -v OFS='\t' '{ if ($1<0) {print $1,"-1",NR} else {print $1,"1",NR} }'
# pour conserver l'odre correct expertise et corrige dans Juicebox grace aux data HiC
# colonne 4: nom des scaffolds

## mise en forme fin fichier assembly
grep -v '>' /home/palasser/results/3d-dna/3d-dna_visualize_candidate_assembly/RENANV2_SCAFFOLDS_PSEUDO_v1_order.4mars_review_FINAL.assembly |sed -E 's/ /\n/g' \
|gawk -v OFS='\t' '{ if ($1<0) {print $1,"-1",NR} else {print $1,"1",NR} }' |sed 's/^-//' > /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order_tmp.txt

## mise en forme debut fichier assembly et join du debut et de la fin
join -1 1 -2 2 <(sort -n -k1,1 /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order_tmp.txt) \
<(grep '>' /home/palasser/results/3d-dna/3d-dna_visualize_candidate_assembly/RENANV2_SCAFFOLDS_PSEUDO_v1_order.4mars_review_FINAL.assembly |sed 's/>//' |tr ' ' '\t' |cut -f1,2,3 |sort -n -k2,2) \
|sort -n -k3,3 > /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order.txt

#################################################################################################################################
##### join des corrections apportees grace aux data HiC et des info de construction des PSEUDOV1 (fichier csv: scaff, length, chrom, strand)
# travail sur les scaff fragmentes
grep -v 'debris' /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order.txt |grep 'fragment' |gawk -v OFS='\t' '{ print $4,$2,$3,$5 }' |sed 's/:::/\t/' \
|sed 's/fragment_3/fragment_2/' |sed 's/fragment_5/fragment_3/' |sed 's/fragment//' |gawk -v OFS='\t' '{print $1,$1$2,$3,$4,$5}' \
> /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order_fragmented_tmp.txt

join -1 1 -2 1 <(sort -k1,1 /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order_fragmented_tmp.txt) \
<(sort -k1,1 /home/palasser/data/RENAN_v2_pseudo/scaffAssignment.50000.csv.order.orientation.csv) |gawk -v OFS='\t' '{ print $4,$2,$3,$5,$7,$8 }' \
> /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order_fragmented.txt

#################################################################################################################################
##### join des corrections apportees grace aux data HiC et des info de construction des PSEUDOV1 (fichier csv: scaff, length, chrom, strand)
# travail sur les scaff non fragmentes (le join exclut les scaffolds fragmentes de fait, n'existant pas dans fichier .csv)
join -1 1 -2 1 <(gawk -v OFS='\t' '{ print $4,$2,$3 }' /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order.txt |sort -k1,1) \
<(sort -k1,1 /home/palasser/data/RENAN_v2_pseudo/scaffAssignment.50000.csv.order.orientation.csv) |gawk -v OFS='\t' '{ print $3,$1,$2,$4,$5,$6 }' \
> /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order_NOTfragmented.txt

##### VERIF
# grep '>' /home/palasser/results/3d-dna/3d-dna_visualize_candidate_assembly/RENANV2_SCAFFOLDS_PSEUDO_v1_order.4mars_review_FINAL.assembly \
# |cut -d'_' -f4 |sort -n |uniq |wc -l
## 2566 scaffolds dans le fichier .assembly

# cut -f2 /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order_fragmented.txt |cut -d'_' -f4 |sort -n |uniq |wc -l
## 18 scaffolds fragmentes
# wc -l /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order_fragmented.txt
## 38 fragments de scaffolds

# wc -l /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order_NOTfragmented.txt
## 2548 scaffolds non fragmentes
## (2548 + 18 = 2566)

#################################################################################################################################
##### mise en forme fichier final ordre et strand des scaffolds dans PSEUDOV2

# info scaff fragmentes et non fragmentes dans un seul fichier, avec comme orientation finale (derniere colonne):
# si -1 -1 => 1 ; si -1 1 => -1 ; si 1 => strand PSEUDOV1 conserve
cat /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order_fragmented.txt /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order_NOTfragmented.txt \
|sort -n -k1,1 |gawk -v OFS='\t' '{ if ($3==-1 && $6==-1) {print $0,"1"} if ($3==-1 && $6==1) {print $0,"-1"} if ($3==1) {print $0,$6} }' \
> /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER_info_tmp.txt

rm /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER_info.txt

## header
echo -e "order\t#scaff\tjuicebox_strand\tlength\tPSEUDOV2_chrom\tPSEUDOV1_strand\tPSEUDOV2_strand" > /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER_info.txt
## MISE A JOUR de PSEUDOV1_chrom en PSEUDOV2_chrom ($8=chrom_before $9=chrom_after)
paste <(awk -v OFS="\t" 'NR!=1 {chrom_before=a[5]} {split($0,a,FS)} {print $0,chrom_before}' /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER_info_tmp.txt) \
<(awk 'f{print $5;f=0} /chr/{f=1}' /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER_info_tmp.txt) \
|gawk -v OFS='\t' '{ if ($8==$9) {print $1,$2,$3,$4,$8,$6,$7} else {print $1,$2,$3,$4,$5,$6,$7} }' \
>> /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER_info.txt

### /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
## §§§§ !!!!!!!!!  §§§§§§§§§ MISE A JOUR NECESSAIRE A LA MAIN QUAND PLUSIEURS SCAFF A LA SUITE SONT ASSIGNES A UN CHROM DIFFERENT

####################################################################
### Creation fichier bed des scaffolds
for chr in "1A" "1B" "1D" "2A" "2B" "2D" "3A" "3B" "3D" "4A" "4B" "4D" "5A" "5B" "5D" "6A" "6B" "6D" "7A" "7B" "7D";
do
  gawk -v chrom=$chr '{ if ($5~chrom) { sum+=$4; count+=1; print $5,sum,$2,(count-1)*100,count,$4,$7} }' /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER_info.txt
done |gawk -v OFS="\t" 'NR!=1 {chrom_before=a[1]; end_before=a[2]} {split($0,a,FS)} { if ($1==chrom_before) {print $1,end_before+$4,$2+$4,$3,$5,$6,$7} else {print $1,0,$2,$3,$5,$6,$7} }' \
|gawk -v OFS="\t" '{ if ($7==-1) {print $1,$2,$3,$4,$5,$6,"-"} else {print $1,$2,$3,$4,$5,$6,"+"} }' > /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0_scaffold.bed



##### Nb de changements de chromosome entre PSEUDOV1 et PSEUDOV2
# paste <(awk -v OFS="\t" 'NR!=1 {chrom_before=a[5]} {split($0,a,FS)} {print $0,chrom_before}' /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER_info_tmp.txt) \
# <(awk 'f{print $5;f=0} /chr/{f=1}' /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER_info_tmp.txt) \
# |gawk -v OFS='\t' '{ if ($5!=$8 || $5!=$9) {print} }'
#comptage a la main: 16

##### VERIF
# grep -c 'Taestivum_RENAN_scaffold_' /home/palasser/data/RENAN_v2_pseudo/scaffAssignment.50000.csv.order.orientation.csv
## 2566 scaffolds dans les pseudomolecules v1

# cut -f2 /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER.txt |grep -v '#scaff' |cut -d'_' -f4 |sort -n |uniq |wc -l
## 2566 scaffolds a integrer dans les pseudomolecules v2, dont 18 decoupes en 38 fragments

# grep -c -v '#' /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER_info.txt
## 2586 (= 2548 + 38)

##### COMPARAISON TAILE CUMULEE PSEUDOV versus PSEUDOV2
# gawk '{ print sum+=$2 }' /home/palasser/data/RENAN_v2_pseudo/scaffAssignment.50000.csv.order.orientation.csv |tail -n1
## 14198186086
# gawk '{ print sum+=$4 }' /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER_info.txt |tail -n1
## 14198136386
# echo "$((14198186086-14198136386))"
## 49700
# grep 'debris' /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order.txt |gawk '{ print sum+=$5 }' |tail -n1
## 49700

#################################################################################################################################
##### fasta pour chaque sacffold a partir de l'assemblage v13/v2 fourni par le Genoscope
mkdir /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fasta
ml exonerate
fastaexplode /storage/groups/gdec/shared/triticum_aestivum/wheatomics/renan/assembly/v2/Triticum_aestivum_RENAN_v2.fasta \
-d /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fasta

#################################################################################################################################
##### fichier decrivant les gaps presents dans les scaffolds
ml bioperl gdecTools
for i in $(cut -f2 /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order_fragmented.txt |cut -d'_' -f4 |sort -n |uniq);
do
  #last sed to supress empty lines
  nCount.pl /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fasta/Taestivum_RENAN_scaffold_${i}.fa |sed 's/://' |sed 's/x//' |sed 's/\.\./ /' |sed '/^$/d' \
  |cut -d' ' -f1,2,3,4 |tr ' ' '\t' |grep -v 'Taes' > /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fasta/Taestivum_RENAN_scaffold_${i}.fa_Ngap.txt
done

#################################################################################################################################
##### creation INPUT pour script closestGAP.py contenant les coord des fragments sur fasta source de l'assemblage fourni par Genoscope
# but: redefinir les coordonnees des fragments de scaff en utilisant le Ngap le plus proche, extraire la seq avec subseq.pl a partir du fasta source
join -t$'\t' -1 1 -2 1 <(grep 'fragment' /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_4mars_HiC_reviewed_order.txt |gawk -v OFS='\t' '{ print $4,$2,$3,$5 }' |sed 's/:::/\t/' |sed 's/:::debris//' |sort -k1,1 -k2,2) \
<(cut -f1,2,4 /home/palasser/data/RENAN_v2_pseudo/scaffAssignment.50000.csv.order.orientation.csv |sort -k1,1 ) \
> /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fragments_coordinates_tmp.txt

rm /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fragments_coordinates.txt

# 1er cmd gawk: cumul tailles -> coord end
# 2eme cmd gawk: split the line to an array and when processing the next line we use the values stored of $8+1 from previous line (coord end) -> coord start of following line
# 3eme cmd gawk: mise en forme colonne start/stop pour fragments non revcom ($7==1) et recalcul des coord sur fasta source pour fragments revcom ($7==-1)
# 4eme cmd gawk: suppression lignes "debris"
# derniere ligne de cmd: mise en forme tel que tableau final pour buildPseudomol.pl
for s in $(cut -f1 /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fragments_coordinates_tmp.txt |sort |uniq);
do
  grep -w $s /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fragments_coordinates_tmp.txt \
  |gawk -v OFS='\t' '{sum+=$5; print $0,sum,1}' |gawk -v OFS='\t' 'NR!=1 { $9=a[8]+1 } { split($0,a,FS) } {print}' \
  |gawk -v OFS='\t' '{ if ($7==1) print $1,$2,$3,$4,$5,$6,$7,$9,$8; if ($7==-1) print $1,$2,$3,$4,$5,$6,$7,$6-$8+1,$6-$9+1 }' \
  |gawk -v OFS='\t' '{ if ($5>5000) print $0 }'\
  |sed 's/fragment_3/fragment_2/' |sed 's/fragment_5/fragment_3/' |sed 's/fragment//' |gawk -v OFS='\t' '{print $1,$1$2,$6,$8,$9}' \
  >> /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fragments_coordinates.txt
done

rm /home/palasser/data/RENAN_v2_pseudo/*_tmp.txt

#################################################################################################################################
##### redefinition des coordonnees start/end des fragment avec closestGAP.py
ml python
/home/palasser/bin/closestGAP.py -i /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fragments_coordinates.txt \
-d /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fasta/ -s .fa_Ngap.txt

#################################################################################################################################
##### MISE A JOUR de scaff length dans /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER_info.txt pour les scaffolds fragmentes
for f in $(cut -f2 /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fragments_coordinates_closestGAP.out);
do
  length=$(grep $f RENANV2_SCAFFOLDS_fragments_coordinates_closestGAP.out |gawk '{ print $4-$3+1 }')
  gawk -v frag=$f -v ln=$length '{ if ($2==frag) {print $1"\t"$2"\t"$3"\t"ln"\t"$5"\t"$6"\t"$7} else {print $0} }' /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER_info.txt \
  > /home/palasser/data/RENAN_v2_pseudo/tmp && mv /home/palasser/data/RENAN_v2_pseudo/tmp /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER_info.txt
done

#################################################################################################################################
##### subseq.pl
ml bioperl gdecTools
for f in $(cut -f2 /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fragments_coordinates_closestGAP.out);
do
  s=$(grep $f /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fragments_coordinates_closestGAP.out |cut -f1)
  start=$(grep $f /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fragments_coordinates_closestGAP.out |cut -f3)
  end=$(grep $f /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fragments_coordinates_closestGAP.out |cut -f4)

  subseq.pl -if FASTA /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fasta/${s}.fa -s $start -e $end -of FASTA > /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fasta/${f}.fa
  sed -i "s/>${s}/>${f}/" /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fasta/${f}.fa
done

##### VERIF
# ls RENANV2_SCAFFOLDS_fasta/Taestivum_RENAN_scaffold_*_*.fa |wc -l
## 38 fichiers .fa

#################################################################################################################################
##### buildPseudomol.pl

## fichier input seqOrderFile pour buildPseudomol.pl
gawk -v OFS='\t' '{ print $5,$2,$7 }' /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER_info.txt |grep -v '#' \
> /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER.txt

rm /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fasta/Triticum_aestivum_RENAN_v2_to_buildPseudomol.fa

## fichier input multifasta des sacffolds a assembler dans les pseudo v2 par buildPseudomol.pl (n'a pas a respecter d'ordre particulier)
for s in $(cut -f2 /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER.txt);
do
  cat /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fasta/${s}.fa >> /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fasta/Triticum_aestivum_RENAN_v2_to_buildPseudomol.fa
done

ml bioperl gdecTools
buildPseudomol.pl -g 100 -s /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER.txt \
/home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fasta/Triticum_aestivum_RENAN_v2_to_buildPseudomol.fa \
> /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0.fa

## creation .fai
ml samtools
samtools faidx /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0.fa

#################################################################################################################################
## multi fasta avec tous les scaffolds non assignes a une pseudo: chromosome unknown
for s in $(cut -f2 /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER.txt |sed -E 's/(Taestivum_RENAN_scaffold_[0-9]*).*/\1/');
do
  rm /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fasta/${s}.fa
done
for s in $(cut -f2 /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER.txt);
do
  rm /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fasta/${s}.fa
done

##### VERIF
# find Tae*.fa |wc -l
## 338
# echo $((2566+338))
## 2904
# wc -l /home/celmonat/Rn_v13_v2_scf/Triticum_aestivum_RENAN_v2.fasta.fai
## 2904 

#################################################################################################################################
## Fichier fasta des scaff non ancres
unanch_scaff=$(\ls -1 /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fasta/Taestivum_RENAN_scaffold_*.fa |sort -V |tr -s '\n' ' ')
cat $unanch_scaff >> /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0_unanchored_scaffolds.fa
# grep '>' TaeRenan_refseq_v2.0_unanchored_scaffolds.fa |wc -l
## 338

#################################################################################################################################
## fichier tsv avec liste des scaffolds, assignation aux chrom, ordre, orientation
rm /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0_scaffold_order.tsv
## liste des scaff dans les pseudo
cat /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER.txt > /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0_scaffold_order.tsv
## ajout des scaff unanchored
grep '>' /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0_unanchored_scaffolds.fa |sed 's/>//' |awk -v OFS='\t' '{print "chrUn",$0,1}' \
>> /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0_scaffold_order.tsv

##### VERIF
## nb scaffolds avec nom de depart dans assemblage Renan_v13_v2 fourni par le Genoscope en novembre 2020
# cut -f2 /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0_scaffold_order.tsv |sed -E 's/(Taestivum_RENAN_scaffold_[0-9]*).*/\1/' |sort |uniq |wc -l
##2904
## nb de scaffolds fragmentes 
# cut -f2 /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0_scaffold_order.tsv |grep -E 'Taestivum_RENAN_scaffold_[0-9]*_[0-9]' |wc -l
## 38

#################################################################################################################################
##### VERIF taille cumulee des scaff par chrom avant/apres buildPseudomol.pl => OK 
##tailles chrom PSEUDOV2 avant buildPseudomol.pl (cumul scaff length + GAP 100N)
# grep -v '#' /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER_info.txt |gawk -v OFS='\t' '{ a[$5]+=$4+100 } END{for (i in a) print i, a[i]-100 }' |sort -k1,1 \
# > /home/palasser/data/RENAN_v2_pseudo/scaff_cumul_size_per_chrom_pseudo_v2.txt
# diff /home/palasser/data/RENAN_v2_pseudo/scaff_cumul_size_per_chrom_pseudo_v2.txt <(cut -f1,2 /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0.fa.fai)

#################################################################################################################################
##### VERIF PSEUDOV1 versus PSEUDO V2  => OK :))
##### COMPARAISON TAILLE CUMULEE DES SCAFF PAR CHROM PSEUDOV1 versus PSEUDOV2 avant contruction des 21 pseudo

# gawk '{ print sum+=$2 }' /home/palasser/data/RENAN_v2_pseudo/scaffAssignment.50000.csv.order.orientation.csv |tail -n1
## 14198186086
# gawk '{ print sum+=$4 }' /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_HiC_REVIEWED_ORDER_info.txt |tail -n1
## 14195387115

## calcul taille cumulee des gaps retires des scaffolds fragmentes suite a closestGAP.py
# /home/palasser/bin/closestGAP.py -i /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fragments_coordinates.txt -s .fa_Ngap.txt \
# -d /home/palasser/data/RENAN_v2_pseudo/RENANV2_SCAFFOLDS_fasta |grep ',' |sed -E 's/\[|\]|,//g' |sed 's/'$"'"'//g' |uniq |gawk -v OFS=' ' '{ print sum+=$2 }' |tail -n1
## 2798971

## difference "taille genome" PSEUDOV1 - PSEUDOV2 - taille cumulee des gaps retires (debris ne sont plus a prendre en compte)
# echo $((14198186086-14195387115-2798971))
## 0

#################################################################################################################################
## TAILLE PSEUDOV1: taille cumulee des 21 pseudo
# gawk '{ print sum+=$2 }' /home/celmonat/Renan_v13_v2.pseudo.v1.fa.fai |tail -n1
## 14200731086

## TAILLE PSEUDOV2: taille cumulee des 21 pseudo
# gawk '{ print sum+=$2 }' /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0.fa.fai |tail -n1
# 14195643615
## diff TAILLE PSEUDOV1 - TAILLE PSEUDOV2 - taille cumulee des gaps retires - taille cumulee des gap 100N ajoutes
# echo $((14200731086-14195643615-2798971-2586*100+21*100))
## 2032000


######### PARTAGE FICHIERS AVEC EQUIPE GENOSCOPE
## ex avec le fichier fasta du genome: connexion au serveur 7serv dans terminal de MobaXterm
## scp de hpc2 vers 7serv pour rapatrier les fichiers:
# scp hpc2:/home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0.fa /home/plasserre/data/RENAN_v2_pseudo/
## puis toujours connecte au 7serv, rsync de 7serv vers cobalt
# rsync -auv /home/plasserre/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0.fa lasserrp@cobalt.ccc.cea.fr:/ccc/store/cont007/fg0118/fg0118/TaeRenan_refseq_v2.0/TaeRenan_refseq_v2.0.fa

## transfert fichiers magatt d'Helene:
#scp hpc2:/home/herimbert/work/wheatomics/wp1-renan/annot/pseudo_v2/magatt_triannot_merged/TaeRenan_refseq_v2.0_annotation.tar.gz /home/plasserre/data/RENAN_v2_pseudo/
#rsync -auv /home/plasserre/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0_annotation.tar.gz lasserrp@cobalt.ccc.cea.fr:/ccc/store/cont007/fg0118/fg0118/TaeRenan_refseq_v2.0/annotation/TaeRenan_refseq_v2.0_annotation.tar.gz


#################################################################################################################################
## FICHIER AGP
#https://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp/agp_validate.cgi

echo -e "##agp-version 2.1\n#ORGANISM: Triticum aestivum; cv. Renan\n#ASSEMBLY NAME: TaeRenan_refseq_v2.0\n#ASSEMBLY DATE: 22-March-2021" > /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0.agp

for chr in "1A" "1B" "1D" "2A" "2B" "2D" "3A" "3B" "3D" "4A" "4B" "4D" "5A" "5B" "5D" "6A" "6B" "6D" "7A" "7B" "7D";
do
  end_chrom=$(grep $chr /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0.fa.fai |cut -f2)
  gawk -v OFS='\t' -v chrom=$chr -v nd_chrm=$end_chrom '{ if ($1~chrom && $3<nd_chrm) {print $1,$2+1,$3,"W",$4,"1",$6,$7"\n"$1,".",".","U",100,"scaffold","yes","align_genus;proximity_ligation"}; if ($1~chrom && $3==nd_chrm) {print $1,$2+1,$3,"W",$4,"1",$6,$7} }' /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0_scaffold.bed \
  |gawk -v OFS='\t' '{ print $1,$2,$3,NR,$4,$5,$6,$7,$8 }'
done |gawk -v OFS='\t' -v chrom=$chr -v nd_chrm=$end_chrom 'NR!=1 {end_before=a[3]} {split($0,a,FS)} { if ($5=="U") {print $1,end_before+1,end_before+100,$4,$5,$6,$7,$8,$9} else {print $0} }' \
>> /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0.agp

#################################################################################################################################
##### VERIF taille chrom identique taille chrom .fai
## AGP ok
diff <(cut -f2 /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0.fa.fai) \
<(for chr in "1A" "1B" "1D" "2A" "2B" "2D" "3A" "3B" "3D" "4A" "4B" "4D" "5A" "5B" "5D" "6A" "6B" "6D" "7A" "7B" "7D";
do grep $chr /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0_AGP.tsv |cut -f3 |tail -n1; done)

##bed ok
diff <(cut -f2 /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0.fa.fai) \
<(for chr in "1A" "1B" "1D" "2A" "2B" "2D" "3A" "3B" "3D" "4A" "4B" "4D" "5A" "5B" "5D" "6A" "6B" "6D" "7A" "7B" "7D";
do grep $chr /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0_scaffold.bed |cut -f3 |tail -n1; done)