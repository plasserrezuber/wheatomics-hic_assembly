# Description  

## Work directory:
dir=/home/fchoulet/data/projects/wheatomics/ontRenan/v11/qa/isbp/
Quality assignment and pseudomolecule construction based on CS-ISBP BLAST versus Renan_v11
5M ISBPs from CS were BLASTED against Renan_v11 by Cécile Monat
Results were filtered: 1 HIT 1 HSP >98%ID >90%ov
file=/home/cmonat/Renan_v11/OK_ISBP_99_ov_max_no_N_no_MAPQ_equal_0_CS_sorted_vs_Renan_v11_i80_qov90_mHSP1_mTS11_bestHit
Only_iSup98.bltn
We want to use these ISBPs to order and orientate (and check chimeric Renan scaff)
Development of assigner.pl
_ The script analyze ISBP matches between CS and Renan
_ Foreach scaffold -> order the matches according to CS-chr position of ISBPs
_ Consider pairs of adjacent ISBPs (that have matched)
_ Calculate the distance separating neighbor ISBPs on CS + on Renan scaff
_ Use a distance threshold (default=50kb) to consider a valid pair
_ Classify valid versus nonvalid pairs
_ Count, for each scaff, the % of valid pairs by chromosome (ex: 90% on chr1A, 10% on chr2A)
_ Calculate the median position of CS-ISBPs that are within valid pairs
_ Determine the orientation of Renan scaff based on majority rule of ISBP matches (+/+ or +/-)
_ Print data in CSV format
I tried 3 distance thresholds: 50kb 100kb 200kb (to consider a valid ISBP pair)
for i in 50000 100000 200000 ; do ./assigner.pl -d $i $file > scaffAssignment.$i.csv 2>scaffAssignment.$i.log ; done
Example for scaffold_1
#scaff length chr nb_valid_pairs percent_valid_pairs median ori ori%
Taestivum_RENAN_scaffold_1 251657134 chr2A 73075 99 571205030 1 90
Taestivum_RENAN_scaffold_1 251657134 chr2D 75 0 439979109 1 85
Taestivum_RENAN_scaffold_1 251657134 chr2B 15 0 555845586 1 100
Taestivum_RENAN_scaffold_1 251657134 chr3A 14 0 557616969 1 50
Taestivum_RENAN_scaffold_1 251657134 chr1A 9 0 159337445 1 77
Taestivum_RENAN_scaffold_1 251657134 chr4B 6 0 66595872 -1 83
Taestivum_RENAN_scaffold_1 251657134 chr4A 5 0 165046690 -1 80
Taestivum_RENAN_scaffold_1 251657134 chr5A 4 0 566194195 1 100
Taestivum_RENAN_scaffold_1 251657134 chr6B 3 0 48530366 -1 66
Taestivum_RENAN_scaffold_1 251657134 chr6A 3 0 422941003 1 66
Taestivum_RENAN_scaffold_1 251657134 chr7B 3 0 201577472 -1 100
Taestivum_RENAN_scaffold_1 251657134 chr7A 2 0 9507074 1 50
Taestivum_RENAN_scaffold_1 251657134 chr6D 2 0 1213543 -1 100
Taestivum_RENAN_scaffold_1 251657134 chr4D 1 0 30691766 1 100
Taestivum_RENAN_scaffold_1 251657134 chr3D 1 0 239096369 1 100
Taestivum_RENAN_scaffold_1 251657134 chr3B 1 0 731535664 1 100
ex: 99% of valid ISBP pairs (<50kb) fall into chr2A - orientation=1 (+ strand) because 90% of valid ISBPs have
matches on + strand

# Consider a minimum % of valid pairs to assign a scaffold to a chromosome
# with >50%
# + order the scaffolds based on the median value
```bash
for i in 50000 100000 200000 ; do cat scaffAssignment.$i.csv | awk '$5>50 {print}' | sort -k3,3 -k6,6n | cut -f1,2,3,7
> scaffAssignment.$i.csv.order.orientation.csv; done
#scaff length chr ori
Taestivum_RENAN_scaffold_630 2084696 chr1A -1
Taestivum_RENAN_scaffold_544 3582496 chr1A 1
Taestivum_RENAN_scaffold_1210 107967 chr1A 1
Taestivum_RENAN_scaffold_409 7116397 chr1A 1
Taestivum_RENAN_scaffold_1680 86333 chr1A -1
Taestivum_RENAN_scaffold_3198 37844 chr1A 1
Taestivum_RENAN_scaffold_415 6884954 chr1A 1
Taestivum_RENAN_scaffold_596 2631102 chr1A -1
Taestivum_RENAN_scaffold_833 310854 chr1A 1
wc -l *csv
2964 scaffAssignment.50000.csv.order.orientation.csv
2967 scaffAssignment.100000.csv.order.orientation.csv
2974 scaffAssignment.200000.csv.order.orientation.csv
```
(/!\ there is a header)
=> The number of assigned scaffolds are almost the same whatever the threshold used
=> KEEP data with 50kb
2963 scaffolds are part of pseudomolecules
14,091,179,560 bp assigned to chromosomes / 14,448,184,776 bp in total = 97.5% anchored

## Analyze chimeric scaffolds  
```bash
diff <(cut -f1 scaffAssignment.50000.csv|egrep -v '^#'|sort -u) <(cut -f1
scaffAssignment.50000.csv.order.orientation.csv |egrep -v '^#'|sort -u) | egrep '<' | cut -d ' ' -f2|sort -u |wc -l
#264 scaffolds remained unassigned => 110,227,146 bps
#= 167 scaffolds with no_valid_pairs
#+ 97 scaffolds with valid pairs but on several chromosomes
```

## Get length of chimeric scaffolds:  
```bash
for scaff in $(diff <(cut -f1 scaffAssignment.50000.csv|egrep -v '^#'|sort -u) <(cut -f1
scaffAssignment.50000.csv.order.orientation.csv |egrep -v '^#'|sort -u) | egrep '<' | cut -d ' ' -f2|sort -u ); do
egrep -w $scaff scaffAssignment.50000.csv ; done | egrep -v no_valid | sort -k2,2rn | cut -f1,2 |uniq | head
#+ 1 scaff of 37 Mb (scaffold_100)
#+ 10 scaff betw 1 and 10 Mb
#+ all others are <1 Mb
Taestivum_RENAN_scaffold_100 37394307
Taestivum_RENAN_scaffold_323 10113916
Taestivum_RENAN_scaffold_383 7932307
Taestivum_RENAN_scaffold_402 7420997
Taestivum_RENAN_scaffold_456 5859505
Taestivum_RENAN_scaffold_499 4481595
Taestivum_RENAN_scaffold_509 4225830
Taestivum_RENAN_scaffold_615 2306412
Taestivum_RENAN_scaffold_679 1533960
Taestivum_RENAN_scaffold_700 1245812
Taestivum_RENAN_scaffold_713 1129245
```
Length of pseudomolecules compared to CS RefSeq_V2.1:

```bash
paste <(for c in $chr; do len=$(egrep -w $c scaffAssignment.50000.csv.order.orientation.csv | sumcalc.pl -f 2|cut -f1);
echo $c $len ; done) <(cat /home/share/iwgsc_refseq_v2.1/CS_pesudo_v2.1.fa.fai) | awk '{OFS="\t"; print
$1,$2,$4,int(($4-$2)/1000000)}'
#chr len_RENAN len_CS diff_Mb
chr1A 590924755 598660471 7
chr1B 697078618 700547350 3
chr1D 491806405 498638509 6
chr2A 775407558 787782082 12
chr2B 804704093 812755788 8
chr2D 630263957 656544405 26
chr3A 750783825 754128162 3
chr3B 854567599 851934019 -2
chr3D 624053163 619618552 -4
chr4A 745415580 754227511 8
chr4B 665160430 673810255 8
chr4D 516965750 518332611 1
chr5A 708451606 713360525 4
chr5B 730055536 714805278 -15
chr5D 533718974 569951140 36
chr6A 623358468 622669697 0
chr6B 717330016 731188232 13
chr6D 491372808 495380293 4
chr7A 745542366 744491536 -1
chr7B 748687958 764081788 15
chr7D 645530095 642921167 -2
```

# Create PSEUDOMOLECULES  

## EXTRACT scaff per chr + revcom when ori=-1  
```bash
db=/home/share/wheatomics/renan/assemblies/Triticum_aestivum_RENAN_v1.fasta
for c in $chr ; do if [[ $c == "chrUn" ]]; then continue; fi; egrep -w $c
scaffAssignment.50000.csv.order.orientation.csv | cut -f1,4 | while IFS=$'\t' read -r scaff ori; do blastdbcmd -entry
$scaff -db $db > indivScaff/$c/$scaff.fa; if [[ $ori == -1 ]]; then revcom.pl indivScaff/$c/$scaff.fa >
indivScaff/$c/$scaff.fa.rc && mv indivScaff/$c/$scaff.fa.rc indivScaff/$c/$scaff.fa ; fi ; done; done
## CREATE multi fasta of scaffolds in the correct order + orientation => and use scaffolder.pl  
for c in $chr ; do if [[ $c == "chrUn" ]]; then continue; fi; for scaff in $(egrep -w $c
scaffAssignment.50000.csv.order.orientation.csv | cut -f1); do cat indivScaff/$c/$scaff.fa >> indivScaff/$c.scaff.fa ;
done; done
mkdir pseudomolecules
for c in $chr ; do if [[ $c == "chrUn" ]]; then continue; fi; scaffolder.pl -g 1000 indivScaff/$c.scaff.fa >
pseudomolecules/$c.fa ; done
rm indivScaff/*.scaff.fa
cd pseudomolecules/
renameId.pl -r -ne *.fa
cat *.fa > renan.pseudo.v1.fa
makeblastdb -in renan.pseudo.v1.fa -dbtype nucl -parse_seqids
samtools faidx renan.pseudo.v1.fa
rm chr*.fa
```

## CREATE GFF files for scaffolds and gaps + GFF with lengths divided by 10 (for ACT visualization)  
```bash
for c in $chr ; do if [[ $c == "chrUn" ]]; then continue; fi; cat ../scaffAssignment.50000.csv.order.orientation.csv |
egrep -w $c | awk 'BEGIN{OFS="\t"; k=1; kgap=1}; {strand="+"; if($4=="-1"){strand="-"}; print
$3,"AGP","scaffold",k,k+$2-1,0,strand,0,"ID="$1 ; k+=$2; print $3,"AGP","gap",k,k+999,0,"+",0,"ID=gap"kgap++";color=2";
k+=1000; }' > gff_scaffolds/renan.pseudo.v1.$c.gff ; done
## REMOVE the last line of GFFs which is a gap  
sed -i '$d' gff_scaffolds/*.gff
## DIVIDE LENGTHS by 10 (for ACT viewing)  
for gff in gff_scaffolds/*.gff; do awk 'BEGIN{OFS="\t"}; { start=int($4/10); if(start==0){start=1}; print $1, $2, $3,
start, int($5/10), $6, $7, $8, $9}' $gff > $gff.div10.gff; done
## ADD SEQUENCE to the virtual GFFs  
for c in $chr; do file=gff_scaffolds/renan.pseudo.v1.$c.gff.div10.gff; max=$(tail -1 $file | cut -f5); export MAX=$max;
export CHR=$c ; perl -e 'print ">$ENV{CHR}\n", "N" x $ENV{MAX}, "\n"' | fold -120 >> $file; done
## CREATE virtal sequence of CS RefSeq_v1.0  
csfai=/home/share/3bseq/csNrgene/pseudomolecules_v1.0/161010_Chinese_Spring_v1.0_pseudomolecules.fai
for c in $chr; do
cslen=$(cat $csfai | egrep -iw $c | cut -f3);
export CSLEN=$(($cslen/10))
export CHR=$c
perl -e 'print ">$ENV{CHR}\n", "N" x $ENV{CSLEN}, "\n"' | fold -120 > refseq_v1.0.$c.div10.fa
done
```

## STUDY colinearity CS versus renan.pseudo.v1 with bwa ISBPs  
```bash
sbatch --mem=20G --wrap "bwa index renan.pseudo.v1.fa"
dir=/storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/isbp/split/
for file in $dir/*sorted_*.fas ; do bn=$(basename $file); sbatch --cpus-per-task=4 --mem=40G --wrap "bwa mem -t 4
renan.pseudo.v1.fa $file | samtools view -bS - > bwaIsbp/$bn.bam" ; done
```
## RETRIEVE the CS_ISBP coordinates on renan.pseudo.v1 using parseBam.pl  
```bash
for file in $(find bwaIsbp -name '*.bam' |sort); do samtools view -q 1 $file | parseBam.pl -nm 1 | egrep -vw M | tr -s
':-' '\t' | cut -f1-3,5,6 ; done >isbpCS_vs_renan.pseudo.v1.coords.bwa.csv
```
## CREATE ACT files of ISBPs while dividing all lengths by 10 (for ACT viewing)  
```bash
file=isbpCS_vs_renan.pseudo.v1.coords.bwa.csv
for c in $chr ; do
if [[ $c == "chrUn" ]]; then continue; fi;
cat $file | egrep -w $c | awk 'BEGIN{OFS="\t"}; { if($1==$4){ print $1, $4, 100, int(150/10), 0, 0, int($2/10),
int($3/10), int($5/10), int(($5+150)/10), 1e-25, 200} }' >act/$c.isbpCS_vs_renan.pseudo.v1.act;
done
```
## RUN ACT with virtual seq  
```bash
for c in $chr; do echo act renan.pseudo.v1.$c.gff.div10.gff $c.isbpCS_vs_renan.pseudo.v1.act refseq_v1.0.$c.div10.fa ;
done
```

# Manual CURATIONS  