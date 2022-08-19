# hic_assembly

hic assembly of Renan genome, WHEATOMICS project, WP1
from /home/palasser/wheatomics_wp1/hic_assembly_TaeRenan_v2.0 on HPC2 cluster

## Git repository initialization

```
cd /home/palasser/wheatomics_wp1/hic_assembly
git remote add hic https://forgemia.inra.fr/gdec-bioinfo/wheatomics/hic_assembly.git
git status
git add bin/*
git commit -m "initial commit"
git push -u hic --all
git branch -M main
git branch -a
git fetch hic main
git pull hic main
```

## Description  
The goal of this part of the project is to assemble the genome of the bread wheat French culivar Renan.
Input data are Oxford Nanopore Technology sequences assembled into 2904 scaffolds:

Input data are on TGCC cluster: /ccc/genostore/cont007/fg0118/fg0118/assemblies
Triticum_aestivum_RENAN_v2.fasta
Triticum_aestivum_RENAN_v2.stats

#-------------------- GLOBAL STATISTICS -------------------#
N50 size= 48460660  number= 79
N80 size= 14449349  number= 241
N90 size= 7948299  number= 376
Cumulative size= 14259284930 number= 2904 minSize= 30021 maxSize= 253607342 averageSize= 4910222 auN= 67017177
#----------------------------------------------------------#
#-------------------- SIZE REPARTITION --------------------#
Size= >= 10000000       Number= 315        (10.85)      CumulativeSize= 12293093327     (86.21)
Size= >= 5000000        Number= 467        (16.08)      CumulativeSize= 13418180023     (94.10)
Size= >= 1000000        Number= 706        (24.31)      CumulativeSize= 14032844541     (98.41)
Size= >= 100000         Number= 1283       (44.18)      CumulativeSize= 14154985870     (99.27)
Size= >= 50000          Number= 2451       (84.40)      CumulativeSize= 14241435579     (99.87)
Size= >= 10000          Number= 2904       (100.00)     CumulativeSize= 14259284930     (100.00)
#----------------------------------------------------------#
#-------------------- BASE COMPOSITION --------------------#
NumberOfN= 258162674 (1.81%) NumberOfGC= 6480316918 (46.28%)
#----------------------------------------------------------#

These scaffolds have been oredered and assembled into 21 pseudomolecules based on ISBPs colinearity with iwgsc_refseqv1.
HiC data sequence (35X, Arima) have then been mapped on such ordered scaffolds to produce a HiC map with juicer and 3d-dna.
The draft assembly has been improved using Juicebox Assembly Tools software.
Final TaeRenan_v2.0 reference genome has been finaly built with scripts developped at GDEC research unit (nCount.pl, closestGAP.py, subseq.pl, buildPseudomol.pl).

## Issue
It was far enough to use 40Gb of raw HiC sequences as input in order to obtain reliable results. 
More than 40Gb were not necessary and make the analysis technically complicated and very long.

## HiC map  
![HiCmap](/Renan_v13_v2.pseudo.v2.0.svg)

## Knowledge transmission
See [GDEC's Wiki](https://wiki.inra.fr/wiki/umr1095/Project+Bioinfo/howto-hic)

## Support  
frederic.choulet@inrae.fr, pauline.lasserre-zuber@inrae.fr, helene.rimbert@inrae.fr, jmaury@genoscope.cns.fr

## Roadmap  
Assemble unanchored yet scaffolds
Introgressions studies
UTR and alternative transcripts annotation: in progress, contact helene.rimbert@inrae.fr

## Authors and acknowledgment  
Jean-Marc AURY (Genoscope), Frederic CHOULET (INRAe), Cecile MONAT, Helene RIMBERT (INRAe), Philippe LEROY (INRAe), Nathan PAPON (INRAe)

## Project status  
Finished.
Paper on line: [Long-read and chromosome-scale assembly of the hexaploid wheat genome achieves high resolution for research and breeding](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac034/6575388)
