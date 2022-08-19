library(readr)
library(dplyr)
library(ggplot2)


#localisation donnees entree
dirin<-"Y:/ANALYSES/results/3d-dna/"
setwd(dirin)

d=read_tsv(paste(dirin, "Comparaison_scaffolding_3ddna_vs_ISBPpaires.txt", sep=""), col_names=T)


# barplot_nb_scaff_per_put_3ddna_chrom
d%>%
  group_by(putat_3ddna_chr)%>%
  summarise(count=n())%>%
  ggplot(aes(putat_3ddna_chr, count))+geom_col(fill="darkblue", colour="darkblue")+
  ggtitle("Nb of scaff per putative 3ddna chrom (total=2904)")+
  theme(plot.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15))


# barplot_perc_same_put_chrom_assign
d%>%
  group_by(putat_3ddna_chr)%>%
  summarise(count=n(),
            perc_same_put_chrom=sum(Nb_same_put_chrom, na.rm=T)/count*100)%>%
  ggplot(aes(putat_3ddna_chr, perc_same_put_chrom))+geom_col(fill="darkblue", colour="darkblue")+
  ggtitle("Percentage of same putative chrom assignation for 3ddna vs ISBPpaires")+
  theme(plot.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15))

d_mean=data.frame(d%>%
  group_by(putat_3ddna_chr)%>%
  summarise(count=n(),
            perc_same_put_chrom=sum(Nb_same_put_chrom, na.rm=T)/count*100)%>%
  filter(putat_3ddna_chr!="Unchr"))
