library(readr)
library(tidyr)
library(stringr)
library(dplyr)
library(ggplot2)
library(tibble)
library(scales)
#library(HiTC)
#browseVignettes("HiTC")


#localisation donnees entree
dirin<-"Y:/ANALYSES/results/juicer/"
setwd(dirin)

################################## SCAFF SIZES ###################################
dsize=read_tsv(paste(dirin, "Triticum_aestivum_RENAN_v2.fasta.fai", sep=""), col_names=F)
dsize=dsize[,1:2]
colnames(dsize)=c("scaff1", "size1")
dsize=dsize%>%mutate(scaff1=gsub("Taestivum_RENAN_scaffold_", "scaff_", scaff1),
                     scaff2=scaff1,
                     size2=size1)

###### MAKE A BED
# colnames(dsize)=c("scaff", "start")
# dsize=dsize%>%mutate(scaff=gsub("Taestivum_RENAN_scaffold_", "scaff_", scaff),
#                      end=cumsum(start))
# vec_start=dsize$end+1
# vec_start[2904]=1
# dsize=dsize%>%mutate(start=sort(vec_start),
#                      bin=scaff)
# 
# write_tsv(dsize, paste(dirin,"xgi.bed", sep=""), quote_escape=FALSE, col_names=F)

################################## scaff pseudo assignment ###################################
pseudo=read_tsv(paste(dirin, "RENAN_V2_scaffAssignment.TAB", sep=""), col_names=T)

pseudo=pseudo%>%mutate(scaff1=gsub("Taestivum_RENAN_scaffold_", "scaff_", scaff),
                       scaff2=scaff1,
                       rank1=row_number(),
                       rank2=rank1,
                       chr2=chr)

################################## scaffolding 3DDNA ###################################
ddna=read_tsv(paste(dirin, "RENAN_V2_3DDNA_scaff_order.txt", sep=""), col_names=T)
colnames(ddna)=c("scaffolding_3ddna","scaff_order_3ddna","scaff1","size","strand_3ddna")

ddna=ddna%>%mutate(scaff1=paste("scaff_",scaff1, sep=""),
                   rank1=row_number(),
                   scaff2=scaff1,
                   rank2=rank1)

########################### Contacts list #####################################

d=read_tsv(paste(dirin, "hic_RENAN_V2_mapQ30_scaff_paires.txt", sep=""), col_names=F)
colnames(d)=c("scaff1", "scaff2", "contacts")

d=d%>%mutate(scaff1=as.numeric(gsub("Taestivum_RENAN_scaffold_", "", scaff1)),
             scaff2=as.numeric(gsub("Taestivum_RENAN_scaffold_", "", scaff2)))
d=d%>%arrange(scaff1, scaff2)

#pour rendre la matrice symetrique  --> nul besoin finalement
d2=data.frame(scaff1=setdiff(unique(d$scaff1), unique(d$scaff2)),
              scaff2=setdiff(unique(d$scaff1), unique(d$scaff2)),
              contacts=rep(0, times=length(setdiff(unique(d$scaff1), unique(d$scaff2))))
              )

d3=data.frame(scaff1=setdiff(unique(d$scaff2), unique(d$scaff1)),
              scaff2=setdiff(unique(d$scaff2), unique(d$scaff1)),
              contacts=rep(0, times=length(setdiff(unique(d$scaff2), unique(d$scaff1))))
              )
d=rbind(d, d2, d3)

d=d%>%arrange(scaff1, desc(contacts))
d=d%>%mutate(scaff2=paste("scaff_",scaff2, sep=""),
             scaff1=paste("scaff_",scaff1, sep=""))

################################################################################
#################### HiC MAP VC NORMALIZATION & ORDER PSEUDO ISBP ##############

################## Contacts MATRICE ############################################
mat=d[,c("scaff1","scaff2","contacts")]%>%
  pivot_wider(names_from = scaff2, values_from = contacts)

mat=column_to_rownames(mat, var="scaff1")

mat=as.matrix(mat)

######### VANILLA COVERAGE NORMALIZATION ###########################

# facteur de correction des lignes (normalization 1/square root() en plus ?): 
L=1/sqrt(apply(data.frame(mat),1,sum, na.rm=T))
# facteur de correction des colonnes (normalization 1/square root() en plus ?): 
C=1/sqrt(apply(data.frame(mat),2,sum, na.rm=T))

# test=matrix(c(4,2,8,9,4,1,2,3,6), nrow=3, ncol=3)
# vec=c(3,1,10)
# test*vec

# matrice ajustee pour les lignes
mat=mat*L

# transposition puis matrice ajustee pour les colonnes
mat=t(mat)
mat=mat*C
mat=t(mat)

mat_to_sort=mat

scaff1=row.names(mat)
mat_tibb=as_tibble(mat)
mat_tibb=cbind(scaff1, mat_tibb)

dplot_VCN=mat_tibb%>%
  pivot_longer(cols=starts_with("scaff_") ,names_to="scaff2", values_to ="contacts",
               names_repair="minimal", values_drop_na=TRUE)

##### HiC MAP: ORDER SCAFF ACCORDING PSEUDO RENAN_V2 (ISBP) ####################

dplot_VCN_pseudo=right_join(pseudo[,c("scaff2","chr2","rank2")], dplot_VCN, by="scaff2")
dplot_VCN_pseudo=right_join(pseudo[,c("scaff1","chr","rank1")], dplot_VCN_pseudo, by="scaff1")


pdf(paste(dirin, "hic_maps/HiCmaps_RENAN_V2_pseudo_order_mapQ30.pdf", sep=""), width = 12, height = 8)

#HiC MAP
g1=ggplot(dplot_VCN_pseudo, aes(rank1, rank2, fill=contacts))+
  geom_tile(show.legend=T, na.rm=T)+
  scale_fill_gradient(limits=c(0,mean(dplot_VCN$contacts, na.rm=T)+sd(dplot_VCN$contacts, na.rm=T)),
                      low="white", high="red")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), trans="reverse")+
  theme_minimal()+
  labs(x="scaff ordered / pseudo", y="scaff ordered / pseudo")+
  ggtitle("HiCmap _ pseudomol scaff order _ mapQ30", sub="VC normalization")

plot(g1)

for (i in (1:(length(unique(dplot_VCN_pseudo$chr))-1))){
  
  chr_c=unique(dplot_VCN_pseudo$chr)[i]
  dplot_VCN_chr=dplot_VCN_pseudo[dplot_VCN_pseudo$chr==chr_c & dplot_VCN_pseudo$chr2==chr_c,]
  limitmax=mean(dplot_VCN_chr$contacts, na.rm=T)+sd(dplot_VCN_chr$contacts, na.rm=T)

  g=ggplot(dplot_VCN_chr, aes(rank1, rank2, fill=contacts))+
    geom_tile(show.legend=F, na.rm=T)+
    scale_fill_gradient(limits=c(0,limitmax), low="white", high="red")+
    scale_x_continuous(expand=c(0,0), breaks=dplot_VCN_chr$rank1 ,labels=dplot_VCN_chr$scaff1)+
    scale_y_continuous(expand=c(0,0), breaks=dplot_VCN_chr$rank2, labels=dplot_VCN_chr$scaff2, trans="reverse")+
    theme_classic()+
    theme(axis.text.x=element_text(size=6, angle=90, vjust=0 ))+
    theme(axis.text.y=element_text(size=6, hjust=0))+
    labs(x=element_blank(), y=element_blank())+
    ggtitle(paste("HiCmap",chr_c,"_ pseudomol scaff order _ mapQ30"), sub="VC normalization")
  plot(g)
}
dev.off()

################################################################################
############ HiC MAP SCAFF_SIZE NORMALIZATION & ORDER PSEUDO ISBP ##############

######### SHORTER SCAFF SIZE NORMALIZATION #####################################

dplot=d%>%
  left_join(dsize[,c("scaff1","size1")], by="scaff1")%>%
  left_join(dsize[,c("scaff2","size2")], by="scaff2")%>%
  mutate(contacts_norm=if_else(size1<=size2, contacts/size1, contacts/size2))

##### HiC MAP: ORDER SCAFF ACCORDING PSEUDO RENAN_V2 (ISBP) ####################

dplot_pseudo=right_join(pseudo[,c("scaff2","chr2","rank2")], dplot, by="scaff2")
dplot_pseudo=right_join(pseudo[,c("scaff1","chr","rank1")], dplot_pseudo, by="scaff1")


pdf(paste(dirin, "hic_maps/HiCmaps_RENAN_V2_pseudo_order_size_norm_mapQ30.pdf", sep=""), width = 12, height = 8)

# HiC MAP
g1=ggplot(dplot_pseudo, aes(rank1, rank2))+
  geom_tile(aes(fill=contacts_norm), show.legend=T, na.rm=T, stat="identity")+
  scale_fill_gradient(limits=c(0,mean(dplot_pseudo$contacts_norm, na.rm=T)+sd(dplot_pseudo$contacts_norm, na.rm=T)),
                          low="white", high="red")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), trans="reverse")+
  theme_minimal()+
  labs(x="scaff ordered / pseudo", y="scaff ordered / pseudo")+
  ggtitle("HiC map _ pseudomol scaff order _ mapQ30", sub="smaller scaff size normalization")
plot(g1)

for (i in (1:(length(unique(dplot_pseudo$chr))-1))){
  
  chr_c=unique(dplot_pseudo$chr)[i]
  dplot_chr=dplot_pseudo[dplot_pseudo$chr==chr_c & dplot_pseudo$chr2==chr_c,]
  limitmax=mean(dplot_chr$contacts_norm, na.rm=T)+sd(dplot_chr$contacts_norm, na.rm=T)
  
  g=ggplot(dplot_chr, aes(rank1, rank2))+
    geom_tile(aes(fill=contacts_norm), show.legend=F, na.rm=T, stat="identity")+
    scale_fill_gradient(limits=c(0,limitmax), low="white", high="red")+
    scale_x_continuous(expand=c(0,0), breaks=dplot_chr$rank1 ,labels=dplot_chr$scaff1)+
    scale_y_continuous(expand=c(0,0), breaks=dplot_chr$rank2, labels=dplot_chr$scaff2, trans="reverse")+
    theme_classic()+
    theme(axis.text.x=element_text(size=6, angle=90, vjust=0 ))+
    theme(axis.text.y=element_text(size=6, hjust=0))+
    labs(x=element_blank(), y=element_blank())+
    ggtitle(paste("HiCmap",chr_c,"_ pseudomol scaff order _ mapQ30"), sub="smaller scaff size normalization")
  plot(g)
  
}
dev.off()

################################################################################
############ HiC MAP VC NORMALIZATION & 3DDNA ORDER ############################

######### VANILLA COVERAGE NORMALIZATION #####################################
# it is dplot_VCN, already done

######### HiC MAP: ORDER SCAFF ACCORDING 3DDNA #################################

dplot_VCN_3ddna=right_join(ddna[ddna$scaffolding_3ddna=="megascaff1",c("scaff2","rank2")], dplot_VCN, by="scaff2")
dplot_VCN_3ddna=right_join(ddna[ddna$scaffolding_3ddna=="megascaff1",c("scaff1","rank1")], dplot_VCN_3ddna, by="scaff1")

pdf(paste(dirin, "hic_maps/HiCmaps_RENAN_V2_3ddna_order_mapQ30.pdf", sep=""), width = 12, height = 8)

# HiC MAP
g2=ggplot(dplot_VCN_3ddna, aes(rank1, rank2))+
  geom_tile(aes(fill=contacts), show.legend=T, na.rm=T, stat="identity")+
  scale_fill_gradient(limits=c(0,mean(dplot_VCN_3ddna$contacts, na.rm=T)+sd(dplot_VCN_3ddna$contacts, na.rm=T)),
                      low="white", high="red")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), trans="reverse")+
  theme_minimal()+
  labs(x="scaff ordered / 3DDNA", y="scaff ordered / 3DDNA")+
  ggtitle("HiC map _ 3DDNA scaff order _ mapQ30", sub="VC normalization")
plot(g2)


################################################################################
############ HiC MAP SCAFF_SIZE NORMALIZATION & 3DDNA ORDER ####################

######### SHORTER SCAFF SIZE NORMALIZATION #####################################
# it is dplot, already done

######### HiC MAP: ORDER SCAFF ACCORDING 3DDNA #################################

dplot_3ddna=right_join(ddna[ddna$scaffolding_3ddna=="megascaff1",c("scaff2","rank2")], dplot, by="scaff2")
dplot_3ddna=right_join(ddna[ddna$scaffolding_3ddna=="megascaff1",c("scaff1","rank1")], dplot_3ddna, by="scaff1")

# HiC MAP
g3=ggplot(dplot_3ddna, aes(rank1, rank2))+
  geom_tile(aes(fill=contacts_norm), show.legend=T, na.rm=T, stat="identity")+
  scale_fill_gradient(limits=c(0,mean(dplot_3ddna$contacts_norm, na.rm=T)+sd(dplot_3ddna$contacts_norm, na.rm=T)),
                      low="white", high="red")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), trans="reverse")+
  theme_minimal()+
  labs(x="scaff ordered / 3DDNA", y="scaff ordered / 3DDNA")+
  ggtitle("HiC map _ 3DDNA scaff order _ mapQ30", sub="smaller scaff size normalization")
plot(g3)

dev.off()


##############################################
##############################################

# remplacer les NA par des 0
for (i in 1:dim(mat)[2]){
  mat[is.na(mat[,i]),i]=0
}

#test tri matrice
require(clue)
pMatrix <- function(A, B) { 
  #finds the permutation P of A such that ||PA - B|| is minimum in Frobenius norm 
  # Uses the linear-sum assignment problem (LSAP) solver in the "clue" package 
  
  # Returns P%*%A and the permutation vector `pvec' such that 
  # A[pvec, ] is the permutation of A closest to B 
  n <- nrow(A) 
  D <- matrix(NA, n, n) 
  for (i in 1:n) {
    print(i)
    for (j in 1:n) { 
      D[j, i] <- (sum((B[j, ] - A[i, ])^2)) 
    }
  } 
  vec <- c(solve_LSAP(D, maximum=T)) 
  list(A=A[vec,], pvec=vec) 
} 

B=diag(1, nrow(mat_to_sort))

mat_sorted=pMatrix(mat_to_sort, B)


# #tri lignes
# for (i in 1:ncol(mat_to_sort) ){
#   mat_to_sort=mat_to_sort[order(mat_to_sort[,i]),]
# }
# 
# #transposition
# mat_to_sort=t(mat_to_sort)
# 
# t_mat_lignes=mat_to_sort
# mat_to_sort=t_mat_lignes
# 
# #tri colonnes
# for (i in sort(1:ncol(mat_to_sort), decreasing=T) ){
#   mat_to_sort=mat_to_sort[order(mat_to_sort[,i]),]
# }

scaff1=as_tibble(row.names(mat_to_sort))
scaff2=as_tibble(colnames(mat_to_sort))
colnames(scaff1)="scaff1"
colnames(scaff2)="scaff2"
mat_sorted=as_tibble(mat_to_sort)
mat_sorted=cbind(scaff1, mat_sorted)

dplot_VCN_mat_sort=mat_sorted%>%
  pivot_longer(cols=starts_with("scaff_") ,names_to="scaff2", values_to ="contacts",
               names_repair="minimal", values_drop_na=TRUE)

#JOIN RANK INFO
scaff1=scaff1%>%mutate(rank1=row_number())
scaff2=scaff2%>%mutate(rank2=row_number())

dplot_VCN_mat_sort=right_join(scaff2, dplot_VCN_mat_sort, by="scaff2")
dplot_VCN_mat_sort=right_join(scaff1, dplot_VCN_mat_sort, by="scaff1")

#HiC MAP
ggplot(dplot_VCN_mat_sort, aes(rank1, rank2, fill=contacts))+
  geom_tile(show.legend=F, na.rm=T)+
  scale_x_continuous(expand=c(0.01,0.01))+
  scale_y_continuous(expand=c(0.01,0.01))+
  ggtitle("HiC map CONTACTS SORTED _ mapQ30", sub="VC normalization")

#######################
#######################



############################## test HiTC package ################################

# #dtsv=d%>%mutate(contacts=paste("x", contacts, sep=""))
# write_tsv(d, paste(dirin,"hic_contacts_RENAN_V2.tsv", sep=""), quote_escape=FALSE, col_names = F)
# 
# importC(con="Y:/ANALYSES/results/juicer/hic_contacts_RENAN_V2.tsv",
#         xgi="Y:/ANALYSES/results/juicer/xgi.bed", ygi=NULL, allPairwise=T)

#https://rdrr.io/github/mckf111/hiceize/src/R/HiC_map.R#sym-triViewC
#https://cran.r-project.org/web/packages/adjclust/vignettes/hicClust.html#session-information

