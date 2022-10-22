#!/usr/bin/env Rscript

# libs
library("here")
library("ape")
library("spider")
library("tidyverse")
library("magrittr")
library("ips")
library("phangorn")
library("ggtree")
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/hapCollapse.R")

## again for new seqs
tissues.df <- read_csv(here("temp/tissues-master.csv"))
tissues.df %>% print(width=Inf)

#
novos.df <- read_csv(here("temp/Dados_GB_Novos_schom.csv"))
novos.df %>% print(width=Inf)

# setdiff
tissues.df %>% filter(!label %in% pull(novos.df,label))%>% print(n=Inf)
novos.df %>% filter(!label %in% pull(tissues.df,label)) %>% print(n=Inf)


# load table and make a label
tissues.df <- read_csv(here("temp/DwC_Val_200722.csv"))
tissues.df %<>% mutate(label=if_else(is.na(associatedSequences),catalogNumber,associatedSequences))
tissues.df %>% print(width=Inf)

# load fasta and make a table
seqs.fas <- read.FASTA(here("temp/M_schomburgkii_130722.fasta"))
seqs.df <- tibble(id=str_split_fixed(names(seqs.fas),"_",4)[,1],genus=str_split_fixed(names(seqs.fas),"_",4)[,2],species=str_split_fixed(names(seqs.fas),"_",4)[,3],locality=str_split_fixed(names(seqs.fas),"_",4)[,4])

#  labels not in fasta / ids not in labels
tissues.df %>% filter(!label %in% pull(seqs.df,id))
seqs.df %>% filter(!id %in% pull(tissues.df,label)) 
summary(sort(tissues.df$label) == sort(seqs.df$id))


# load up old tissue table from Ota et al. 2020

nigro.df <- read_csv("https://raw.githubusercontent.com/boopsboops/myloplus-spnov/master/data/tissues-master.csv")
pull(nigro.df,associatedSequences)
seqs.df %>% filter(!id %in% pull(nigro.df,associatedSequences)) %>% print(n=Inf)
new.data.df <- tissues.df %>% filter(!label %in% pull(nigro.df,associatedSequences))
new.data.df %>% print(n=Inf)
new.data.df %>% write_csv(here("temp/new-data.csv"))


# read in new data
new.data.rf <- read_csv(here("temp/tissues-master.csv"))

new.data.rf %<>% mutate(label=if_else(is.na(associatedSequences),catalogNumber,associatedSequences))

# check
new.data.rf %>% distinct(scientificName,specificEpithet) %>% arrange(scientificName,specificEpithet) %>% print(n=Inf)

new.data.rf %>% mutate(sciNameShort=paste(str_split_fixed(scientificName," ",3)[,1],str_split_fixed(scientificName," ",3)[,2]), genusspp=paste(genus,specificEpithet)) %>% 
    mutate(tf=if_else(sciNameShort==genusspp,TRUE,FALSE)) %>%
    filter(tf==FALSE) %>%
    distinct(label,scientificName,taxonRank,genus,identificationQualifier,specificEpithet) %>%
    print(n=Inf)

#
new.data.rf.red <- new.data.rf %>% mutate(sciNameOta=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier)),waterBodyNew=waterBody) %>% distinct(label,sciNameOta,waterBodyNew)
old.data.rf.red <- tissues.df %>% mutate(sciNameVal=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>% distinct(label,sciNameVal,waterBody)

seqs.df.names <- seqs.df %>% mutate(sciNameFasta=paste(genus,species)) %>% rename(label=id) %>% select(label,sciNameFasta)


mismatch.names <- full_join(old.data.rf.red,new.data.rf.red) %>% 
   mutate(mismatchName=if_else(sciNameOta==sciNameVal,TRUE,FALSE),mismatchWater=if_else(waterBody==waterBodyNew,TRUE,FALSE)) %>% 
   select(label,sciNameOta,sciNameVal,mismatchName) %>%
   arrange(sciNameOta,label) %>%
   #left_join(seqs.df.names) %>%
   filter(mismatchName==FALSE) %>%
   select(-mismatchName)

mismatch.names %>% print(n=Inf)

mismatch.names %>% write_csv(here("temp/mismatch-names.csv"))


#  labels not in fasta / ids not in labels
new.data.rf %>% filter(!label %in% pull(seqs.df,id))
seqs.df %>% filter(!id %in% pull(new.data.rf,label)) 

###############

# load cleaned fasta
myloplus.unaligned <- read.FASTA(here("temp/myloplus-unaligned.fasta"))

# align
myloplus.aligned <- as.matrix(mafft(myloplus.unaligned,exec="mafft"))

# make a tree
myloplus.tr <- nj(dist.dna(myloplus.aligned,model="raw",pairwise.deletion=TRUE))

# ladderize tree
myloplus.tr <- ladderize(midpoint(myloplus.tr))

myloplus.tr$edge.length[which(myloplus.tr$edge.length < 0)] <- 0

new.data.rf.red.labs <- new.data.rf.red %>% mutate(labs=paste(label,sciNameOta,waterBodyNew,sep=" | "))

# plot tree
p <- myloplus.tr %>% 
    ggtree(color="grey50") %<+% new.data.rf.red.labs +
    geom_tiplab(aes(label=labs),geom="text",offset=0.0005,align=TRUE,color="grey50") + #,align=TRUE
    theme(legend.position="none") +
    #scale_color_manual(values=brewer.pal(9,"Set1")) + 
    xlim(0,0.15)

# plot
ggsave(filename=here("temp/myloplus.tr.pdf"),plot=p,width=297,height=1800,units="mm",limitsize=FALSE)







###########################################
# read data
tissues.df <- read_csv(here("assets/tissues-master.csv"))
# get seqs from GenBank
pull(tissues.df,associatedSequences)[!is.na(pull(tissues.df,associatedSequences))]
pacus.gb <- read.GenBank(access.nb=pull(tissues.df,associatedSequences)[!is.na(pull(tissues.df,associatedSequences))],species.names=FALSE)

# join with new seqs
pacus.new <- read.FASTA("../data/pacus-COI-new.fasta")
pacus.coi <- c(pacus.gb,pacus.new)
seqStat(DNAbin=pacus.coi, thresh=500)


# align
pacus.coi.mat <- as.matrix(mafft(pacus.coi,exec="mafft"))
#write.FASTA(pacus.coi.mat,"../temp/pacus.coi.fasta")

# trim after checking in geneious
pacus.coi.mat.trimmed <- pacus.coi.mat[,70:690]
summary(dim(pacus.coi.mat.trimmed)[2] - checkDNA(DNAbin=pacus.coi.mat.trimmed, gapsAsMissing=TRUE))

# make a tree
pacus.tr <- nj(dist.dna(pacus.coi.mat.trimmed,model="raw",pairwise.deletion=TRUE))

# ladderize tree
pacus.tr <- ladderize(midpoint(pacus.tr))

# make nice labels
tissues.df %<>% mutate(identifier=if_else(!is.na(associatedSequences),associatedSequences,otherCatalogNumbers)) %>%
    mutate(label=if_else(taxonRank!="species" | !is.na(identificationQualifier), paste0(identifier," ",genus," ",identificationQualifier," (",waterBody,")"), paste0(identifier," ",genus," ",specificEpithet," (",waterBody,")"))) 

# how many spp?
tissues.df %>% 
    mutate(label=if_else(taxonRank!="species", paste0(genus," ",identificationQualifier), paste0(genus," ",specificEpithet))) %>% 
    dplyr::select(label) %>% 
    distinct() %>%
    print(n=Inf)

# how many nigrolineatus
tissues.df %>% filter(grepl("nigrolineatus",label))
# how many nigro haps (run haps first)
#tissues.df %>% filter(identifier %in% labels(pacus.coi.haps.mat.trimmed)) %>% filter(grepl("nigrolineatus",label))
# how many localities
tissues.df %>% filter(grepl("nigrolineatus",label)) %>% mutate(loc=paste(decimalLatitude,decimalLongitude)) %>% dplyr::select(loc) %>% distinct()


# colour by specimen
cols <- rep("#BBBBBB",length(pacus.tr$tip.label))
source <- tissues.df$basisOfRecord[match(pacus.tr$tip.label,tissues.df$identifier)]
cols[which(source=="PreservedSpecimen")] <- "#4477AA"
gb <- tissues.df$associatedSequences[match(pacus.tr$tip.label,tissues.df$identifier)]
cols[which(is.na(gb))] <- "#ff5f19"

# match
pacus.tr$tip.label <- tissues.df$label[match(pacus.tr$tip.label,tissues.df$identifier)]

# remove neg branches
pacus.tr$edge.length[which(pacus.tr$edge.length<0)] <- 0

# write out the tree
pdf(file="../temp/nj.tree.406.cols.pdf",width=30,height=80)
plot.phylo(pacus.tr,no.margin=TRUE,font=1,label.offset=0.0003,tip.color=cols,edge.color="#CCBB44",edge.width=3)
dev.off()

# write out a nex
write.nexus.data(pacus.coi.mat.trimmed,file="../temp-local-only/pacus.trimmed.nex",format="dna",interleaved=FALSE)


# to same for haplotypes only
pacus.coi.haps <- hapCollapse(data=pacus.coi,cores=8)
pacus.coi.haps.mat <- as.matrix(mafft(pacus.coi.haps,exec="mafft"))
pacus.coi.haps.mat.trimmed <- pacus.coi.haps.mat[,70:690]
write.nexus.data(pacus.coi.haps.mat.trimmed,file="../temp-local-only/pacus.haps.trimmed.nex",format="dna",interleaved=FALSE)

# sample trees from beast 
trees.combined <- read.nexus(file="../temp-local-only/pacus.haps.trimmed.combined.trees")

# remove burnin
trees.combined <- trees.combined[337:3336]

# trees subsampled
set.seed(42)
trees.subsampled <- sample(trees.combined,1000)

# write out 
write.nexus(trees.subsampled, file="../data/pacus.COI.trees")

# quickly rename the new seqs with GenBank accs
dat <- read.nexus.data(file="../temp-local-only/pacus.trimmed.nex")
new <- tissues.df %>% filter(associatedSequences %in% setdiff(pull(tissues.df,associatedSequences),names(dat)))
names(dat)[which(names(dat) %in% new$otherCatalogNumbers)] <- new$associatedSequences[match(names(dat),new$otherCatalogNumbers)][!is.na(new$associatedSequences[match(names(dat),new$otherCatalogNumbers)])]
write.nexus.data(dat,file="../temp-local-only/pacus.trimmed.bg.nex",format="dna",interleaved=FALSE)


# old checl code
?setdiff

sup.df <- read_csv("../../serrasalmids/data/supplementary_table1.csv")

# species not in 
sup.df.my <- sup.df %>% filter(genus=="Tometes" | genus=="Myloplus" | genus=="Mylesinus" | genus=="Myleus" | genus=="Ossubtus" | genus=="Utiaritichthys" | genus=="Acnodon")
tissues.df.my <- tissues.df %>% filter(genus=="Tometes" | genus=="Myloplus" | genus=="Mylesinus" | genus=="Myleus" | genus=="Ossubtus" | genus=="Utiaritichthys" | genus=="Acnodon")

# read fas
pacus.coi <- read.FASTA("../data/pacus-COI.fasta")

# make a label col
#tissues.df %<>% mutate(label=if_else(taxonRank=="species",paste(catalogNumber,genus,specificEpithet,waterBody),paste(catalogNumber,genus,identificationQualifier,waterBody)))

# check names
setdiff(pull(sup.df.my,otherCatalogNumbers),pull(tissues.df.my,catalogNumber))
setdiff(pull(tissues.df.my,catalogNumber),pull(sup.df.my,otherCatalogNumbers))

setdiff(names(pacus.coi), pull(tissues.df.my,catalogNumber))
setdiff(pull(tissues.df.my,catalogNumber),names(pacus.coi))

setdiff(names(pacus.coi), pull(sup.df.my,otherCatalogNumbers))# in fasta, not in machado2019
setdiff(pull(sup.df.my,otherCatalogNumbers),names(pacus.coi))# in machado2019, not in fasta

# subset the new ones from the fasta and df
pacus.coi.new <-  pacus.coi[setdiff(names(pacus.coi), pull(sup.df.my,otherCatalogNumbers))]
tissues.df.my.new <- tissues.df.my %>% filter(catalogNumber %in% setdiff(names(pacus.coi), pull(sup.df.my,otherCatalogNumbers)))
# write out
write.FASTA(pacus.coi.new,file="../data/pacus-COI-new.fasta")
write_csv(tissues.df.my.new,path="../data/tissues-new.csv")
write_csv(sup.df.my,path="../data/tissues-master.csv")

