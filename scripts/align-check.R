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
library("traits")
library("lubridate")
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/hapCollapse.R")
# renv::install("boopsboops/traits@b62d0431b34d119038e23cec125bb46e8df52e13")

## again for new seqs
tissues.df <- read_csv(here("temp/tissues-master.csv"))
tissues.df %>% print(width=Inf)

#
novos.df <- read_csv(here("temp/Dados_GB_Novos_schom.csv"))
novos.df %>% print(width=Inf)

# setdiff
tissues.df %>% filter(!label %in% pull(novos.df,label))%>% print(n=Inf)
novos.df.new <- novos.df %>% filter(!label %in% pull(tissues.df,label)) 
#novos.df.new %>% write_csv("temp/additions.csv")


# get fresh copies

# FUNCTION TO RUN PARALLEL NCBI_BYID WITH TIMEOUT AND REPEAT
ncbi_byid_parallel <- function(accs){
    start_time <- Sys.time()
    Sys.sleep(time=runif(n=1,min=0,max=3))
    crul::set_opts(http_version=2)
    ncbi.tab <- traits::ncbi_byid(accs,verbose=FALSE)
    if(class(ncbi.tab)!="data.frame") {
        #writeLines("Error found! Repeating ...")
        Sys.sleep(time=3)
        crul::set_opts(http_version=2)
        ncbi.tab <- traits::ncbi_byid(accs,verbose=FALSE)
    } else {
        ncbi.tab <- ncbi.tab
    }
    if(class(ncbi.tab)!="data.frame") {
        stop(writeLines("Searches failed ... aborted")) 
    } else {
        end_time <- Sys.time()
        writeLines(paste0("Metadata for ",length(accs)," accessions downloaded (starting ",accs[1],"). Download took ",round(as.numeric(end_time-start_time),digits=2)," seconds."))
        return(ncbi.tab)
    }
}

# run
novos.df.gb <- ncbi_byid_parallel(pull(novos.df.new,label))

novos.df.gb.clean <- novos.df.gb %>% filter(gi_no!="NCBI_GENOMES") %>% #glimpse()
    distinct(gi_no, .keep_all=TRUE) %>% 
    mutate(lat=paste(str_split_fixed(lat_lon, " ", 4)[,1], str_split_fixed(lat_lon, " ", 4)[,2]), lon=paste(str_split_fixed(lat_lon, " ", 4)[,3], str_split_fixed(lat_lon, " ", 4)[,4])) %>%
    mutate(lat=if_else(grepl(" N",lat), true=str_replace_all(lat," N",""), false=if_else(grepl(" S",lat), true=paste0("-",str_replace_all(lat," S","")), false=lat))) %>%
    mutate(lon=if_else(grepl(" E",lon), true=str_replace_all(lon," E",""), false=if_else(grepl(" W",lon), true=paste0("-",str_replace_all(lon," W","")), false=lon))) %>% 
    mutate(lat=str_replace_all(lat,"^ ", NA_character_), lon=str_replace_all(lon,"^ ", NA_character_)) %>%
    mutate(lat=as.numeric(lat), lon=as.numeric(lon)) %>% 
    select(acc_no,taxon,specimen_voucher,country,lat,lon,collected_by,identified_by,collection_date) %>%
    rename(catalogNumber=specimen_voucher,associatedSequences=acc_no,scientificName=taxon,decimalLatitude=lat,decimalLongitude=lon,verbatimLocality=country,eventDate=collection_date,recordedBy=collected_by,identifiedBy=identified_by) %>%
    mutate(label=str_replace_all(associatedSequences,"\\.[0-9]",""),
        genus=str_split_fixed(scientificName," ",2)[,1],
        specificEpithet=str_split_fixed(scientificName," ",2)[,2],
        country=str_split_fixed(verbatimLocality,":",2)[,1],
        eventDate=str_replace_all(eventDate,"\\..*",". et al."),
        recordedBy=str_replace_all(recordedBy,"\\..*",". et al."),
        eventDate=as.character(lubridate::dmy(eventDate))) %>%
    mutate(otherCatalogNumbers=NA,
        institutionCode=NA,
        collectionCode=NA,
        basisOfRecord=NA,
        typeStatus=NA,
        taxonRank="species",
        kingdom="Animalia",
        phylum="Chordata",
        class="Actinopterygii",
        order="Characiformes",
        family="Serrasalmidae",
        identificationQualifier=NA,
        stateProvince=NA,
        waterBody=NA) %>%
    select(label,associatedSequences,otherCatalogNumbers,catalogNumber,institutionCode,collectionCode,basisOfRecord,typeStatus,scientificName,kingdom,phylum,class,order,family,genus,specificEpithet,taxonRank,identificationQualifier,identifiedBy,decimalLatitude,decimalLongitude,country,stateProvince,waterBody,eventDate,verbatimLocality,recordedBy)

novos.df.gb.clean %>% write_csv("temp/additions-gb.csv")

# reload and check 
tissues.df <- read_csv(here("temp/tissues-master.csv"))
tissues.df %>% print(width=Inf)
unique(pull(tissues.df,scientificName)[-which(paste(str_split_fixed(pull(tissues.df,scientificName)," ",3)[,1], str_split_fixed(pull(tissues.df,scientificName)," ",3)[,2]) == paste(pull(tissues.df,genus),pull(tissues.df,specificEpithet)))])

# check against fasta
old.seqs.fas <- read.FASTA(here("temp/myloplus-unaligned.fasta"))
new.seqs.fas <- read.FASTA(here("temp/myloplus-unaligned-new.fasta"))
renamed.seqs.fas <- read.FASTA(here("temp/sequences_030822.fasta"))

# get diffs
base::setdiff(labels(old.seqs.fas),labels(new.seqs.fas))
base::setdiff(labels(new.seqs.fas),labels(old.seqs.fas))

joined.seqs.fas <- c(old.seqs.fas,new.seqs.fas[base::setdiff(labels(new.seqs.fas),labels(old.seqs.fas))])
base::setdiff(labels(joined.seqs.fas),pull(tissues.df,label))
base::setdiff(pull(tissues.df,label),labels(joined.seqs.fas))

# check joined 
new.joined.seqs.fas <- c(new.seqs.fas,renamed.seqs.fas[base::setdiff(labels(renamed.seqs.fas),labels(new.seqs.fas))])
base::setdiff(labels(joined.seqs.fas),pull(tissues.df,label))
base::setdiff(pull(tissues.df,label),labels(joined.seqs.fas))

# write out
write.FASTA(joined.seqs.fas[sort(labels(joined.seqs.fas))],file="temp/myloplus-unaligned-all.fasta")

#read back in
joined.seqs.fas <- read.FASTA("temp/myloplus-unaligned-all.fasta")

# check with others
base::setdiff(labels(joined.seqs.fas),labels(renamed.seqs.fas))
base::setdiff(labels(renamed.seqs.fas),labels(joined.seqs.fas))

# compare the tables
master.df <- read_csv(here("temp/tissues-master.csv"))
july.df <- read_csv(here("temp/DwC_Val_200722.csv"))
aug.df <- read_csv(here("temp/DwC_Val_200722_v2.csv"))

master.df.spp <- master.df %>% mutate(speciesNameMaster=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>% select(label,speciesNameMaster)
july.df.spp <- july.df %>% mutate(speciesNameJuly=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>% select(label,speciesNameJuly)
aug.df.spp <- aug.df %>% mutate(speciesNameAug=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>% select(label,speciesNameAug)

names.changed <- master.df.spp %>% 
    full_join(july.df.spp) %>% 
    full_join(aug.df.spp) %>% 
    replace_na(list(speciesNameMaster="NA",speciesNameJuly="NA",speciesNameAug="NA")) %>%
    mutate(nameMatch=if_else(speciesNameMaster==speciesNameJuly & speciesNameMaster==speciesNameAug,TRUE,FALSE)) %>%
    mutate(inFasta=if_else(label %in% labels(renamed.seqs.fas) | label %in% labels(new.seqs.fas),TRUE,FALSE)) %>%
    arrange(nameMatch,inFasta,speciesNameMaster,label)# %>% print(n=Inf)
names.changed %>% print(n=Inf)
#names.changed %>% write_csv(here("temp/name-changes.csv"))


# check against fasta
joined.seqs.fas <- read.FASTA("temp/myloplus-unaligned-all.fasta")
master.df <- read_csv(here("temp/tissues-master.csv"))
base::setdiff(labels(joined.seqs.fas),pull(master.df,label))
base::setdiff(pull(master.df,label),labels(joined.seqs.fas))
length(labels(joined.seqs.fas)) 
length(labels(joined.seqs.fas)) == length(pull(master.df,label))

# check missing fasta
latest.fas <- read.FASTA("temp/orig/update_01-11-2022/Myleus_group_new_seqs_311022_red.fasta")
labels(latest.fas)
base::setdiff(labels(latest.fas),labels(joined.seqs.fas))
base::setdiff(labels(joined.seqs.fas),labels(latest.fas))


################################### make tree
joined.seqs.fas <- read.FASTA("temp/myloplus-unaligned-all.fasta")
master.df <- read_csv(here("temp/tissues-master.csv"))
names.changed <- read_csv(here("temp/name-changes.csv"))

# align
myloplus.aligned <- as.matrix(mafft(joined.seqs.fas,exec="mafft"))

# make a tree
myloplus.tr <- nj(dist.dna(myloplus.aligned,model="TN93",pairwise.deletion=TRUE))
myloplus.tr.ml <- phangorn::optim.pml(phangorn::pml(myloplus.tr, phangorn::as.phyDat(myloplus.aligned), k=4, inv=0, model="HKY"), optNni=FALSE, optGamma=TRUE, optInv=FALSE, model="HKY")#rearrangement="stochastic",

# ladderize tree
myloplus.tr <- ladderize(midpoint(myloplus.tr.ml$tree))
myloplus.tr$edge.length[which(myloplus.tr$edge.length < 0)] <- 0

master.df.labs <- master.df %>% 
    mutate(sciName=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>% 
    mutate(labsPhy=paste(label,sciName,waterBody,sep=" | ")) %>%
    mutate(problem=if_else(label %in% pull(names.changed,label),TRUE,FALSE)) %>%
    select(label,labsPhy,problem)


# plot tree
p <- myloplus.tr %>% 
    ggtree(color="grey50") %<+% master.df.labs +
    geom_tiplab(aes(label=labsPhy,color=problem),geom="text",offset=0.0005,align=TRUE) + 
    theme(legend.position="none") +
    scale_color_manual(values=c("grey50","red2")) + 
    xlim(0,0.35)

# plot
ggsave(filename=here("temp/myloplus.tr.nov.pdf"),plot=p,width=297,height=2500,units="mm",limitsize=FALSE)





##################

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

