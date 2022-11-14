#!/usr/bin/env Rscript

### LOAD LIBS ###
source(here::here("scripts/load-libs.R"))


### LOAD DATA ###
seqs.fas <- read.FASTA(here("assets/sequences-master.fasta"))
seqs.fas.red <- read.FASTA(here("temp/alignments/myloplus-209x621-aligned.fasta"))

master.df <- read_csv(here("assets/tissues-master.csv"),show_col_types=FALSE)
master.df.red <- read_csv(here("temp/alignments/myloplus-209.csv"),show_col_types=FALSE)
master.df.red.probs <- read_csv(here("temp/alignments/myloplus-209-mptp.csv"),show_col_types=FALSE)

### BASIC STATS ###

# n samples 
master.df %>% count()
master.df.red %>% count()

# samples per spp (ALL)
master.df %>% 
    mutate(sciNameValid=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>% 
    count(sciNameValid) %>%
    print(n=Inf)

master.df %>% 
    mutate(sciNameValid=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>% 
    count(sciNameValid) %>%
    summarise(min=min(n),max=max(n),mean=mean(n),median=median(n)) %>%
    print(n=Inf)

# samples per spp (HAPS)
master.df.red %>% 
    mutate(sciNameValid=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>% 
    count(sciNameValid) %>% 
    print(n=Inf)

master.df.red %>% 
    mutate(sciNameValid=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>% 
    count(sciNameValid) %>% 
    summarise(min=min(n),max=max(n),mean=mean(n),median=median(n)) %>%
    print(n=Inf)

# count localities
master.df %>% 
    count(verbatimLocality) %>% 
    print(n=Inf)

# count unique gps
master.df %>% 
    mutate(lonlat=paste(decimalLatitude,decimalLongitude)) %>%
    count(lonlat) %>% 
    print(n=Inf)

# count waterbody
master.df %>% 
    count(waterBody) %>% 
    print(n=Inf)

# distances
delimit.names <- master.df.red.probs %>% 
    mutate(sciNameValid=if_else(is.na(identificationQualifier),paste(genus,specificEpithet,sep="_"),paste(genus,identificationQualifier,sep="_"))) %>% 
    mutate(delimitName=paste(sciNameValid,labelHash,sep="_")) %>% 
    select(dbidNex,delimitName) %>%
    arrange(dbidNex[names(seqs.fas.red)])
#delimit.names %>% print(n=Inf)

# get dist matrix and get spider dists
mylo.dist <- ape::dist.dna(seqs.fas.red,model="raw",pairwise.deletion=TRUE)
mylo.spp <- delimit.names %>% pull(delimitName)
mylo.min.inter <- nonConDist(distobj=mylo.dist,sppVector=mylo.spp,propZero=FALSE,rmNA=FALSE)
mylo.max.intra <- maxInDist(distobj=mylo.dist,sppVector=mylo.spp,propZero=FALSE,rmNA=FALSE)

# join and sort
delimit.names %>% 
    mutate(minInter=mylo.min.inter,maxIntra=mylo.max.intra) %>% 
    arrange(delimitName) %>% 
    print(n=Inf)

# collapse spp
delimit.names %>% 
    mutate(minInter=mylo.min.inter,maxIntra=mylo.max.intra) %>% 
    arrange(delimitName) %>%
    group_by(delimitName) %>% 
    summarise(minInterSpp=min(minInter),maxIntraSpp=max(maxIntra)) %>%
    ungroup() %>%
    print(n=Inf)
