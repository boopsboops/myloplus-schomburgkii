#!/usr/bin/env Rscript

### LOAD LIBS ###
source(here::here("scripts/load-libs.R"))


### LOAD DATA ###
seqs.fas <- read.FASTA(here("assets/sequences-master.fasta"))
seqs.fas.red <- read.FASTA(here("temp/alignments/myloplus-209x621-aligned.fasta"))
# 
master.df <- read_csv(here("assets/tissues-master.csv"),show_col_types=FALSE)
master.df.red <- read_csv(here("temp/alignments/myloplus-209.csv"),show_col_types=FALSE)
master.df.red.probs <- read_csv(here("temp/alignments/myloplus-209-mptp.csv"),show_col_types=FALSE)


### MAKE NJ TREE (ALL HAPS) ###

# mafft align seqs
seqs.fas.ali <- as.matrix(mafft(seqs.fas,method="auto",exec="mafft"))
# trim
seqs.fas.ali <- seqs.fas.ali[,44:dim(seqs.fas.ali)[2]]
# make nj and root
seqs.fas.ali.nj <- phangorn::midpoint(ape::nj(ape::dist.dna(seqs.fas.ali,model="raw",pairwise.deletion=TRUE)))
# remove negative edges
seqs.fas.ali.nj$edge.length[which(seqs.fas.ali.nj$edge.length < 0)] <- 0

# make labels
master.df.lab <- master.df %>% 
    mutate(sciNameValid=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>%
    mutate(labsPhy=paste(label,sciNameValid,waterBody,sep="|")) %>%
    select(label,labsPhy)

# plot
p <- seqs.fas.ali.nj %>% 
    ggtree(color="grey40",ladderize=TRUE,right=TRUE) %<+% master.df.lab +
    geom_tiplab(aes(label=labsPhy),geom="text",offset=0.0005,align=TRUE,color="grey40",size=3) + 
    theme(legend.position="none") +
    xlim(0,0.14)
# plot
filename <- glue("temp/trees/myloplus.tr.nj.",as.character(Sys.Date()),".pdf")
ggsave(filename=here(filename),plot=p,width=297,height=2000,units="mm",limitsize=FALSE)


### BASIC STATS ###

# n samples 
master.df %>% count()
master.df.red %>% count()

# samples per spp (ALL)
master.df %>% 
    mutate(sciNameValid=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>% 
    count(sciNameValid) %>%
    print(n=Inf)

# samples per spp (ALL AVGS)
master.df %>% 
    mutate(sciNameValid=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>% 
    count(sciNameValid) %>%
    summarise(min=min(n),max=max(n),mean=mean(n),median=median(n)) %>%
    print(n=Inf)

# samples per spp (HAPS)
master.df.red %>% 
    count(sciNameValid) %>% 
    print(n=Inf)

# samples per spp avg (HAPS)
master.df.red %>% 
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

# summary table
# n samples 
master.df %>% 
    mutate(sciNameValid=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>% 
    mutate(lonlat=paste(decimalLatitude,decimalLongitude)) %>%
    group_by(sciNameValid) %>%
    summarise(nSamples=length(sciNameValid),nBasins=length(unique(waterBody)),nLocalities=length(unique(lonlat))) %>% 
    ungroup() %>% 
    left_join(rename(count(master.df.red,sciNameValid),nHaps=n)) %>%
    print(n=Inf) #%>%
    #write_csv(here("temp/manuscript/tables/sampling-summary.csv"))

# schomb
master.df %>% 
    mutate(sciNameValid=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>% 
    filter(specificEpithet=="schomburgkii") %>%
    distinct(sciNameValid,waterBody)

# published versus not
master.df %>% 
    mutate(sciNameValid=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>% 
    filter(specificEpithet=="schomburgkii") %>%
    distinct(sciNameValid,label,associatedSequences) %>%
    arrange(associatedSequences) %>%
    print(n=Inf)

# seq lengths (ALL)
dim(seqs.fas.ali)
fas2tab(as.list(seqs.fas.ali)) %>% 
    mutate(nucleotides=str_replace_all(nucleotides,"-","")) %>%
    mutate(length=str_length(nucleotides)) %>%
    summarise(mean=mean(length),med=median(length),min=min(length),max=max(length))

# seq lengths (HAPS)
dim(as.matrix(seqs.fas.red))
fas2tab(seqs.fas.red) %>% 
    mutate(nucleotides=str_replace_all(nucleotides,"-","")) %>%
    mutate(length=str_length(nucleotides)) %>%
    summarise(mean=mean(length),med=median(length),min=min(length),max=max(length))


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

# run for only M. schomburgkii delims

# load clades
clades.df <- read_csv(here("temp/alignments/clades-test.csv"),show_col_types=FALSE)
clades.df <- clades.df %>% select(dbidNex,clade)

# subset dna
seqs.fas.red.schomb <- seqs.fas.red[pull(filter(master.df.red.probs,specificEpithet=="schomburgkii"),dbidNex)]

# get dist matrix and get spider dists
mylo.dist.schomb <- ape::dist.dna(seqs.fas.red.schomb,model="raw",pairwise.deletion=TRUE)


get_group_dists <- function(sp) {
    # subset
    delimit.names <- master.df.red.probs %>% 
        filter(specificEpithet=="schomburgkii") %>% 
        select(dbidNex) %>%
        left_join(filter(clades.df,clade==sp)) %>% 
        replace_na(list(clade="other")) %>%
        rename(delimitName=clade) %>%
        arrange(dbidNex[names(seqs.fas.red.schomb)])
    # make dists per group
    mylo.spp.schomb <- delimit.names %>% pull(delimitName)
    mylo.min.inter.schomb <- nonConDist(distobj=mylo.dist.schomb,sppVector=mylo.spp.schomb,propZero=FALSE,rmNA=FALSE)
    mylo.max.intra.schomb <- maxInDist(distobj=mylo.dist.schomb,sppVector=mylo.spp.schomb,propZero=FALSE,rmNA=FALSE)
    # summarise
    delimit.names.dist <- delimit.names %>% 
        mutate(minInter=mylo.min.inter.schomb,maxIntra=mylo.max.intra.schomb) %>% 
        arrange(delimitName) %>%
        group_by(delimitName) %>% 
        summarise(minInterSpp=min(minInter),maxIntraSpp=max(maxIntra)) %>%
        ungroup()
    return(delimit.names.dist)
}

# test
get_group_dists(sp="mapinguari")

# run on all
group.dists <- bind_rows(mapply(function(x) get_group_dists(sp=x), x=pull(distinct(clades.df,clade),clade), USE.NAMES=FALSE,SIMPLIFY=FALSE))

# filter 
group.dists <- group.dists %>% filter(delimitName!="other") %>% rename(clade=delimitName)

# join with delim results and write out
delim.df <- read_csv(here("temp/alignments/groups-results.csv"),show_col_types=FALSE)
names.df <- read_csv(here("temp/alignments/species-names-ranges.csv"),show_col_types=FALSE)

# tidy up 
delim.df.tidy <- delim.df %>% 
    left_join(group.dists) %>% 
    left_join(names.df) %>% 
    mutate(postDelimProb=round(postDelimProb,digits=2),
        postCladeProb=round(postCladeProb,digits=2),
        minInterSpp=round(minInterSpp,digits=3),
        maxIntraSpp=round(maxIntraSpp,digits=3)) %>%
    arrange(sortOrder) %>% 
    select(species,range,postDelimProb,postCladeProb,minInterSpp,maxIntraSpp,clade)

# write out
print(delim.df.tidy)
#delim.df.tidy %>% write_csv(here("temp/alignments/groups-results-distances.csv"))


### MAKE PARSIMONY TREE (REDUCED HAPS) ###

# nj start tree dropping tips
# seqs.fas.ali.nj$tip.label <- paste("dbid",seqs.fas.ali.nj$tip.label,sep="_")
# seqs.fas.ali.nj.red <- ape::drop.tip(phy=seqs.fas.ali.nj,tip=base::setdiff(seqs.fas.ali.nj$tip.label,names(seqs.fas.red)))
# 
# write out reduced nj tree
# write.nexus(seqs.fas.ali.nj.red,file=here("temp/trees/myloplus.tr.nj.nex"))
# treeannotator -target myloplus.tr.nj.nex -burninTrees 0 -heights ca all.runs.joined.rooted.trees all.runs.joined.rooted.mcc.nj.tre
#
# # get parsinomy score on nj tree
# seqs.fas.red.pd <- as.phyDat(seqs.fas.red)
# phangorn::parsimony(tree=seqs.fas.ali.nj.red,data=seqs.fas.red.pd)
# 
# # run parsimony ratchet
# seqs.fas.ali.pars <- phangorn::pratchet(data=seqs.fas.red.pd,start=seqs.fas.ali.nj.red,all=FALSE,rearrangements="stochastic")
# seqs.fas.ali.pars <- phangorn::acctran(seqs.fas.ali.pars,seqs.fas.red.pd)
# seqs.fas.ali.pars <- phangorn::midpoint(seqs.fas.ali.pars)
# phangorn::parsimony(tree=seqs.fas.ali.pars,data=seqs.fas.red.pd)
# seqs.fas.ali.pars
# 
# # reduce db for plot
# master.df.lab.pars <- master.df.lab %>% mutate(label=paste("dbid",label,sep="_")) %>% filter(label %in% seqs.fas.ali.pars$tip.label)
# 
# # plot
# p <- seqs.fas.ali.pars %>% 
#     ggtree(color="grey40",ladderize=TRUE,right=TRUE) %<+% master.df.lab.pars +
#     geom_tiplab(aes(label=labsPhy),geom="text",offset=0.0005,align=TRUE,color="grey40",size=4) + 
#     theme(legend.position="none") + 
#     xlim(0,150)
# # save
# filename <- glue("temp/trees/myloplus.tr.pars.",as.character(Sys.Date()),".pdf")
# ggsave(filename=here(filename),plot=p,width=297,height=1000,units="mm",limitsize=FALSE)
