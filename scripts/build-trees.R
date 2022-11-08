#!/usr/bin/env Rscript

### load libs ###
source(here::here("scripts/load-libs.R"))


#### check names/numbers ###
seqs.fas <- read.FASTA("assets/sequences-master.fasta")
master.df <- read_csv(here("assets/tissues-master.csv"))
base::setdiff(labels(seqs.fas),pull(master.df,label))
base::setdiff(pull(master.df,label),labels(seqs.fas))


### make a quick tree ###

# align
myloplus.aligned <- as.matrix(ips::mafft(seqs.fas,exec="mafft"))

# make a tree
myloplus.tr <- nj(dist.dna(myloplus.aligned,model="TN93",pairwise.deletion=TRUE))
myloplus.tr.ml <- phangorn::optim.pml(phangorn::pml(myloplus.tr, phangorn::as.phyDat(myloplus.aligned), k=4, inv=0, model="HKY"), optNni=FALSE, optGamma=TRUE, optInv=FALSE, model="HKY")#rearrangement="stochastic",

# ladderize tree
myloplus.tr <- ladderize(midpoint(myloplus.tr.ml$tree))
myloplus.tr$edge.length[which(myloplus.tr$edge.length < 0)] <- 0

# 
master.df.labs <- master.df %>% 
    mutate(sciName=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>% 
    mutate(labsPhy=paste(label,sciName,waterBody,sep=" | ")) %>%
    select(label,labsPhy)

# plot tree
p <- myloplus.tr %>% 
    ggtree(color="grey50") %<+% master.df.labs +
    geom_tiplab(aes(label=labsPhy),geom="text",offset=0.0005,align=TRUE,color="grey50") + 
    theme(legend.position="none") +
    xlim(0,0.35)

# plot
filename <- glue("temp/myloplus.tr.",as.character(Sys.Date()),".pdf")
ggsave(filename=here(filename),plot=p,width=297,height=2500,units="mm",limitsize=FALSE)


### dereplicate ###

# convert fas to tab
seqs.fas.df <- fas2tab(seqs.fas) %>% 
    rename(nucleotidesFrag=nucleotides) %>% 
    mutate(lengthFrag=str_length(nucleotidesFrag))

master.df <- read_csv(here("assets/tissues-master.csv"))

# format for derep
master.df %<>% 
    left_join(seqs.fas.df ) %>%
    mutate(sciNameValid=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>%
    mutate(dbid=label,dbidNex=glue("dbid_{label}"))# old way = paste0("dbid_",label)

# collapse
master.df.red <- haps2fas(df=master.df)
glimpse(master.df.red)

# make labels and align
seqs.fas.red <- tab2fas(df=master.df.red,seqcol="nucleotidesFrag",namecol="dbidNex")
seqs.fas.red.ali <- as.matrix(ips::mafft(seqs.fas.red,exec="mafft"))

# write out nexus
write.nexus.data(seqs.fas.red.ali,here("temp-local-only/myloplus-209x664-aligned.nex"),interleaved=FALSE)

# mb commands
#mpirun -np 4 mb myloplus-209x664-aligned-2022-11-07-mb.nex

### load trees and sample ###

# fun to load mrbayes trees
read_t <- function(basename,run,burnin){    
    tree.file.name <- paste0(basename,run,".t")
    mrbayes.trees <- ape::read.nexus(here(tree.file.name))
    mrbayes.trees.burnin <- mrbayes.trees[burnin:length(mrbayes.trees)]
    return(mrbayes.trees.burnin)
}

# test
#read_t(basename="temp-local-only/myloplus-209x664-aligned.nex.run",run=1,burnin=102)

# load four runs
nruns <- 1:4
all.runs <- mapply(function(x) read_t(basename="temp-local-only/myloplus-209x664-aligned.nex.run",run=x,burnin=102),x=nruns,SIMPLIFY=FALSE,USE.NAMES=FALSE)

# join
all.runs.joined <- do.call(c,all.runs)

# sample
set.seed(42)
all.runs.joined.sample <- sample(all.runs.joined,1000)

# name the trees
trees.names <- paste0("tree",str_pad(1:length(all.runs.joined.sample),width=4,pad="0"))

# write out the trees
mapply(function(x,y) ape::write.tree(x,file=here("temp-local-only","mptp",paste0(y,".nwk"))),x=all.runs.joined.sample,y=trees.names)





