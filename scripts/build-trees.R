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
#write_csv(master.df.red,here("temp-local-only/myloplus-209x664-aligned.csv"))

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

# root and ladderise trees
all.runs.joined.sample.rooted <- mapply(function(x) ape::ladderize(phangorn::midpoint(x)), x=all.runs.joined.sample,SIMPLIFY=FALSE,USE.NAMES=FALSE)

#plot(all.runs.joined.sample.rooted[[1]])

# name the trees
trees.names <- paste0("tree",str_pad(1:length(all.runs.joined.sample.rooted),width=4,pad="0"))

# write out the trees
mapply(function(x,y) ape::write.tree(x,file=here("temp-local-only","mptp",paste0(y,".nwk"))),x=all.runs.joined.sample.rooted,y=trees.names)

# also write out a multiphylo to view
class(all.runs.joined.sample.rooted) <- "phylo"
ape::write.tree(all.runs.joined.sample.rooted,file=here("temp-local-only/all.runs.joined.sample.rooted.nwk"))
ape::write.nexus(all.runs.joined.sample.rooted,file=here("temp-local-only/all.runs.joined.sample.rooted.trees"))


# to make mcc tree
#treeannotator -burninTrees 0 -heights ca temp-local-only/all.runs.joined.sample.rooted.trees temp-local-only/all.runs.joined.sample.rooted.mcc.tre

# load up reduced annotated df and mcc tree
master.df.red <- read_csv(here("temp-local-only/myloplus-209x664-aligned.csv"))
mcc.tre <- read.nexus("temp-local-only/all.runs.joined.sample.rooted.mcc.tre")
master.df.plot <- master.df.red %>% mutate(labsPhy=paste(label,sciNameValid,waterBody,nHaps,sep=" | ")) %>%
    select(dbidNex,labsPhy)
# write out the mcc
mcc.tre$tip.label <- pull(master.df.plot,labsPhy)[match(mcc.tre$tip.label,pull(master.df.plot,dbidNex))]
#plot(mcc.tre)
#ape::write.tree(mcc.tre,file=here("temp-local-only/mcc.nwk"))


# plot tree
p <- mcc.tre %>% 
    ggtree(color="grey50") %<+% master.df.plot +
    geom_tiplab(aes(label=labsPhy),geom="text",offset=0.0005,align=TRUE,color="grey50") + 
    theme(legend.position="none") +
    xlim(0,0.25)

# plot
filename <- glue("temp/myloplus.tr.",as.character(Sys.Date()),".pdf")
ggsave(filename=here(filename),plot=p,width=297,height=1000,units="mm",limitsize=FALSE)

# mptp 
run_mptp <- function(file,threshold,minbr) {
    string.mptp <- paste0("mptp --ml --",threshold," --minbr ",minbr," --tree_file ",file," --output_file ",file,".mptp.out")
    system(command=string.mptp,ignore.stdout=FALSE)
    #return(string.mptp)
}

run_mptp(file=here("temp-local-only/mcc.nwk"),threshold="multi",minbr=0.0001)


### RAXML NG ###

# NEW RAXML-NG FUN
raxml_ng <- function(file,verbose) {
    if(verbose == "true") {
        string.mafft <- paste0("mafft --thread -1 --maxiterate 2 --retree 2 ",file," > ",file,".ali")
        system(command=string.mafft,ignore.stdout=FALSE)
        string.parse <- paste0("raxml-ng --parse --msa ",file,".ali --model TN93+F+G --seed 42 --redo --threads auto")
        system(command=string.parse,ignore.stdout=FALSE)
        string.search <- paste0("raxml-ng --search --msa ",file,".ali.raxml.rba --tree pars{1} --lh-epsilon 0.001 --seed 42 --redo --threads auto")
        system(command=string.search,ignore.stdout=FALSE)
        rax.tr <- ape::read.tree(file=paste0(file,".ali.raxml.rba.raxml.bestTree"))
    } else if (verbose == "false") {
        string.mafft <- paste0("mafft --quiet --thread -1 --maxiterate 2 --retree 2 ",file," > ",file,".ali")
        system(command=string.mafft,ignore.stdout=FALSE)
        string.parse <- paste0("raxml-ng --parse --msa ",file,".ali --model TN93+F+G --seed 42 --redo --threads auto")
        system(command=string.parse,ignore.stdout=TRUE)
        string.search <- paste0("raxml-ng --search --msa ",file,".ali.raxml.rba --tree pars{1} --lh-epsilon 0.001 --seed 42 --redo --threads auto")
        system(command=string.search,ignore.stdout=TRUE)
        rax.tr <- ape::read.tree(file=paste0(file,".ali.raxml.rba.raxml.bestTree"))
    } else stop(writeLines("'-v' value must be 'true' or 'false'."))
    return(rax.tr)
}

#
# write alignment
write.FASTA(tab2fas(df=master.df.red,seqcol="nucleotidesFrag",namecol="dbidNex"),here("temp-local-only/myloplus-209x664.fas"))

rax.tr <- raxml_ng(file=here("temp-local-only/myloplus-209x664.fas"),verbose="true")
rax.tr.root <- ladderize(midpoint(rax.tr))
rax.tr.root$tip.label <- pull(master.df.plot,labsPhy)[match(rax.tr.root$tip.label,pull(master.df.plot,dbidNex))]
ape::write.tree(rax.tr.root,file=here("temp-local-only/ml.nwk"))

run_mptp(file=here("temp-local-only/ml.nwk"),threshold="single",minbr=0.0001)
run_mptp(file=here("temp-local-only/ml.nwk"),threshold="single",minbr=0)
run_mptp(file=here("temp-local-only/ml.nwk"),threshold="multi",minbr=0.0001)
run_mptp(file=here("temp-local-only/ml.nwk"),threshold="multi",minbr=0)

# FUNCTION TO READ MPTP OUTPUT FILES
read_mptp <- function(file) {
    mptp.scan <- scan(file=file,what="character",sep="\n",quiet=TRUE)
    skiplines <- grep("Species 1:",mptp.scan)
    skiplines <- skiplines - 1
    writeLines(mptp.scan[1:skiplines])
    mptp.raw <- readr::read_delim(file,skip=skiplines,delim=",",col_names="label",show_col_types=FALSE)
    mptp.tab <- mptp.raw %>% 
        dplyr::mutate(mptpDelim=ifelse(grepl(":",label),label,NA)) %>%
        tidyr::fill(mptpDelim,.direction="down") %>%
        dplyr::filter(!grepl(":",label)) %>%
        dplyr::mutate(mptpDelim=stringr::str_replace_all(mptpDelim,":","")) %>%
        dplyr::mutate(mptpDelim=stringr::str_replace_all(mptpDelim,"Species ","")) %>%
        dplyr::mutate(mptpDelim=paste0("mptp",str_pad(mptpDelim,pad="0",width=4))) %>%
        dplyr::relocate(mptpDelim,.before=label)
    return(mptp.tab)
}

mptp.out <- read_mptp(file=here("temp-local-only/ml.nwk.mptp.out.txt"))

mptp.out %>% arrange(mptpDelim,label) %>% group_by(mptpDelim) %>% mutate(labelCat=paste(label,collapse=","),labelHash=openssl::md5(labelCat))

