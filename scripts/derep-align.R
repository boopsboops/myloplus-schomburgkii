#!/usr/bin/env Rscript

### LOAD LIBS ###
source(here::here("scripts/load-libs.R"))


### DEREPLICATE SEQUENCES ###

# load data
seqs.fas <- read.FASTA("assets/sequences-master.fasta")
master.df <- read_csv(here("assets/tissues-master.csv"),show_col_types=FALSE)

# convert fas to tab
seqs.fas.df <- fas2tab(seqs.fas) %>% 
    rename(nucleotidesFrag=nucleotides) %>% 
    mutate(lengthFrag=str_length(nucleotidesFrag))

# format for derep
master.df %<>% 
    left_join(seqs.fas.df ) %>%
    mutate(sciNameValid=if_else(is.na(identificationQualifier),paste(genus,specificEpithet),paste(genus,identificationQualifier))) %>%
    mutate(dbid=label,dbidNex=glue("dbid_{label}"))# old way = paste0("dbid_",label)

# collapse
master.df.red <- haps2fas(df=master.df)
glimpse(master.df.red)
# write out
#write_csv(master.df.red,here("temp/alignments/myloplus-209.csv"))

# make labels and align
seqs.fas.red <- tab2fas(df=master.df.red,seqcol="nucleotidesFrag",namecol="dbidNex")
seqs.fas.red.ali <- as.matrix(ips::mafft(seqs.fas.red,exec="mafft"))

# trim the alignment
dim(seqs.fas.red.ali)
seqs.fas.red.ali.trim <- seqs.fas.red.ali[,44:664]
dim(seqs.fas.red.ali.trim)

# write out nexus
write.nexus.data(seqs.fas.red.ali.trim,here("temp/alignments/myloplus-209x621-aligned.nex"),interleaved=FALSE)
# write out fasta
write.FASTA(seqs.fas.red.ali.trim,here("temp/alignments/myloplus-209x621-aligned.fasta"))


### MAKE A QUICK TREE ###

# load data
seqs.fas <- read.FASTA("assets/sequences-master.fasta")
master.df <- read_csv(here("assets/tissues-master.csv"),show_col_types=FALSE)
# check
base::setdiff(labels(seqs.fas),pull(master.df,label))
base::setdiff(pull(master.df,label),labels(seqs.fas))

# align
myloplus.aligned <- as.matrix(ips::mafft(seqs.fas,exec="mafft"))

# make a tree
myloplus.tr <- nj(dist.dna(myloplus.aligned,model="TN93",pairwise.deletion=TRUE))
myloplus.tr.ml <- phangorn::optim.pml(phangorn::pml(myloplus.tr, phangorn::as.phyDat(myloplus.aligned), k=4, inv=0, model="HKY"), optNni=FALSE, optGamma=TRUE, optInv=FALSE, model="HKY")#rearrangement="stochastic",

# ladderize tree
myloplus.tr <- ladderize(midpoint(myloplus.tr.ml$tree))
myloplus.tr$edge.length[which(myloplus.tr$edge.length < 0)] <- 0

# make plotting labels
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
