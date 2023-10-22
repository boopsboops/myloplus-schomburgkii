#!/usr/bin/env Rscript

### LOAD LIBS ###
source(here::here("scripts/load-libs.R"))

# set output dir
outdir <- as.character(Sys.Date())
dir.create(here("temp-local-only",outdir))

# copy alignment
file.copy(from=here("temp/alignments/myloplus-209x621-aligned.fasta"),to=here("temp-local-only",outdir,"myloplus-209x621-aligned.fasta"))


### RUN RAXML NG ###

#load df
master.df.red <- read_csv(here("temp/alignments/myloplus-209.csv"),show_col_types=FALSE)

# run rax
rax.tr <- raxml_align(file=here("temp-local-only",outdir,"myloplus-209x621-aligned.fasta"),align=FALSE,model="GTR+G",epsilon=0.01)
rax.tr.root <- phangorn::midpoint(rax.tr)

# plot labels
master.df.ml <- master.df.red %>% 
    mutate(labsPhy=paste(label,sciNameValid,waterBody,nHaps,sep="|")) %>%
    select(dbidNex,labsPhy)

# plot tree
p <- rax.tr.root %>% 
    ggtree(color="grey50") %<+% master.df.ml +
    geom_tiplab(aes(label=labsPhy),geom="text",offset=0.0005,align=TRUE,color="grey50") + 
    theme(legend.position="none") +
    xlim(0,0.4)

# plot
filename <- glue("temp/trees/myloplus.tr.ml.",outdir,".pdf")
ggsave(filename=here(filename),plot=p,width=297,height=1000,units="mm",limitsize=FALSE)


### RUN MRBAYES IN TERMINAL ### 

# mb commands
# mpirun -np 4 mb myloplus-209x621.mb
# chain ((1000000/360)*4)*0.9 = 10k
# burnin ((2779-280)*4) = 10k

### LOAD MRBAYES TREES AND SUBSAMPLE ###

# test
#read_t(basename="temp-local-only/2022-11-10/myloplus-209x621-aligned.nex.run",run=1,burnin=102)

# load four runs
nruns <- 1:4
# burnin
all.runs <- mcmapply(function(x) read_t(basename=paste("temp-local-only",outdir,"myloplus-209x621-aligned.nex.run",sep="/"),run=x,burnin=280),x=nruns,SIMPLIFY=FALSE,USE.NAMES=FALSE,mc.cores=4)

# join
all.runs.joined <- do.call(c,all.runs)

# root
all.runs.joined.rooted <- mcmapply(function(x) phangorn::midpoint(x), x=all.runs.joined,SIMPLIFY=FALSE,USE.NAMES=FALSE,mc.cores=4)

# write out 10k trees for annotate
ape::write.nexus(all.runs.joined.rooted,file=here("temp-local-only",outdir,"all.runs.joined.rooted.trees"))

# get paths for consensus
trees.path <- here("temp-local-only",outdir,"all.runs.joined.rooted.trees")
mcc.path <- here("temp-local-only",outdir,"all.runs.joined.rooted.mcc.tre")

# run
tree_annotator(infile=trees.path,outfile=mcc.path)

# plot to check mcc

# read mcc tree
master.df.red <- read_csv(here("temp/alignments/myloplus-209.csv"),show_col_types=FALSE)
mcc.tre <- read.nexus(here("temp-local-only",outdir,"all.runs.joined.rooted.mcc.tre"))
# root
#mcc.tre <- ape::ladderize(phangorn::midpoint(mcc.tre))

# format labels
master.df.mcc <- master.df.red %>% 
    mutate(labsPhy=paste(label,sciNameValid,waterBody,nHaps,sep="|")) %>%
    select(dbidNex,labsPhy)

# ggtree to check topology
p <- mcc.tre %>% 
    ggtree(color="grey50",ladderize=TRUE) %<+% master.df.mcc +
    geom_tiplab(aes(label=labsPhy),geom="text",align=TRUE,color="grey50",offset=0.005,size=3.5) +
    theme(legend.position="none") +
    xlim(0,0.25)

# plot
filename <- glue("temp/trees/myloplus.tr.mcc.",outdir,".default.pdf")
ggsave(filename=here(filename),plot=p,width=297,height=1000,units="mm",limitsize=FALSE)

# sample
set.seed(42)
all.runs.joined.rooted.sample <- sample(all.runs.joined.rooted,1000)

# name the trees
trees.names <- paste0("tree",str_pad(1:length(all.runs.joined.rooted.sample),width=4,pad="0"))

# write out the trees
dir.create(here("temp-local-only",outdir,"mptp"))
mapply(function(x,y) ape::write.tree(x,file=here("temp-local-only",outdir,"mptp",paste0(y,".nwk"))),x=all.runs.joined.rooted.sample,y=trees.names)

# also write out a multiphylo to view
class(all.runs.joined.rooted.sample) <- "multiPhylo"
ape::write.tree(all.runs.joined.rooted.sample,file=here("temp-local-only",outdir,"all.runs.joined.rooted.sample.nwk"))
ape::write.nexus(all.runs.joined.rooted.sample,file=here("temp-local-only",outdir,"all.runs.joined.rooted.sample.trees"))



### GET MPTP FOR MCC TREE ###

# load up reduced annotated df and mcc tree
master.df.red <- read_csv(here("temp/alignments/myloplus-209.csv"),show_col_types=FALSE)
mcc.tre <- read.nexus(here("temp-local-only",outdir,"all.runs.joined.rooted.mcc.tre"))

# save as nwk
ape::write.tree(mcc.tre,file=here("temp-local-only",outdir,"mcc.nwk"))

# get mptp for consensus tree and load data
run_mptp(file=here("temp-local-only",outdir,"mcc.nwk"),threshold="single",minbr=0)#minbr=0.0001
mcc.mptp.out <- read_mptp(file=here("temp-local-only",outdir,"mcc.nwk.mptp.out.txt"))

# get hashes of sets
mcc.mptp.out.hash <- mcc.mptp.out %>% 
    arrange(rep,mptpDelim,label) %>% 
    group_by(rep,mptpDelim) %>% 
    mutate(labelCat=paste(label,collapse="|"),labelHash=openssl::md5(labelCat)) %>%
    ungroup() %>%
    select(-labelCat)


### GET MPTP FOR POSTERIOR TREE SAMPLE ###

# run mptp in parallel
# get tree names
mptp.nwk.trees <- list.files(here("temp-local-only",outdir,"mptp"),pattern=".nwk",full.names=TRUE)

# run parallel
mcmapply(function(x) run_mptp(file=x,threshold="single",minbr=0),x=mptp.nwk.trees,mc.cores=4)#minbr=0.0001

# read in parallel
mptp.nwk.out <- list.files(here("temp-local-only",outdir,"mptp"),pattern=".txt",full.names=TRUE)
mptp.nwk.out.df <- mcmapply(function(x) read_mptp(file=x),x=mptp.nwk.out,SIMPLIFY=FALSE,USE.NAMES=FALSE,mc.cores=4)

# join dataframes
mptp.nwk.out.df.joined <- bind_rows(mptp.nwk.out.df)

# hash the labels (slow)
mptp.nwk.out.df.joined.hash <- mptp.nwk.out.df.joined %>% 
    arrange(rep,mptpDelim,label) %>% 
    group_by(rep,mptpDelim) %>% 
    mutate(labelCat=paste(label,collapse="|"),labelHash=openssl::md5(labelCat)) %>%
    ungroup() %>%
    select(-labelCat)

# get summary
mptp.nwk.out.df.joined.hash.summary <- mptp.nwk.out.df.joined.hash %>% 
    distinct(rep,labelHash) %>% 
    group_by(rep) %>% 
    mutate(nDelim=length(unique(labelHash))) %>%
    ungroup()

# plot n delims
p <- mptp.nwk.out.df.joined.hash.summary %>% 
    distinct(rep,nDelim) %>%
    ggplot(aes(nDelim)) + 
        #geom_histogram(binwidth=1)
        geom_density()
# plot(p)


### GET POST CLADE PROBS ###
post.probs <- mptp.nwk.out.df.joined.hash.summary %>% 
     count(labelHash,sort=TRUE) %>% 
     mutate(postProb=n/1000) %>% 
     select(-n)

# join to mcc
post.probs.label <- mcc.mptp.out.hash %>% 
    left_join(post.probs) %>% 
    select(label,labelHash,postProb) %>%
    rename(dbidNex=label)

# join to metadata df
master.df.red.probs <- master.df.red %>% left_join(post.probs.label)
glimpse(master.df.red.probs)
master.df.red.probs %>% write_csv(here("temp/alignments/myloplus-209-mptp.csv"))


### PLOT WITH APLOT ###

# load up previously saved data
master.df.red.probs <- read_csv(here("temp/alignments/myloplus-209-mptp.csv"),show_col_types=FALSE)

# load tree
mcc.tre.beast <- treeio::read.beast(here("temp-local-only",outdir,"all.runs.joined.rooted.mcc.tre"))

# make df for plotting
master.df.plot <- master.df.red.probs %>% 
    mutate(labsPhy=paste(label,sciNameValid,waterBody,nHaps,sep="|")) %>%
    select(dbidNex,labsPhy,labelHash,postProb)

# colours with randomcoloR
ncols <- master.df.plot %>% distinct(labelHash) %>% pull() %>% length()
set.seed(11)
cols <- randomcoloR::distinctColorPalette(k=ncols)

# tree
# https://yulab-smu.top/treedata-book/
mcc.p <- mcc.tre.beast %>% 
    ggtree(color="grey40",ladderize=TRUE,right=TRUE) %<+% master.df.plot +
    geom_tiplab(aes(label=labsPhy),geom="text",align=TRUE,color="grey40",offset=0.005,size=4) +
    geom_tippoint(aes(color=as.character(labelHash)),size=2.5) +
    geom_nodepoint(aes(subset=posterior > 0.95),color="grey40",size=1) +
    scale_color_manual(values=cols) +
    theme(legend.position="none") +
    xlim(0,0.24)

# bar plot
post.probs.plot <- master.df.plot %>% 
    ggplot(aes(y=postProb,x=dbidNex)) + 
    geom_col(fill="grey90") + 
    geom_point(aes(color=as.character(labelHash)),shape=15,size=2.5) + 
    scale_color_manual(values=cols) +
    coord_flip() +
    theme_minimal() +
    scale_y_continuous(sec.axis=dup_axis()) +
    ylab("Posterior probability") +
    theme(legend.position="none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_line(size=0.5,linetype=2)
        )

# join plots
plots.joined <- post.probs.plot %>% aplot::insert_left(mcc.p,width=5)

filename <- glue("temp/trees/myloplus.tr.mptp.",outdir,".pdf")
ggsave(filename=here(filename),plot=plots.joined,width=297,height=1000,units="mm",limitsize=FALSE)


### HYPOTHESES SUPPORT ###

# add species clade groups by hand
# load table of clades
clades.df <- read_csv(here("assets/clades-test.csv"),show_col_types=FALSE)

# add post probs
clades.pp <- clades.df %>%
    arrange(clade,dbidNex) %>% 
    group_by(clade) %>% 
    mutate(labelCat=paste(dbidNex,collapse="|"),labelHash=openssl::md5(labelCat)) %>%
    ungroup() %>%
    distinct(clade,labelHash) %>%
    left_join(post.probs) %>%
    select(clade,postProb)

print(clades.pp)
clades.pp %>% write_csv(here("temp/alignments/clades-results.csv"))

# get post clade probs for each group 
# load trees
trees.path <- here("temp-local-only",outdir,"all.runs.joined.rooted.trees")
all.runs.joined.rooted <- read.nexus(trees.path)

# test on one
#group_monophyly(df=clades.df,trs=all.runs.joined.rooted,topo="guiana-shield-subset")

# run over all trees and clades
clades.tibs <- mapply(function(x) group_monophyly(df=clades.df,trs=all.runs.joined.rooted,topo=x,cores=4), x=pull(distinct(clades.df,clade),clade), SIMPLIFY=FALSE,USE.NAMES=FALSE)

# combine
groups.pp <- bind_rows(clades.tibs)
post.delim.clade.probs <- clades.pp %>% left_join(groups.pp) %>% rename(postDelimProb=postProb)
print(post.delim.clade.probs)
post.delim.clade.probs %>% write_csv(here("temp/alignments/groups-results.csv"))
