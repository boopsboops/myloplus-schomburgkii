#!/usr/bin/env Rscript

library("here")
library("ape")
library("spider")
library("tidyverse")
library("magrittr")
library("ips")
library("phangorn")
library("ggtree")
library("traits")
# renv::install("boopsboops/traits@b62d0431b34d119038e23cec125bb46e8df52e13")
library("lubridate")
library("glue")
library("parallel")

# load funs
source("https://raw.githubusercontent.com/boopsboops/UTILITIES/master/RScripts/hapCollapse.R")
source("https://raw.githubusercontent.com/boopsboops/UTILITIES/master/RScripts/tab2fas.R")
source("https://raw.githubusercontent.com/boopsboops/UTILITIES/master/RScripts/haps2fas.R")
source("https://raw.githubusercontent.com/boopsboops/UTILITIES/master/RScripts/hap_collapse_df.R")
source("https://raw.githubusercontent.com/boopsboops/UTILITIES/master/RScripts/get_sames.R")
