#!/usr/bin/env Rscript
library(methods)
library(data.table)
library(optparse)

source("featurecountsFUN.R")

#Check the version of Rsubread
checkRsubreadVersion()

option_list = list(
  make_option(c("-b", "--bamfile"), type="character", default=NULL, 
              help="bam file", metavar="character"),
  make_option(c("-g", "--gtf"), type="character", default=NULL, 
              help="gtf file", metavar="character"),
  make_option(c("-c", "--n_cores"), type="integer", default=NULL, 
              help="number of cores", metavar="integer")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



##############################################################
##### featureCounts
print(Sys.time())
print("I am making SAF for featureCounts...")
saf <- makeSAF(opt$gtf)
abamfile = opt$bamfile

print(Sys.time())
print("I am counting exonic reads...")
fnex<-runFeatureCount(abamfile,
                       saf=saf$exons,
                       strand=0,
                       type="ex",
                       primaryOnly = TRUE,
                       cpu = opt$n_cores,
                       mem = 4)

print(Sys.time())
print("I am counting intronic reads...")

  fnin  <-runFeatureCount(abamfile,
                           saf=saf$introns,
                           strand=0,
                           type="in",
                           primaryOnly = TRUE,
                           cpu = opt$n_cores,
                           mem = 4)
