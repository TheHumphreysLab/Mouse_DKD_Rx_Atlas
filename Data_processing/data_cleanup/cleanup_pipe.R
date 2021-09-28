#!/usr/bin/env Rscript
library(Seurat)
library(ggplot2)
library(biomaRt)
library(DoubletFinder)
library(methods)
library(data.table)
library(optparse)
library(Matrix)
library(hdf5r)
source('data_process_funx.R')
option_list = list(
  make_option(c("-r", "--rdsfile"), type="character", default=NULL, 
              help="The rds file from CellBender", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="The output directory", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=NULL, 
              help="The sample ID", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
data_processing(cellbender_h5 = opt$rdsfile, sampleID = opt$name, out_dir = opt$outdir)
