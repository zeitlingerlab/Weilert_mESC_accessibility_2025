#Melanie Weilert
#March 2020
#Purpose: Merge bigwig files using the "sum" function
#Note: Bws with overlapping regions will error out of this function.

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(testit))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(parallel))

option_list <- list(
  make_option(c("-f", "--files"),
              type="character",
              help="Path of bw files with wildcard (*) representing the samples to merge. e.g. test_chr*.Dl.contrib.counts.bw"),
  make_option(c("-o", "--output"),
              type="character",
              help="Output bw file name."),
  make_option(c("-g", "--bsgenome"),
              type="character",
              help="BSGenome object"),
  make_option(c("-c", "--cores"),
              type="integer", default = 6,
              help="Number of cores to use when loading multiple bws."))

opt <- parse_args(OptionParser(option_list=option_list))

# opt$files<-"/l/Zeitlinger/ZeitlingerLab/data/bigwigs/chip_next_development/mesc_fixed_homemade_500k_*vxl_atac_[0-9].bw"
# opt$bsgenome<-"BSgenome.Dmelanogaster.UCSC.dm6"
# opt$output<-"/n/projects/mw2098/analysis/chipnext_development/bw/combined/mesc_fixed_homemade_500k_atac_combined.bw"

message("Collecting files...")
bw.files<-system(paste0("ls ", opt$files), intern = TRUE)

message("Assigning genome...")
genome<-getBSgenome(opt$bsgenome)

message("Importing files and formatting to .cov...")
gr.list<-mclapply(bw.files, function(x){
  bw.gr<-import(x)
  seqlevels(bw.gr)<-seqlevels(genome)
  seqinfo(bw.gr)<-seqinfo(genome)
  return(bw.gr)
}, mc.cores = opt$cores)

#
# gr<-gr.list %>% as(Class = "GRangesList") %>% unlist
# assert("Input BW files have overlapping ranges. Please resolve these overlapping ranges, as they will be incorrectly summed upon execution of this script.",
#        length(reduce(gr, min.gapwidth = 0L))==length(gr))

cov.list<-mclapply(gr.list, function(x){
  bw.cov<-coverage(x, weight = x$score)
  return(bw.cov)
}, mc.cores = opt$cores)

message("Writing to .bw...")
merged_cov<-Reduce(f = "+", x = cov.list)
export(merged_cov, opt$output, format = "BigWig")
