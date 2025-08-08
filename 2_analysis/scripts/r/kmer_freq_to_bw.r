#Melanie Weilert
#March 3030
#Purpose: Extracting a kmer frequency across a genome and export as bigwig.
#kmer: any combination of nucleotides

#Computational setup
library(rtracklayer) ; library(GenomicRanges); library(magrittr) ; library(Biostrings)
library(plyranges); library(Rsamtools); library(data.table); library(optparse)
library(plyr); library(dplyr); library(stringr); library(parallel)

option_list <- list(
  make_option(c("-k", "--kmer"),
              type="character",
              help="kmer to search genome across."),
  make_option(c("-g", "--genome"),
              type="character",
              help="BSgenome object/package to look across"),
  make_option(c("-o", "--output_bw_filepath"),
              type="character",
              help="Output filepath for BW denoting hits across the genome"),
  make_option(c("-c", "--cores"),
              type="integer",
              default=4,
              help="Number of cores"))
opt <- parse_args(OptionParser(option_list=option_list))

#Test params
# opt$kmer = "ATTA"
# opt$genome = "BSgenome.Mmusculus.UCSC.mm10"
# opt$output_bw_filepath = "/n/projects/mw2098/genomes/mm10/kmers/test.bw"
# opt$cores=4

message("Make sure genome is installed...")
is.installed <- function(genome){
  is.element(genome, installed.packages()[,1])
}
if (!is.installed(opt$genome)){
  install.packages(opt$genome)
}
library(opt$genome,character.only = TRUE)

#Reassign genome as the actual object, not just the string
genome<-getBSgenome(opt$genome, masked=FALSE)

message("Get the locations of the kmer across the genome...")
locations_of_kmer_gr<-mclapply(seqnames(genome), function(x){
  chr_gr<-GRanges(seqnames=x, IRanges(start=1, end=seqlengths(genome)[x]), strand="*") #200M reads will take about 30seconds
  chr_seq<-getSeq(genome, chr_gr)
  chr_seq_str<-as.character(chr_seq)
  chr_locs<-gregexpr(pattern = paste0("(?=", opt$kmer, ")"), text = chr_seq_str, perl=TRUE) %>% unlist
    # str_locate_all(pattern = opt$kmer, string = chr_seq_str)[[1]][,"start"]
  if(length(chr_locs)>0){
    hits_gr<-GRanges(seqnames = x, IRanges(start=chr_locs, end=chr_locs), strand="+") #translate each of the nucleotide hits as reads and then generate coverage
    return(hits_gr)
  }
  else{ #Return null if no matches are found. This is really only applicable for small genomes and complex kmers.
    return(NULL)
  }
}, mc.cores=opt$cores)

#This GRanges represents the start site of the kmer (1-based coordinates) of interest.
locations_of_kmer_gr<-locations_of_kmer_gr[lapply(locations_of_kmer_gr, function(x) !is.null(x)) %>% unlist] %>% as(Class = "GRangesList") %>% unlist
seqlevels(locations_of_kmer_gr)<-seqlevels(genome)
seqinfo(locations_of_kmer_gr)<-seqinfo(genome)

message("Generate coverage...")
locations_of_kmer_cov<-coverage(locations_of_kmer_gr)

message("Export to bw...")
export(locations_of_kmer_cov, opt$output_bw_filepath, format="BigWig")

message("Goodbye! ;D")
