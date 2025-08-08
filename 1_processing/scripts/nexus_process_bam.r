suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

pn <- function(value) {
  prettyNum(value, big.mark=",")
}

option_list <- list(
  make_option(c("-f", "--file"),
              type="character",
              help="Path of BAM file to process"),
  make_option(c("-n", "--rdata_output"),
              type="character",
              help="Name for resulting output file objects"),
  make_option(c("-b", "--bigwig_output"),
              type="character",
              help="Name for resulting output file of .bigwig"),
  make_option(c("-p", "--paired"),
              action="store_true",
              default=FALSE,
              help="Paired-end mode"),
  make_option(c("-u", "--unique"),
              type="logical",
		          default=TRUE,
		          help="remove duplicates based on random barcode"))

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(rtracklayer))

id <- paste0("[", opt$rdata_output, "] ")

stopifnot(file.exists(opt$file))

if(opt$paired) {
  message(id, "Reading paired-end BAM: ", opt$file)
  bam <- readGAlignmentPairs(opt$file, param=ScanBamParam(what="qname"))
  bam.gr <- granges(bam)
  bam.gr$barcode <- gsub("^(.*)_pair_.*$", "\\1", mcols(first(bam))$qname)
} else {
  message(id, "Reading BAM: ", opt$file)
  bam.gr <- granges(readGAlignments(opt$file, use.names=TRUE))
  bam.gr$barcode <- names(bam.gr)
  names(bam.gr) <- NULL
}

message(id, pn(length(bam.gr)), " fragments")

if(opt$unique){
	message(id, "Removing barcode duplicates...")
	bam_uniq.gr <- GenomicRanges::split(bam.gr, bam.gr$barcode) %>%
	               unique %>%
	               unlist(use.names=FALSE)
}else{
	bam_uniq.gr <- bam.gr
}


message(id, "Saving ", pn(length(bam_uniq.gr)), " fragments...")
saveRDS(bam_uniq.gr, file=paste(opt$rdata_output, ".granges.rds", sep=""))

message(id, "Processing to .bw ...")
gr <- resize(bam_uniq.gr, 1, "start")
gr.pos <- gr[strand(gr) == "+"]
gr.neg <- gr[strand(gr) == "-"]
gr.pos.cov <- coverage(gr.pos)
gr.neg.cov <- coverage(gr.neg) * (-1)

message(id, "Saving to .bw ...")
export(gr.pos.cov, paste0(opt$bigwig_output, "_positive.bw"))
export(gr.neg.cov, paste0(opt$bigwig_output, "_negative.bw"))
