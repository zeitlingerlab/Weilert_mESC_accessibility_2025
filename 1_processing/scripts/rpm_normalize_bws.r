suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-f", "--file"),
              type="character",
              help="Path to granges file"),
  make_option(c("-n", "--name"),
              type="character",
              help="name of the output data"),
  make_option(c("-t", "--type"),
              default = "atac",
              type="character",
              help="experiment type: nexus, atac or mnase"))

opt <- parse_args(OptionParser(option_list=option_list))

total_reads_normalization <- function(gr, type){

  #Standard RPM
  if(type %in% c('atac', 'mnase')){
    total_reads <- length(gr)
    gr.cov <- coverage(gr) / total_reads * 1000000
    gr.cov
  }

  #Read-specific RPM
  if(type == "chipnexus"){
    total_reads <- length(gr)
    gr.pos <- gr[strand(gr) == "+"]
    gr.neg <- gr[strand(gr) == "-"]
    gr.pos.cov <- coverage(gr.pos) / total_reads * 1000000
    gr.neg.cov <- coverage(gr.neg) / total_reads * 1000000 *(-1)
    cov.list <- list(pos = gr.pos.cov, neg = gr.neg.cov )
  }
}

if(opt$type %in% c('atac', 'mnase')){

  message("sample being processed is ATAC-seq")
  message("reading the granges file")
  gr <- readRDS(opt$file)

  message("generating the normalized coverage files")
  gr.cov <- total_reads_normalization(gr)

  message("exporting ATAC-seq bigwig files")
  export(gr.cov, paste0(opt$name, ".bw"))

}

if(opt$type == "nexus"){

	message("sample being processed is ChIP-nexus")
	message("reading the granges file")
	gr  <- readRDS(opt$file)
	gr <- resize(gr, 1, "start")

	message("generating the normalized coverage files")
	cov.list <- total_reads_normalization(gr, "nexus")

	message("exporting bigwig files")
	export(cov.list$pos, paste0(opt$name, "_positive.bw"))
	export(cov.list$neg, paste0(opt$name, "_negative.bw"))

}
