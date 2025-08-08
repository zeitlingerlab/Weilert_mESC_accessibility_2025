suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(GenomicAlignments, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(parallel, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-f", "--file"),
              type="character",
              default=NA,
              help="Path of BAM file to process"),
  make_option(c("-e", "--extension"),
              type="integer",
              default="native",
              help="Extension length (bp)"),
  make_option(c("-n", "--name"),
              type="character",
              help="R prefix for resulting coverage and ranges objects"),
  make_option(c("-b", "--bigwig"),
              type="character",
              default=NA,
              help="File prefix for resulting BigWig"),
  make_option(c("-p", "--paired"),
              action="store_true",
              default=FALSE,
              help="Paired-end data"),
  make_option(c("-d", "--deduplicate"),
              type="logical",
              default=TRUE,
              help="Deduplicate and apply filter for ChIP-seq artifacts"),
  make_option(c("-c", "--cores"),
              type="integer",
              default=2,
              help="Number of processor cores to use when performing artifact filtering")
  )

write_bigwig <- function(cov, filename) {
  message(id, "Writing bigWig...")
  export(cov, filename, format = "BigWig")
}

pn <- function(value) {
  prettyNum(value, big.mark=",")
}

# ------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list=option_list))

# opt$file<-"/l/Zeitlinger/ZeitlingerLab/pipeline/data/lab_data/bam/zelda_and_nucleosomes/dme_h3_orer_mbt_seq_1.bam"
# opt$extension<-183
# opt$name<-"test"
# opt$bigwig<-"~/tmp/test.bw"

if(is.na(opt$file)) {
  message("No BAM file specified. Use --help to list available options.")
  q(status=1)
}

suppressPackageStartupMessages(library(Rsamtools, warn.conflicts=F, quietly=T))
# suppressPackageStartupMessages(library(chipseq, warn.conflicts=F))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F))
suppressPackageStartupMessages(library(GenomicAlignments))

# used as a prefix for all output messages
id <- "[unknown] "

bam_file      <- opt$file
ext_length    <- opt$extension
var_name      <- opt$name

id <- paste("[", var_name, "] ", sep="")

bigwig_file   <- opt$bigwig
paired        <- opt$paired

if(is.na(bigwig_file)) bigwig_file <- paste(var_name, "bw", sep=".")

message(id, "Converting BAM to ranges object:")
message(id, "Input BAM: ", bam_file)
message(id, "Object name: ", var_name)

if(!file.exists(bam_file)) {
	stop("Could not open BAM file: ", bam_file)
}

if(paired) {
  if(opt$deduplicate){
    bam.gr <- granges(readGAlignmentPairs(opt$file, param=ScanBamParam(what = "qname", flag=scanBamFlag(isDuplicate=FALSE))))
  } else {
    bam.gr <- granges(readGAlignmentPairs(opt$file, param=ScanBamParam(what=scanBamWhat())))
  }
} else {
  if(opt$deduplicate){
    bam.gr <- granges(readGAlignments(opt$file, param=ScanBamParam(what = "qname", flag=scanBamFlag(isDuplicate=FALSE))))
  } else {
    bam.gr <- granges(readGAlignments(opt$file, param=ScanBamParam(what=scanBamWhat())))
  }

  message(id, "Resizing to fragment length...")
  bam.gr<-GenomicRanges::trim(GenomicRanges::resize(bam.gr, width = opt$extension, 'start'))
  # bam.gr<-GenomicRanges::resize(bam.gr, opt$extension, 'start')
}

message(id, "Saving ranges object...")
bam.gr <- bam.gr[order(bam.gr)]
saveRDS(bam.gr, paste(var_name, ".granges.rds", sep=""))

message(id, "Generating coverage object...")
sample.cov <- coverage(bam.gr)

nothing <- write_bigwig(sample.cov, bigwig_file)
