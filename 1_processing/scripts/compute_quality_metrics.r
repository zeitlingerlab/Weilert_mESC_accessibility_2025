suppressPackageStartupMessages(library(optparse, warn.conflicts=F))
suppressPackageStartupMessages(library(readr, warn.conflicts=F))
suppressPackageStartupMessages(library(RCurl, warn.conflicts=F))
suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts=F))
suppressPackageStartupMessages(library(Rsamtools, warn.conflicts=F))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyranges))

option_list <- list(
  make_option(c("-f", "--fastq_reads_file"),
              type="character",
              help="Original FASTQ"),
  make_option(c("-s", "--bowtie2_statistics"),
              type="character",
              help=".txt file with bowtie2 alignment statistics (obtained by using `2>bowtie2_stats.txt` at the end of the bowtie2 commmand)."),
  make_option(c("-r", "--granges"),
             type="character",
             help="GRanges RData of unique reads"),
  make_option(c("-b", "--nexus_bam"),
              type="character",
              default = NA,
              help="ChIP-nexus .bam file for computing % sequenced"),
  make_option(c("-o", "--output_table"),
              type="character",
              help="Output filepath for .tsv file of performance metrics."))

opt <- parse_args(OptionParser(option_list=option_list))

# opt$fastq_reads_file<-'/l/Zeitlinger/ZeitlingerLab/pipeline/data/lab_data/logs/fastq_counts/fly_shavenbaby_collab/dme_svb_orer_12to14_nexus_4_raw_reads.txt'
# opt$bowtie2_statistics<-'/l/Zeitlinger/ZeitlingerLab/pipeline/data/lab_data/logs/bowtie_metrics/fly_shavenbaby_collab/dme_svb_orer_12to14_nexus_4_bowtie2_stats.txt'
# opt$granges<-'/l/Zeitlinger/ZeitlingerLab/pipeline/data/lab_data/rdata/fly_shavenbaby_collab/dme_svb_orer_12to14_nexus_4.granges.rds'
# opt$nexus_bam<-'/l/Zeitlinger/ZeitlingerLab/pipeline/data/lab_data/bam/fly_shavenbaby_collab/dme_svb_orer_12to14_nexus_4.bam'

short_name <- gsub("\\_bowtie2_stats.txt", "", basename(opt$bowtie2_statistics)) #get sample name
id <- paste0("[", short_name, "] ")
sid <- short_name

message(id, "Counting original FASTQ reads...")
original_fastq_count  <- as.integer(read.table(opt$fastq_reads_file)[[1]])

#Import GRanges.rds
gr<-opt$granges %>% readRDS(.)
deduplicated_read_count <- gr %>% length

message(id, "Parsing bowtie2 statistics for quality metrics...")
stats<-read_delim(opt$bowtie2_statistics, '\n', col_names = FALSE)$X1 %>% .[!grepl('Warning', .)]
processed_fastq_count<-stats[1] %>%
  gsub(pattern = " reads; of these:", replacement = "", .) %>%
  as.integer(.)
unaligned_read_count<-stats[3] %>%
  gsub(pattern = "\\(.*", replacement = "", .) %>%
  gsub(" ", "", .) %>%
  as.integer(.)
unique_aligned_read_count<-stats[4] %>%
  gsub(pattern = "\\(.*", replacement = "", .) %>%
  gsub(" ", "", .) %>%
  as.integer(.)

#Note that these multi-aligned reads are integrated into the reads kept for analysis.
multi_aligned_read_count<-stats[5] %>%
  gsub(pattern = "\\(.*", replacement = "", .) %>%
  gsub(" ", "", .) %>%
  as.integer(.)

#Import discordant reads that are kept
discordant_aligned_read_count<-stats[8] %>%
  gsub(pattern = "\\(.*", replacement = "", .) %>%
  gsub(" ", "", .) %>%
  as.integer(.)

total_aligned_read_count<- unique_aligned_read_count + multi_aligned_read_count + discordant_aligned_read_count
duplicated_read_count<-total_aligned_read_count-deduplicated_read_count

message("Calculating quality metrics...")
#Note: By design we will calculate quality metrics of coverage from a "seq-like" approach
#Note: In spite of ChIP-nexus being base resolution, these metrics will give better consensus for fragment-like features.
gr.cov<-coverage(gr)
dm <- as.data.frame(as(gr.cov, "GRanges"))

if (class(gr) == "CompressedIRangesList") {
  genome_size <- sum(sapply(gr.cov, function(x) as.numeric(length(x))))
} else {
  genome_size <- sum(as.numeric(seqlengths(gr)))
}

message("Calculating signal at 1pc...")
dm <- dm[order(dm$score, decreasing=T),]
dm$running_width_total <- cumsum(as.numeric(dm$width))
dm$running_score_total <- cumsum(dm$score * as.numeric(dm$width))
target_width <- as.integer(genome_size * 0.01)

message("Calculating genome coverage...")
iv <- which(dm$score == 0)
coverage <- sum(as.numeric(dm[-iv,]$width)) / sum(as.numeric(dm$width))

message(id, "Calculating uniqueness...")
uniqueness <- (1 - (duplicated_read_count / total_aligned_read_count))

#Format metrics
uniqueness<- round(uniqueness*100, 1)
coverage<- round(coverage*100, 1)

message(id, "Uniq: ", uniqueness)
message(id, "Coverage: ", coverage)


percent_sequenced<-NA
if(!is.na(opt$nexus_bam)){
  #Reading BAM file
  message(id, "Calculating % sequenced for ChIP-nexus samples...")
  message(id, "Reading BAM: ", opt$nexus_bam)
  bam.gr <- granges(readGAlignments(opt$nexus_bam, use.names=TRUE))
  bam.gr$barcode <- names(bam.gr)
  names(bam.gr) <- NULL
  bam_count <- length(bam.gr)

  #Obtain unique barcodes
  message("Obtaining unique barcodes:")
  unique_labeling<-paste(start(bam.gr), seqnames(bam.gr), bam.gr$barcode) #Identify unique labels

  #Find PCR Frequency Table
  message("Finding PCR frequency tables...")
  unique.labeling.freq<-table(unique_labeling)
  PCR_dup.freq<-table(unique.labeling.freq)
  PCR_dup.df<-data.frame(PCR_dup.freq, dups=as.numeric(names(PCR_dup.freq)), stringsAsFactors = F)[,2:3];
  names(PCR_dup.df)<-c("x_k","k")
  PCR_dup.df<-PCR_dup.df[,c(2,1)]

  #Extrapolate probability from poisson distribution
  message("Extrapolating metrics from frequency data:")
  x_1ton<-sum(PCR_dup.df$x_k)
  lambda<-2*PCR_dup.df$x_k[which(PCR_dup.df$k==2)]/PCR_dup.df$x_k[which(PCR_dup.df$k==1)] #obtain lambda
  p_0<-exp(-lambda)
  percent_sequenced<-1-p_0 #obtain probability

}
if(length(percent_sequenced)==0){percent_sequenced<-NA}

message("Metrics done uploading. Writing summary and exiting...")
output.df<-data.frame(
  sample_name = short_name,
  original_fastq_count = original_fastq_count,
  processed_fastq_count = processed_fastq_count,
  unaligned_read_count = unaligned_read_count,
  deduplicated_read_count = deduplicated_read_count,
  duplicated_read_count = duplicated_read_count,
  uniqueness = uniqueness,
  genome_coverage = coverage,
  percent_sequenced = percent_sequenced
)
readr::write_tsv(output.df, opt$output_table)
