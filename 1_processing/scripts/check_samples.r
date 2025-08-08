suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(parallel, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(readxl, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(testit, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-i", "--input_samples_xlsx"),
              type="character",
              help="Path of XLSX file to process")
  )

opt <- parse_args(OptionParser(option_list=option_list))

# #Test scenario
# opt$input_samples_xlsx<-'/l/Zeitlinger/ZeitlingerLab/pipeline/users/ke2488/ke2488_samples.xlsx'

sample_types<-c('ChIP-seq','ChIP-nexus','MNase','ATAC-seq', 'scATACmulti')

samples.df<-readxl::read_xlsx(opt$input_samples_xlsx) %>%
  dplyr::mutate(bowtie2_path = paste0('/l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/', reference_genome, '/', reference_genome),
                starting_files_r2 = ifelse(is.na(starting_files_r2), 'X', starting_files_r2),
                nexus_barcodes = ifelse(is.na(nexus_barcodes), 'none', nexus_barcodes))

##########################################################################################
#Delete sample_assertion_check.txt
##########################################################################################
sample_assert_file<-paste0(dirname(opt$input_samples_xlsx), '/sample_assertion_check.txt')
system(paste0('rm -f ', sample_assert_file))

##########################################################################################
#Create assertions to make sure pipeline is run properly
##########################################################################################

testit::assert('There should be no ChIP-seq samples that are missing extension lengths',
               (samples.df %>% dplyr::filter(experiment=='ChIP-seq', sequencing_type=='single', is.na(extension_length)) %>% nrow(.))==0)

testit::assert('Your deduplicate section is not a valid value. Type either T or F.',
               (samples.df %>% dplyr::filter(!(deduplicate %in% c('T','F'))) %>% nrow(.))==0)

testit::assert('You have sequencing types that are not single or paired. Please check',
               (samples.df %>% dplyr::filter(!(sequencing_type %in% c('single','paired'))) %>% nrow(.))==0)

testit::assert('There should be no ChIP-nexus samples processed in paired-end format.',
               (samples.df %>% dplyr::filter(experiment=='ChIP-nexus', sequencing_type=='paired') %>% nrow(.))==0)

testit::assert('There should be no ATAC-seq samples processed in single-end format.',
               (samples.df %>% dplyr::filter(experiment=='ATAC-seq', sequencing_type=='single') %>% nrow(.))==0)

testit::assert('Your bowtie2 indexes do not exist based on the reference genome prefix listed. Please check to make sure they exist under `supplemental_data/indexes/bowtie2`.',
               all(samples.df$bowtie2_path %>% paste0(., '.1.bt2') %>% file.exists(.)))

testit::assert('Your sample names are not all unique. Please address.',
               length(samples.df$sample_name)==length(unique(samples.df$sample_name)))

testit::assert('Your read1 and read2 .fastq.gz files are the same. Please check your samples.df.',
               all(samples.df %>% dplyr::mutate(r1_vs_r2 = (starting_files_r1!=starting_files_r2)) %>% .$r1_vs_r2))

testit::assert('`output_sample` needs to have a _positive or a _negative under every ChIP-nexus field.',
               (samples.df %>% dplyr::filter(experiment=='ChIP-nexus', !grepl('positive|negative', output_sample)) %>% nrow(.))==0)

testit::assert('ChIP-nexus barcodes are required for pipeline processing, please fill in correctly.',
               (samples.df %>% dplyr::filter(experiment=='ChIP-nexus', nexus_barcodes == 'none') %>% nrow(.))==0)

testit::assert('Check your `starting_files_r2` for correct naming of your paired-end samples.df, there is likely missing starting files.',
               (samples.df %>% dplyr::filter(starting_files_r2=='X', sequencing_type == 'paired') %>% nrow(.))==0)

testit::assert('Experiment type is invalid in one of your rows, please check.',
               all(samples.df$experiment %in% sample_types))

#Check for trailing white spaces
testit::assert('You have trailing white spaces in your spreadsheet--this could lead to upload errors and therefore should be checked.',
               !any(samples.df %>% as.matrix() %>% as.character() %>% .[!is.na(.)] %>% stringr::str_sub(., -1, -1) ==' '))

#Check for WCE samples that aren't in spreadsheet.
testit::assert('The WCE samples listed under `control_name` are not actually in the spreadsheet. Please add to ensure pipeline runs properly.',
               (nrow(samples.df %>% dplyr::filter(!is.na(control_name)))==0) |
                 all(samples.df %>% dplyr::filter(!is.na(control_name)) %>% .$control_name %in% samples.df$sample_name))

testit::assert('Your sample read1 files need to exist',
               nrow(samples.df %>% dplyr::filter(!grepl('fastq|txt', starting_files_r1)))==0)

#scATACmulti assertions
testit::assert('You must only have one .fastq.gz file for R1 for this experiment type.',
               all((samples.df %>% dplyr::filter(experiment=='scATACmulti') %>% .$starting_files_r1 %>%
                 lapply(., function(x) length(Sys.glob(x))) %>% unlist())==1))
testit::assert('You must only have one .fastq.gz file for R2 for this experiment type.',
               all((samples.df %>% dplyr::filter(experiment=='scATACmulti') %>% .$starting_files_r2 %>%
                      lapply(., function(x) length(Sys.glob(x))) %>% unlist())==1))
testit::assert('This experiment requires _R1_ to have an equivalent naming convention for _I1_ indexes in the fastq files',
               all((samples.df %>% dplyr::filter(experiment=='scATACmulti') %>% .$starting_files_r1 %>%
                      lapply(., function(x) length(Sys.glob(gsub('_R1_', '_I1_', x)))) %>% unlist())==1))
testit::assert('This experiment requires _R2_ to have an equivalent naming convention for _I2_ indexes in the fastq files',
               all((samples.df %>% dplyr::filter(experiment=='scATACmulti') %>% .$starting_files_r2 %>%
                      lapply(., function(x) length(Sys.glob(gsub('_R2_', '_I2_', x)))) %>% unlist())==1))

#Allow that the samples pass asserts.
system(paste0('touch ', sample_assert_file))

##########################################################################################
# Check if samples.df need to be re-run and if so, delete samples.df
##########################################################################################

pipeline_columns <- c('sample_name', 'reference_genome', 'sequencing_type',
                    'experiment', 'project', 'deduplicate', 'extension_length', 'output_sample',
                    'read_1_3p_adapter', 'read_2_3p_adapter', 'nexus_barcodes', 'control_name',
                    'starting_files_r1', 'starting_files_r2')
key_dir <- '/l/Zeitlinger/ZeitlingerLab/pipeline/data/lab_data/logs/sample_records/'

#Generate pipeline-specific strings to denote each sample.
sample_keys<-lapply(1:nrow(samples.df), function(x){
  samples.df[pipeline_columns][x,] %>%
    as.character() %>%
    paste0(., collapse = '+') %>%
    gsub('/', '_', .)
}) %>% unlist

#Mark the sample key path.
samples.df<-samples.df %>%
  dplyr::mutate(key = sample_keys,
                key_path = paste0(key_dir, '/', sample_name))
samples.df$key <- sample_keys

#Are there already-existing samples.df which have been marked? If so, delete all files generated by them.
all_current_samples = Sys.glob(paste0(key_dir, '*'))
all_current_names = gsub(key_dir, '', all_current_samples)

filler<-lapply(1:nrow(samples.df), function(x){
  sample_status<-samples.df$sample_name[x] %in% all_current_names
  if(sample_status){
    old_sample_key<-samples.df$key_path[x] %>% read.table(.) %>% .$V1
    if(old_sample_key != samples.df$key[x]){
      message('Sample ', samples.df$sample_name[x], ' metadata changed...rerunning.')

      #Find all files
      sample_data_files <- Sys.glob(paste0('/l/Zeitlinger/ZeitlingerLab/pipeline/data/lab_data/*/*/', samples.df$output_sample[x], '.*'))
      output_sample_data_files <- Sys.glob(paste0('/l/Zeitlinger/ZeitlingerLab/pipeline/data/lab_data/*/*/', samples.df$sample_name[x], '.*'))
      log_files <- Sys.glob(paste0('/l/Zeitlinger/ZeitlingerLab/pipeline/data/lab_data/logs/*/*/', samples.df$sample_name[x], '_*'))

      if(length(c(sample_data_files, output_sample_data_files, log_files))!=0){
        all_files <-paste0(c(sample_data_files, output_sample_data_files, log_files) %>% unique(.), collapse = ' ')
        cmd_rm<-paste0('rm -f ', all_files, collapse = '')
        system(cmd_rm)
      }

    }
  }
  #Mark new/redundant/changed sample information
  cmd_echo<-paste0('echo ', samples.df$key[x], ' > ', samples.df$key_path[x], collapse = '')
  system(cmd_echo)
})

##########################################################################################
# Create slurm script
##########################################################################################

cmds<-c(
  '#!/bin/bash',
  '#SBATCH --job-name=zeitlinger_lab_pipeline',
  '#SBATCH --nodes=1',
  '#SBATCH --ntasks=1',
  '#SBATCH --cpus-per-task=8',
  '#SBATCH --mem=10gb',
  '#SBATCH --time=288:00:00',
  '#SBATCH --output=slurm_%j.log',
  'module load zeitlinger',
  paste0('export SNAKEMAKE_INPUT=', opt$input_samples_xlsx),
  'echo $SNAKEMAKE_INPUT',
  'if [ -f sample_assertion_check.txt ]; then',
  'snakemake --snakefile /l/Zeitlinger/ZeitlingerLab/pipeline/scripts/Snakefile --latency-wait 120 --profile slurm --resources disk_mb=500000',
  'fi '
)
writeLines(cmds, paste0(dirname(opt$input_samples_xlsx), '/run_snakemake.slurm'))
