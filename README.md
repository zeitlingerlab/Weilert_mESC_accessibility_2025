---
Author: Melanie Weilert
Date: August 2025
Purpose: Setup directions for accessibility cooperativity/low-affinity paper
---

# Introduction

For the work titled "Widespread low-affinity motifs enhance chromatin accessibility through intra-nucleosomal cooperativity in mESCs" from the Zeitlinger Lab, we give setup instructions for someone interested in reproducing analysis from this manuscript. 

For raw sequencing data, access the GEO referenced in the manuscript.
For intermediate files like trained models, TF-MoDISco results, motif coordinates, access the Zenodo referenced in the manuscript.
For code used to generated the intermediate files and final figures, refer to the scripts and code here. 

All raw data/intermediate data/code/etc is readily available upon request to the Zeitlinger lab and Melanie Weilert.

The steps we take are:

1. Environment setup (Anaconda3)
2. Process sequencing data (Snakemake)
3. Model optimization, training, and validation
4. Analysis (Python and R)

Please keep in mind that the Stowers Institute infrastructure relies on an HPC that runs SLURM jobs as well as allows for interactive sessions hosted using both CPUs and GPUs. Because of this, much of the work reported below will be designed under that workflow.

# Step 1: Acquire published raw data

To acquire GEO raw data to reprocess previously published samples, we used the following code:

```
sbatch code/0_setup/1_get_GSE174774_fastq.slurm
```

# Step 2: Build custom CRISPR genomes

In addition to setting up environments and processing samples from raw sequencing data, we also need to build our CRISPR genomes based on the custom edits imposed on the R1 mESC cell lines. See `code/0_setup/2_build_crispr_genomes.Rmd/html` for documentation on how these genomes were built for subsequent alignment and editing.

Keep in mind, this requires the `mm10.fa` file from UCSC to be imported and indexed at `code/0_setup/mm10.fa*`.

# Step 3: Set up modeling base environment

## GPU usage

This is the environment used to perform modeling and analysis in Python. We performed all deep learning modeling using a Keras/TensorFlow framework on NVIDIA® A100 TensorCore 40GB GPUs (CUDA v11.8).

## Clone the BPReveal repo

First, clone the BPReveal (v.4.0.4) software implementation of BPNet (https://github.com/mmtrebuchet/bpreveal):

`git clone https://github.com/mmtrebuchet/bpreveal -b 4.0.4`

## Create environment

Here, we create an Anaconda3 environment (`bpreveal_404`) by which we can process our scripts, train our models, and analyze our code. To do this, we will run the following script to set up our environment (originally created by Charles McAnany, modified to incorporate pipeline processing features).

`sbatch 0_setup/3_build_conda_env.slurm`

To run this script the following software needs to be in your environment:

+ python=3.11.7
+ conda==24.1.1
+ mamba==1.5.0

In brief, the versions of key packages are as follows for this environment:

+ pysam==0.22.0
+ pybedtools==0.9.1
+ pybigwig==0.3.22
+ numpy==1.26.4
+ pandas==2.2.0
+ matplotlib==3.8.3
+ plotnine==0.12.4
+ pathos==0.3.2
+ jupyterlab==4.1.1
+ nb_conda==2.2.1
+ cuda-toolkit==11.8.0
+ keras==2.15.0
+ tensorflow==2.15.0
+ tensorflow-probability==0.23.0
+ cmake==3.28.3
+ h5py==3.10.0
+ tqdm==4.66.2
+ gxx_linux-64==13.2.0
+ gfortran==13.2.0

For thermodynamics/kinetics work, an additional package was required:

+ sympy==1.13.3

# Step 4-5: Processing environment

This is the environment used to align and process genomics data based on the given `Snakefile`. This is part of a container set up by Jonathon Russell from the Stowers Institute to manage the custom Zeitlinger lab pipeline environment requirements.

All experimental samples and reprocessed samples from GEO are listed under: `/n/projects/mw2098/publications/2024_weilert_acc/code/0_setup/4_define_samples.xlsx`. We ran a `snakemake` pipeline to process all samples. The command was run to allocated SLURM queued jobs using the following script:

```
sbatch code/0_setup/5_process_samples.slurm
```

Upon processing and storage in `code/1_processing`, data was then modeled and analyzed in `code/2_analysis`.

## Versions

In brief, the versions of key packages are as follows for this environment:

+ python=3.11.5
  + os
  + math
  + itertools
  + csv=1.0
  + numpy=1.24.3
  + pandas=2.0.3
  + openpyxl=3.0.10
+ R=4.2.3
  + optparse_1.7.3
  + testit_0.13
  + parallel_4.2.3
  + dplyr_1.1.2
  + plyranges_1.18.0
  + readxl_1.4.3
  + readr_2.1.4
  + Rsamtools_2.14.0
  + BiocParallel_1.32.6
  + RCurl_1.98-1.12
  + ShortRead_1.56.1
  + GenomicAlignments_1.34.1
  + GenomicRanges_1.50.2
  + magrittr_2.0.3
  + rtracklayer_1.58.0
  + Rcpp_1.0.11
  + stringr_1.5.1
+ ucsc
+ fastx-toolkit=0.0.14
+ cutadapt=4.2
+ macs2=2.2.7.1
+ samtools=1.15.1
+ bowtie2=2.3.5.1
+ zcat=1.9
+ snakemake=7.32.4
+ STAR=2.7.3a

### Supplementary software

Other packages are required for processing the pipeline:
+ Java OpenJDK==1.8.0_191
+ PICARD==2.23.8
+ meme==5.5.3

Additionally, IDR was performed for reproducible peak identification in the downstream work in R markdowns.
+ idr==2.0.3

SciPy was used for null_matrix resolution in our kinetics models:
+ scipy=1.12.0

# Step 6: Analysis and modeling

Analysis was conducted in both Python and R under `code/2_analysis`. For Python, analysis environment is the same as the modeling environment described below. For R, `sessionInfo()` is reported with relevant packages and versions at the bottom of each markdown. Markdown files are numbered in the order by which they were run. Raw figures can be found here as well as code and associated scripts to run analysis.