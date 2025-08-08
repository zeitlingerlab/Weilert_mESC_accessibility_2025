#!/usr/bin/env zsh
#SBATCH --job-name dumpFASTQ
#SBATCH --ntasks=65
#SBATCH --nodes=1
#SBATCH --mem=500gb
#SBATCH --time=5:00:00
#SBATCH --output=logs/dumpFASTQ_%A_%a.out
#SBATCH --partition=compute
#SBATCH --array=1-25%4

#Kudos: Charles McAnany's code

source /home/mw2098/.bashrc
module load bpreveal
module load bedtools
module load meme
module load sratoolkit

cd /n/projects/mw2098/publications/2024_weilert_acc/public/published/GSE174774

if [[ ${SLURM_ARRAY_TASK_ID} == 1 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607374 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 2 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607363 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 3 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607364 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 4 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607365 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 5 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607366 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 6 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607367 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 7 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607368 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 8 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607369 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 9 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607370 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 10 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607371 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 11 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607372 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 12 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607373 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 14 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607351 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 15 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607352 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 16 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607353 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 17 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607354 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 18 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607355 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 19 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607356 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 20 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607357 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 21 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607358 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 22 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607359 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 23 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607360 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 24 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607361 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi

if [[ ${SLURM_ARRAY_TASK_ID} == 25 ]] ; then
    /opt/apps/dev/containers/sratoolkit/3.0.1/bin/fasterq-dump SRR14607362 -O fastq --threads 65 --progress --mem 50000MB --bufsize 50000MB --curcache 50000MB 
fi










