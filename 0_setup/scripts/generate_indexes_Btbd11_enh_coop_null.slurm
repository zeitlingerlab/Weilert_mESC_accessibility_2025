#!/bin/bash
#SBATCH --job-name=bowtie_indexes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100gb
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%j.log
source /home/mw2098/.bashrc
module load zeitlinger
cd /n/projects/mw2098/publications/2024_weilert_acc/code/0_setup/
bowtie2-build fa/Btbd11_enh_coop_null.fa ../1_processing/indexes/bowtie2/Btbd11_enh_coop_null
