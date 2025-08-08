#!/bin/bash
#SBATCH --job-name=process_data
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
#SBATCH --time=288:00:00
#SBATCH --output=slurm_%j.log
source ~/.bashrc
module load zeitlinger
module load star/2.7.3a
cd /n/projects/mw2098/publications/2024_weilert_acc/code/1_processing
ulimit -n 8092
snakemake --latency-wait 120 --profile slurm --resources disk_mb=500000
