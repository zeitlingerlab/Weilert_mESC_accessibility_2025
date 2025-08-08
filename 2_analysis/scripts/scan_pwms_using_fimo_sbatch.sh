#!/bin/bash
#SBATCH --job-name=mw2098-2024-02-26
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32GB
#SBATCH --output=mw2098-2024-02-26.out

#!bin/bash
module load meme
fimo --oc memesuite --skip-matched-sequence --parse-genomic-coord --max-strand --thresh 0.001 --max-stored-scores 10000000 memesuite/custom_meme.txt memesuite/all_islands_curated_seqs.fa > memesuite/fimo_matches.tsv
