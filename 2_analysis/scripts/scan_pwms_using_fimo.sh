#!bin/bash
module load meme
fimo --oc memesuite --skip-matched-sequence --parse-genomic-coord --max-strand --thresh 0.001 --max-stored-scores 10000000 memesuite/custom_meme.txt memesuite/all_islands_curated_seqs.fa > memesuite/fimo_matches.tsv
