#!/bin/bash
#SBATCH --job-name=modiscolite
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=70
#SBATCH --mem=100gb
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%j.log
source /home/mw2098/.bashrc
conda deactivate
ml meme
conda activate bpreveal_404
cd /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis
/home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bpnet_osknz_fold1_sox2_counts.h5 \
--seqs npz/bpnet_osknz_fold1_sox2_counts_seqs \
--scores npz/bpnet_osknz_fold1_sox2_counts_scores
mkdir -p modiscolite/bpnet_osknz_fold1_sox2_counts
modisco motifs \
-s npz/bpnet_osknz_fold1_sox2_counts_seqs.npy \
-a npz/bpnet_osknz_fold1_sox2_counts_scores.npy \
-n 50000 \
-w 1000 \
--output modiscolite/bpnet_osknz_fold1_sox2_counts/modisco.h5
modisco report -i modiscolite/bpnet_osknz_fold1_sox2_counts/modisco.h5 \
-o modiscolite/bpnet_osknz_fold1_sox2_counts \
-m ../../public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
