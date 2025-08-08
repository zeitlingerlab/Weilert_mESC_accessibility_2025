#!/bin/bash
#SBATCH --job-name=modiscolite
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%j.log
source /home/mw2098/.bashrc
conda deactivate
ml meme
conda activate acc_vs_bind
cd /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis
/home/mw2098/anaconda3/envs/acc_vs_bind/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal//src/shapToNumpy.py \
--h5 shap/atac_bias_bias_atac_counts.h5 \
--seqs npz/atac_bias_bias_atac_counts_seqs \
--scores npz/atac_bias_bias_atac_counts_scores
mkdir -p modiscolite/atac_bias_bias_atac_counts
modisco motifs \
-s npz/atac_bias_bias_atac_counts_seqs.npy \
-a npz/atac_bias_bias_atac_counts_scores.npy \
-n 50000 \
-w 1000 \
--output modiscolite/atac_bias_bias_atac_counts/modisco.h5
modisco report -i modiscolite/atac_bias_bias_atac_counts/modisco.h5 \
-o modiscolite/atac_bias_bias_atac_counts \
-m ../../public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
