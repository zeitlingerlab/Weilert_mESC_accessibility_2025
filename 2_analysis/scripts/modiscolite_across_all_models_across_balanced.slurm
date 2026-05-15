#!/usr/bin/bash
#SBATCH --array=1-18%15
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
if [[ ${SLURM_ARRAY_TASK_ID} == 1 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_genome_fold1_atac_counts_across_balanced.h5 \
--seqs npz/atac_wt_across_genome_fold1_atac_counts_seqs_across_balanced \
--scores npz/atac_wt_across_genome_fold1_atac_counts_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_genome_fold1_atac_counts
modisco motifs \
-s npz/atac_wt_across_genome_fold1_atac_counts_seqs_across_balanced.npy \
-a npz/atac_wt_across_genome_fold1_atac_counts_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_genome_fold1_atac_counts_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_genome_fold1_atac_counts_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_genome_fold1_atac_counts_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 2 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_genome_fold1_atac_profile_across_balanced.h5 \
--seqs npz/atac_wt_across_genome_fold1_atac_profile_seqs_across_balanced \
--scores npz/atac_wt_across_genome_fold1_atac_profile_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_genome_fold1_atac_profile
modisco motifs \
-s npz/atac_wt_across_genome_fold1_atac_profile_seqs_across_balanced.npy \
-a npz/atac_wt_across_genome_fold1_atac_profile_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_genome_fold1_atac_profile_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_genome_fold1_atac_profile_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_genome_fold1_atac_profile_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 3 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_genome_fold2_atac_counts_across_balanced.h5 \
--seqs npz/atac_wt_across_genome_fold2_atac_counts_seqs_across_balanced \
--scores npz/atac_wt_across_genome_fold2_atac_counts_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_genome_fold2_atac_counts
modisco motifs \
-s npz/atac_wt_across_genome_fold2_atac_counts_seqs_across_balanced.npy \
-a npz/atac_wt_across_genome_fold2_atac_counts_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_genome_fold2_atac_counts_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_genome_fold2_atac_counts_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_genome_fold2_atac_counts_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 4 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_genome_fold2_atac_profile_across_balanced.h5 \
--seqs npz/atac_wt_across_genome_fold2_atac_profile_seqs_across_balanced \
--scores npz/atac_wt_across_genome_fold2_atac_profile_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_genome_fold2_atac_profile
modisco motifs \
-s npz/atac_wt_across_genome_fold2_atac_profile_seqs_across_balanced.npy \
-a npz/atac_wt_across_genome_fold2_atac_profile_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_genome_fold2_atac_profile_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_genome_fold2_atac_profile_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_genome_fold2_atac_profile_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 5 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_genome_fold3_atac_counts_across_balanced.h5 \
--seqs npz/atac_wt_across_genome_fold3_atac_counts_seqs_across_balanced \
--scores npz/atac_wt_across_genome_fold3_atac_counts_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_genome_fold3_atac_counts
modisco motifs \
-s npz/atac_wt_across_genome_fold3_atac_counts_seqs_across_balanced.npy \
-a npz/atac_wt_across_genome_fold3_atac_counts_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_genome_fold3_atac_counts_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_genome_fold3_atac_counts_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_genome_fold3_atac_counts_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 6 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_genome_fold3_atac_profile_across_balanced.h5 \
--seqs npz/atac_wt_across_genome_fold3_atac_profile_seqs_across_balanced \
--scores npz/atac_wt_across_genome_fold3_atac_profile_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_genome_fold3_atac_profile
modisco motifs \
-s npz/atac_wt_across_genome_fold3_atac_profile_seqs_across_balanced.npy \
-a npz/atac_wt_across_genome_fold3_atac_profile_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_genome_fold3_atac_profile_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_genome_fold3_atac_profile_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_genome_fold3_atac_profile_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 7 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_split_fold1_atac_counts_across_balanced.h5 \
--seqs npz/atac_wt_across_split_fold1_atac_counts_seqs_across_balanced \
--scores npz/atac_wt_across_split_fold1_atac_counts_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_split_fold1_atac_counts
modisco motifs \
-s npz/atac_wt_across_split_fold1_atac_counts_seqs_across_balanced.npy \
-a npz/atac_wt_across_split_fold1_atac_counts_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_split_fold1_atac_counts_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_split_fold1_atac_counts_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_split_fold1_atac_counts_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 8 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_split_fold1_atac_profile_across_balanced.h5 \
--seqs npz/atac_wt_across_split_fold1_atac_profile_seqs_across_balanced \
--scores npz/atac_wt_across_split_fold1_atac_profile_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_split_fold1_atac_profile
modisco motifs \
-s npz/atac_wt_across_split_fold1_atac_profile_seqs_across_balanced.npy \
-a npz/atac_wt_across_split_fold1_atac_profile_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_split_fold1_atac_profile_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_split_fold1_atac_profile_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_split_fold1_atac_profile_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 9 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_split_fold2_atac_counts_across_balanced.h5 \
--seqs npz/atac_wt_across_split_fold2_atac_counts_seqs_across_balanced \
--scores npz/atac_wt_across_split_fold2_atac_counts_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_split_fold2_atac_counts
modisco motifs \
-s npz/atac_wt_across_split_fold2_atac_counts_seqs_across_balanced.npy \
-a npz/atac_wt_across_split_fold2_atac_counts_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_split_fold2_atac_counts_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_split_fold2_atac_counts_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_split_fold2_atac_counts_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 10 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_split_fold2_atac_profile_across_balanced.h5 \
--seqs npz/atac_wt_across_split_fold2_atac_profile_seqs_across_balanced \
--scores npz/atac_wt_across_split_fold2_atac_profile_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_split_fold2_atac_profile
modisco motifs \
-s npz/atac_wt_across_split_fold2_atac_profile_seqs_across_balanced.npy \
-a npz/atac_wt_across_split_fold2_atac_profile_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_split_fold2_atac_profile_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_split_fold2_atac_profile_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_split_fold2_atac_profile_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 11 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_split_fold3_atac_counts_across_balanced.h5 \
--seqs npz/atac_wt_across_split_fold3_atac_counts_seqs_across_balanced \
--scores npz/atac_wt_across_split_fold3_atac_counts_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_split_fold3_atac_counts
modisco motifs \
-s npz/atac_wt_across_split_fold3_atac_counts_seqs_across_balanced.npy \
-a npz/atac_wt_across_split_fold3_atac_counts_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_split_fold3_atac_counts_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_split_fold3_atac_counts_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_split_fold3_atac_counts_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 12 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_split_fold3_atac_profile_across_balanced.h5 \
--seqs npz/atac_wt_across_split_fold3_atac_profile_seqs_across_balanced \
--scores npz/atac_wt_across_split_fold3_atac_profile_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_split_fold3_atac_profile
modisco motifs \
-s npz/atac_wt_across_split_fold3_atac_profile_seqs_across_balanced.npy \
-a npz/atac_wt_across_split_fold3_atac_profile_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_split_fold3_atac_profile_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_split_fold3_atac_profile_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_split_fold3_atac_profile_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 13 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_acc_fold1_atac_counts_across_balanced.h5 \
--seqs npz/atac_wt_across_acc_fold1_atac_counts_seqs_across_balanced \
--scores npz/atac_wt_across_acc_fold1_atac_counts_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_acc_fold1_atac_counts
modisco motifs \
-s npz/atac_wt_across_acc_fold1_atac_counts_seqs_across_balanced.npy \
-a npz/atac_wt_across_acc_fold1_atac_counts_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_acc_fold1_atac_counts_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_acc_fold1_atac_counts_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_acc_fold1_atac_counts_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 14 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_acc_fold1_atac_profile_across_balanced.h5 \
--seqs npz/atac_wt_across_acc_fold1_atac_profile_seqs_across_balanced \
--scores npz/atac_wt_across_acc_fold1_atac_profile_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_acc_fold1_atac_profile
modisco motifs \
-s npz/atac_wt_across_acc_fold1_atac_profile_seqs_across_balanced.npy \
-a npz/atac_wt_across_acc_fold1_atac_profile_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_acc_fold1_atac_profile_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_acc_fold1_atac_profile_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_acc_fold1_atac_profile_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 15 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_acc_fold2_atac_counts_across_balanced.h5 \
--seqs npz/atac_wt_across_acc_fold2_atac_counts_seqs_across_balanced \
--scores npz/atac_wt_across_acc_fold2_atac_counts_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_acc_fold2_atac_counts
modisco motifs \
-s npz/atac_wt_across_acc_fold2_atac_counts_seqs_across_balanced.npy \
-a npz/atac_wt_across_acc_fold2_atac_counts_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_acc_fold2_atac_counts_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_acc_fold2_atac_counts_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_acc_fold2_atac_counts_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 16 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_acc_fold2_atac_profile_across_balanced.h5 \
--seqs npz/atac_wt_across_acc_fold2_atac_profile_seqs_across_balanced \
--scores npz/atac_wt_across_acc_fold2_atac_profile_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_acc_fold2_atac_profile
modisco motifs \
-s npz/atac_wt_across_acc_fold2_atac_profile_seqs_across_balanced.npy \
-a npz/atac_wt_across_acc_fold2_atac_profile_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_acc_fold2_atac_profile_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_acc_fold2_atac_profile_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_acc_fold2_atac_profile_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 17 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_acc_fold3_atac_counts_across_balanced.h5 \
--seqs npz/atac_wt_across_acc_fold3_atac_counts_seqs_across_balanced \
--scores npz/atac_wt_across_acc_fold3_atac_counts_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_acc_fold3_atac_counts
modisco motifs \
-s npz/atac_wt_across_acc_fold3_atac_counts_seqs_across_balanced.npy \
-a npz/atac_wt_across_acc_fold3_atac_counts_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_acc_fold3_atac_counts_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_acc_fold3_atac_counts_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_acc_fold3_atac_counts_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 18 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_acc_fold3_atac_profile_across_balanced.h5 \
--seqs npz/atac_wt_across_acc_fold3_atac_profile_seqs_across_balanced \
--scores npz/atac_wt_across_acc_fold3_atac_profile_scores_across_balanced
mkdir -p modiscolite/atac_wt_across_acc_fold3_atac_profile
modisco motifs \
-s npz/atac_wt_across_acc_fold3_atac_profile_seqs_across_balanced.npy \
-a npz/atac_wt_across_acc_fold3_atac_profile_scores_across_balanced.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_acc_fold3_atac_profile_across_balanced/modisco.h5
modisco report -i modiscolite/atac_wt_across_acc_fold3_atac_profile_across_balanced/modisco.h5 \
-o modiscolite/atac_wt_across_acc_fold3_atac_profile_across_balanced \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


