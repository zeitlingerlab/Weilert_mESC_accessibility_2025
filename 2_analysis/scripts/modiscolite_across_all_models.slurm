#!/usr/bin/bash
#SBATCH --array=1-36%15
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
--h5 shap/bunina_atac_day0_across_combined_fold1_atac_counts.h5 \
--seqs npz/bunina_atac_day0_across_combined_fold1_atac_counts_seqs \
--scores npz/bunina_atac_day0_across_combined_fold1_atac_counts_scores
mkdir -p modiscolite/bunina_atac_day0_across_combined_fold1_atac_counts
modisco motifs \
-s npz/bunina_atac_day0_across_combined_fold1_atac_counts_seqs.npy \
-a npz/bunina_atac_day0_across_combined_fold1_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day0_across_combined_fold1_atac_counts/modisco.h5
modisco report -i modiscolite/bunina_atac_day0_across_combined_fold1_atac_counts/modisco.h5 \
-o modiscolite/bunina_atac_day0_across_combined_fold1_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 2 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day0_across_combined_fold1_atac_profile.h5 \
--seqs npz/bunina_atac_day0_across_combined_fold1_atac_profile_seqs \
--scores npz/bunina_atac_day0_across_combined_fold1_atac_profile_scores
mkdir -p modiscolite/bunina_atac_day0_across_combined_fold1_atac_profile
modisco motifs \
-s npz/bunina_atac_day0_across_combined_fold1_atac_profile_seqs.npy \
-a npz/bunina_atac_day0_across_combined_fold1_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day0_across_combined_fold1_atac_profile/modisco.h5
modisco report -i modiscolite/bunina_atac_day0_across_combined_fold1_atac_profile/modisco.h5 \
-o modiscolite/bunina_atac_day0_across_combined_fold1_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 3 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day0_across_combined_fold2_atac_counts.h5 \
--seqs npz/bunina_atac_day0_across_combined_fold2_atac_counts_seqs \
--scores npz/bunina_atac_day0_across_combined_fold2_atac_counts_scores
mkdir -p modiscolite/bunina_atac_day0_across_combined_fold2_atac_counts
modisco motifs \
-s npz/bunina_atac_day0_across_combined_fold2_atac_counts_seqs.npy \
-a npz/bunina_atac_day0_across_combined_fold2_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day0_across_combined_fold2_atac_counts/modisco.h5
modisco report -i modiscolite/bunina_atac_day0_across_combined_fold2_atac_counts/modisco.h5 \
-o modiscolite/bunina_atac_day0_across_combined_fold2_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 4 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day0_across_combined_fold2_atac_profile.h5 \
--seqs npz/bunina_atac_day0_across_combined_fold2_atac_profile_seqs \
--scores npz/bunina_atac_day0_across_combined_fold2_atac_profile_scores
mkdir -p modiscolite/bunina_atac_day0_across_combined_fold2_atac_profile
modisco motifs \
-s npz/bunina_atac_day0_across_combined_fold2_atac_profile_seqs.npy \
-a npz/bunina_atac_day0_across_combined_fold2_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day0_across_combined_fold2_atac_profile/modisco.h5
modisco report -i modiscolite/bunina_atac_day0_across_combined_fold2_atac_profile/modisco.h5 \
-o modiscolite/bunina_atac_day0_across_combined_fold2_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 5 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day0_across_combined_fold3_atac_counts.h5 \
--seqs npz/bunina_atac_day0_across_combined_fold3_atac_counts_seqs \
--scores npz/bunina_atac_day0_across_combined_fold3_atac_counts_scores
mkdir -p modiscolite/bunina_atac_day0_across_combined_fold3_atac_counts
modisco motifs \
-s npz/bunina_atac_day0_across_combined_fold3_atac_counts_seqs.npy \
-a npz/bunina_atac_day0_across_combined_fold3_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day0_across_combined_fold3_atac_counts/modisco.h5
modisco report -i modiscolite/bunina_atac_day0_across_combined_fold3_atac_counts/modisco.h5 \
-o modiscolite/bunina_atac_day0_across_combined_fold3_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 6 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day0_across_combined_fold3_atac_profile.h5 \
--seqs npz/bunina_atac_day0_across_combined_fold3_atac_profile_seqs \
--scores npz/bunina_atac_day0_across_combined_fold3_atac_profile_scores
mkdir -p modiscolite/bunina_atac_day0_across_combined_fold3_atac_profile
modisco motifs \
-s npz/bunina_atac_day0_across_combined_fold3_atac_profile_seqs.npy \
-a npz/bunina_atac_day0_across_combined_fold3_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day0_across_combined_fold3_atac_profile/modisco.h5
modisco report -i modiscolite/bunina_atac_day0_across_combined_fold3_atac_profile/modisco.h5 \
-o modiscolite/bunina_atac_day0_across_combined_fold3_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 7 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day4_across_combined_fold1_atac_counts.h5 \
--seqs npz/bunina_atac_day4_across_combined_fold1_atac_counts_seqs \
--scores npz/bunina_atac_day4_across_combined_fold1_atac_counts_scores
mkdir -p modiscolite/bunina_atac_day4_across_combined_fold1_atac_counts
modisco motifs \
-s npz/bunina_atac_day4_across_combined_fold1_atac_counts_seqs.npy \
-a npz/bunina_atac_day4_across_combined_fold1_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day4_across_combined_fold1_atac_counts/modisco.h5
modisco report -i modiscolite/bunina_atac_day4_across_combined_fold1_atac_counts/modisco.h5 \
-o modiscolite/bunina_atac_day4_across_combined_fold1_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 8 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day4_across_combined_fold1_atac_profile.h5 \
--seqs npz/bunina_atac_day4_across_combined_fold1_atac_profile_seqs \
--scores npz/bunina_atac_day4_across_combined_fold1_atac_profile_scores
mkdir -p modiscolite/bunina_atac_day4_across_combined_fold1_atac_profile
modisco motifs \
-s npz/bunina_atac_day4_across_combined_fold1_atac_profile_seqs.npy \
-a npz/bunina_atac_day4_across_combined_fold1_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day4_across_combined_fold1_atac_profile/modisco.h5
modisco report -i modiscolite/bunina_atac_day4_across_combined_fold1_atac_profile/modisco.h5 \
-o modiscolite/bunina_atac_day4_across_combined_fold1_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 9 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day4_across_combined_fold2_atac_counts.h5 \
--seqs npz/bunina_atac_day4_across_combined_fold2_atac_counts_seqs \
--scores npz/bunina_atac_day4_across_combined_fold2_atac_counts_scores
mkdir -p modiscolite/bunina_atac_day4_across_combined_fold2_atac_counts
modisco motifs \
-s npz/bunina_atac_day4_across_combined_fold2_atac_counts_seqs.npy \
-a npz/bunina_atac_day4_across_combined_fold2_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day4_across_combined_fold2_atac_counts/modisco.h5
modisco report -i modiscolite/bunina_atac_day4_across_combined_fold2_atac_counts/modisco.h5 \
-o modiscolite/bunina_atac_day4_across_combined_fold2_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 10 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day4_across_combined_fold2_atac_profile.h5 \
--seqs npz/bunina_atac_day4_across_combined_fold2_atac_profile_seqs \
--scores npz/bunina_atac_day4_across_combined_fold2_atac_profile_scores
mkdir -p modiscolite/bunina_atac_day4_across_combined_fold2_atac_profile
modisco motifs \
-s npz/bunina_atac_day4_across_combined_fold2_atac_profile_seqs.npy \
-a npz/bunina_atac_day4_across_combined_fold2_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day4_across_combined_fold2_atac_profile/modisco.h5
modisco report -i modiscolite/bunina_atac_day4_across_combined_fold2_atac_profile/modisco.h5 \
-o modiscolite/bunina_atac_day4_across_combined_fold2_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 11 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day4_across_combined_fold3_atac_counts.h5 \
--seqs npz/bunina_atac_day4_across_combined_fold3_atac_counts_seqs \
--scores npz/bunina_atac_day4_across_combined_fold3_atac_counts_scores
mkdir -p modiscolite/bunina_atac_day4_across_combined_fold3_atac_counts
modisco motifs \
-s npz/bunina_atac_day4_across_combined_fold3_atac_counts_seqs.npy \
-a npz/bunina_atac_day4_across_combined_fold3_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day4_across_combined_fold3_atac_counts/modisco.h5
modisco report -i modiscolite/bunina_atac_day4_across_combined_fold3_atac_counts/modisco.h5 \
-o modiscolite/bunina_atac_day4_across_combined_fold3_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 12 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day4_across_combined_fold3_atac_profile.h5 \
--seqs npz/bunina_atac_day4_across_combined_fold3_atac_profile_seqs \
--scores npz/bunina_atac_day4_across_combined_fold3_atac_profile_scores
mkdir -p modiscolite/bunina_atac_day4_across_combined_fold3_atac_profile
modisco motifs \
-s npz/bunina_atac_day4_across_combined_fold3_atac_profile_seqs.npy \
-a npz/bunina_atac_day4_across_combined_fold3_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day4_across_combined_fold3_atac_profile/modisco.h5
modisco report -i modiscolite/bunina_atac_day4_across_combined_fold3_atac_profile/modisco.h5 \
-o modiscolite/bunina_atac_day4_across_combined_fold3_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 13 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day8_across_combined_fold1_atac_counts.h5 \
--seqs npz/bunina_atac_day8_across_combined_fold1_atac_counts_seqs \
--scores npz/bunina_atac_day8_across_combined_fold1_atac_counts_scores
mkdir -p modiscolite/bunina_atac_day8_across_combined_fold1_atac_counts
modisco motifs \
-s npz/bunina_atac_day8_across_combined_fold1_atac_counts_seqs.npy \
-a npz/bunina_atac_day8_across_combined_fold1_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day8_across_combined_fold1_atac_counts/modisco.h5
modisco report -i modiscolite/bunina_atac_day8_across_combined_fold1_atac_counts/modisco.h5 \
-o modiscolite/bunina_atac_day8_across_combined_fold1_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 14 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day8_across_combined_fold1_atac_profile.h5 \
--seqs npz/bunina_atac_day8_across_combined_fold1_atac_profile_seqs \
--scores npz/bunina_atac_day8_across_combined_fold1_atac_profile_scores
mkdir -p modiscolite/bunina_atac_day8_across_combined_fold1_atac_profile
modisco motifs \
-s npz/bunina_atac_day8_across_combined_fold1_atac_profile_seqs.npy \
-a npz/bunina_atac_day8_across_combined_fold1_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day8_across_combined_fold1_atac_profile/modisco.h5
modisco report -i modiscolite/bunina_atac_day8_across_combined_fold1_atac_profile/modisco.h5 \
-o modiscolite/bunina_atac_day8_across_combined_fold1_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 15 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day8_across_combined_fold2_atac_counts.h5 \
--seqs npz/bunina_atac_day8_across_combined_fold2_atac_counts_seqs \
--scores npz/bunina_atac_day8_across_combined_fold2_atac_counts_scores
mkdir -p modiscolite/bunina_atac_day8_across_combined_fold2_atac_counts
modisco motifs \
-s npz/bunina_atac_day8_across_combined_fold2_atac_counts_seqs.npy \
-a npz/bunina_atac_day8_across_combined_fold2_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day8_across_combined_fold2_atac_counts/modisco.h5
modisco report -i modiscolite/bunina_atac_day8_across_combined_fold2_atac_counts/modisco.h5 \
-o modiscolite/bunina_atac_day8_across_combined_fold2_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 16 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day8_across_combined_fold2_atac_profile.h5 \
--seqs npz/bunina_atac_day8_across_combined_fold2_atac_profile_seqs \
--scores npz/bunina_atac_day8_across_combined_fold2_atac_profile_scores
mkdir -p modiscolite/bunina_atac_day8_across_combined_fold2_atac_profile
modisco motifs \
-s npz/bunina_atac_day8_across_combined_fold2_atac_profile_seqs.npy \
-a npz/bunina_atac_day8_across_combined_fold2_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day8_across_combined_fold2_atac_profile/modisco.h5
modisco report -i modiscolite/bunina_atac_day8_across_combined_fold2_atac_profile/modisco.h5 \
-o modiscolite/bunina_atac_day8_across_combined_fold2_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 17 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day8_across_combined_fold3_atac_counts.h5 \
--seqs npz/bunina_atac_day8_across_combined_fold3_atac_counts_seqs \
--scores npz/bunina_atac_day8_across_combined_fold3_atac_counts_scores
mkdir -p modiscolite/bunina_atac_day8_across_combined_fold3_atac_counts
modisco motifs \
-s npz/bunina_atac_day8_across_combined_fold3_atac_counts_seqs.npy \
-a npz/bunina_atac_day8_across_combined_fold3_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day8_across_combined_fold3_atac_counts/modisco.h5
modisco report -i modiscolite/bunina_atac_day8_across_combined_fold3_atac_counts/modisco.h5 \
-o modiscolite/bunina_atac_day8_across_combined_fold3_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 18 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day8_across_combined_fold3_atac_profile.h5 \
--seqs npz/bunina_atac_day8_across_combined_fold3_atac_profile_seqs \
--scores npz/bunina_atac_day8_across_combined_fold3_atac_profile_scores
mkdir -p modiscolite/bunina_atac_day8_across_combined_fold3_atac_profile
modisco motifs \
-s npz/bunina_atac_day8_across_combined_fold3_atac_profile_seqs.npy \
-a npz/bunina_atac_day8_across_combined_fold3_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day8_across_combined_fold3_atac_profile/modisco.h5
modisco report -i modiscolite/bunina_atac_day8_across_combined_fold3_atac_profile/modisco.h5 \
-o modiscolite/bunina_atac_day8_across_combined_fold3_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 19 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day12_across_combined_fold1_atac_counts.h5 \
--seqs npz/bunina_atac_day12_across_combined_fold1_atac_counts_seqs \
--scores npz/bunina_atac_day12_across_combined_fold1_atac_counts_scores
mkdir -p modiscolite/bunina_atac_day12_across_combined_fold1_atac_counts
modisco motifs \
-s npz/bunina_atac_day12_across_combined_fold1_atac_counts_seqs.npy \
-a npz/bunina_atac_day12_across_combined_fold1_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day12_across_combined_fold1_atac_counts/modisco.h5
modisco report -i modiscolite/bunina_atac_day12_across_combined_fold1_atac_counts/modisco.h5 \
-o modiscolite/bunina_atac_day12_across_combined_fold1_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 20 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day12_across_combined_fold1_atac_profile.h5 \
--seqs npz/bunina_atac_day12_across_combined_fold1_atac_profile_seqs \
--scores npz/bunina_atac_day12_across_combined_fold1_atac_profile_scores
mkdir -p modiscolite/bunina_atac_day12_across_combined_fold1_atac_profile
modisco motifs \
-s npz/bunina_atac_day12_across_combined_fold1_atac_profile_seqs.npy \
-a npz/bunina_atac_day12_across_combined_fold1_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day12_across_combined_fold1_atac_profile/modisco.h5
modisco report -i modiscolite/bunina_atac_day12_across_combined_fold1_atac_profile/modisco.h5 \
-o modiscolite/bunina_atac_day12_across_combined_fold1_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 21 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day12_across_combined_fold2_atac_counts.h5 \
--seqs npz/bunina_atac_day12_across_combined_fold2_atac_counts_seqs \
--scores npz/bunina_atac_day12_across_combined_fold2_atac_counts_scores
mkdir -p modiscolite/bunina_atac_day12_across_combined_fold2_atac_counts
modisco motifs \
-s npz/bunina_atac_day12_across_combined_fold2_atac_counts_seqs.npy \
-a npz/bunina_atac_day12_across_combined_fold2_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day12_across_combined_fold2_atac_counts/modisco.h5
modisco report -i modiscolite/bunina_atac_day12_across_combined_fold2_atac_counts/modisco.h5 \
-o modiscolite/bunina_atac_day12_across_combined_fold2_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 22 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day12_across_combined_fold2_atac_profile.h5 \
--seqs npz/bunina_atac_day12_across_combined_fold2_atac_profile_seqs \
--scores npz/bunina_atac_day12_across_combined_fold2_atac_profile_scores
mkdir -p modiscolite/bunina_atac_day12_across_combined_fold2_atac_profile
modisco motifs \
-s npz/bunina_atac_day12_across_combined_fold2_atac_profile_seqs.npy \
-a npz/bunina_atac_day12_across_combined_fold2_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day12_across_combined_fold2_atac_profile/modisco.h5
modisco report -i modiscolite/bunina_atac_day12_across_combined_fold2_atac_profile/modisco.h5 \
-o modiscolite/bunina_atac_day12_across_combined_fold2_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 23 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day12_across_combined_fold3_atac_counts.h5 \
--seqs npz/bunina_atac_day12_across_combined_fold3_atac_counts_seqs \
--scores npz/bunina_atac_day12_across_combined_fold3_atac_counts_scores
mkdir -p modiscolite/bunina_atac_day12_across_combined_fold3_atac_counts
modisco motifs \
-s npz/bunina_atac_day12_across_combined_fold3_atac_counts_seqs.npy \
-a npz/bunina_atac_day12_across_combined_fold3_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day12_across_combined_fold3_atac_counts/modisco.h5
modisco report -i modiscolite/bunina_atac_day12_across_combined_fold3_atac_counts/modisco.h5 \
-o modiscolite/bunina_atac_day12_across_combined_fold3_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 24 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/bunina_atac_day12_across_combined_fold3_atac_profile.h5 \
--seqs npz/bunina_atac_day12_across_combined_fold3_atac_profile_seqs \
--scores npz/bunina_atac_day12_across_combined_fold3_atac_profile_scores
mkdir -p modiscolite/bunina_atac_day12_across_combined_fold3_atac_profile
modisco motifs \
-s npz/bunina_atac_day12_across_combined_fold3_atac_profile_seqs.npy \
-a npz/bunina_atac_day12_across_combined_fold3_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/bunina_atac_day12_across_combined_fold3_atac_profile/modisco.h5
modisco report -i modiscolite/bunina_atac_day12_across_combined_fold3_atac_profile/modisco.h5 \
-o modiscolite/bunina_atac_day12_across_combined_fold3_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 25 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_split_fold1_atac_counts.h5 \
--seqs npz/atac_wt_across_split_fold1_atac_counts_seqs \
--scores npz/atac_wt_across_split_fold1_atac_counts_scores
mkdir -p modiscolite/atac_wt_across_split_fold1_atac_counts
modisco motifs \
-s npz/atac_wt_across_split_fold1_atac_counts_seqs.npy \
-a npz/atac_wt_across_split_fold1_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_split_fold1_atac_counts/modisco.h5
modisco report -i modiscolite/atac_wt_across_split_fold1_atac_counts/modisco.h5 \
-o modiscolite/atac_wt_across_split_fold1_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 26 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_split_fold1_atac_profile.h5 \
--seqs npz/atac_wt_across_split_fold1_atac_profile_seqs \
--scores npz/atac_wt_across_split_fold1_atac_profile_scores
mkdir -p modiscolite/atac_wt_across_split_fold1_atac_profile
modisco motifs \
-s npz/atac_wt_across_split_fold1_atac_profile_seqs.npy \
-a npz/atac_wt_across_split_fold1_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_split_fold1_atac_profile/modisco.h5
modisco report -i modiscolite/atac_wt_across_split_fold1_atac_profile/modisco.h5 \
-o modiscolite/atac_wt_across_split_fold1_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 27 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_split_fold2_atac_counts.h5 \
--seqs npz/atac_wt_across_split_fold2_atac_counts_seqs \
--scores npz/atac_wt_across_split_fold2_atac_counts_scores
mkdir -p modiscolite/atac_wt_across_split_fold2_atac_counts
modisco motifs \
-s npz/atac_wt_across_split_fold2_atac_counts_seqs.npy \
-a npz/atac_wt_across_split_fold2_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_split_fold2_atac_counts/modisco.h5
modisco report -i modiscolite/atac_wt_across_split_fold2_atac_counts/modisco.h5 \
-o modiscolite/atac_wt_across_split_fold2_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 28 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_split_fold2_atac_profile.h5 \
--seqs npz/atac_wt_across_split_fold2_atac_profile_seqs \
--scores npz/atac_wt_across_split_fold2_atac_profile_scores
mkdir -p modiscolite/atac_wt_across_split_fold2_atac_profile
modisco motifs \
-s npz/atac_wt_across_split_fold2_atac_profile_seqs.npy \
-a npz/atac_wt_across_split_fold2_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_split_fold2_atac_profile/modisco.h5
modisco report -i modiscolite/atac_wt_across_split_fold2_atac_profile/modisco.h5 \
-o modiscolite/atac_wt_across_split_fold2_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 29 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_split_fold3_atac_counts.h5 \
--seqs npz/atac_wt_across_split_fold3_atac_counts_seqs \
--scores npz/atac_wt_across_split_fold3_atac_counts_scores
mkdir -p modiscolite/atac_wt_across_split_fold3_atac_counts
modisco motifs \
-s npz/atac_wt_across_split_fold3_atac_counts_seqs.npy \
-a npz/atac_wt_across_split_fold3_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_split_fold3_atac_counts/modisco.h5
modisco report -i modiscolite/atac_wt_across_split_fold3_atac_counts/modisco.h5 \
-o modiscolite/atac_wt_across_split_fold3_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 30 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_split_fold3_atac_profile.h5 \
--seqs npz/atac_wt_across_split_fold3_atac_profile_seqs \
--scores npz/atac_wt_across_split_fold3_atac_profile_scores
mkdir -p modiscolite/atac_wt_across_split_fold3_atac_profile
modisco motifs \
-s npz/atac_wt_across_split_fold3_atac_profile_seqs.npy \
-a npz/atac_wt_across_split_fold3_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_split_fold3_atac_profile/modisco.h5
modisco report -i modiscolite/atac_wt_across_split_fold3_atac_profile/modisco.h5 \
-o modiscolite/atac_wt_across_split_fold3_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 31 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_acc_fold1_atac_counts.h5 \
--seqs npz/atac_wt_across_acc_fold1_atac_counts_seqs \
--scores npz/atac_wt_across_acc_fold1_atac_counts_scores
mkdir -p modiscolite/atac_wt_across_acc_fold1_atac_counts
modisco motifs \
-s npz/atac_wt_across_acc_fold1_atac_counts_seqs.npy \
-a npz/atac_wt_across_acc_fold1_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_acc_fold1_atac_counts/modisco.h5
modisco report -i modiscolite/atac_wt_across_acc_fold1_atac_counts/modisco.h5 \
-o modiscolite/atac_wt_across_acc_fold1_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 32 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_acc_fold1_atac_profile.h5 \
--seqs npz/atac_wt_across_acc_fold1_atac_profile_seqs \
--scores npz/atac_wt_across_acc_fold1_atac_profile_scores
mkdir -p modiscolite/atac_wt_across_acc_fold1_atac_profile
modisco motifs \
-s npz/atac_wt_across_acc_fold1_atac_profile_seqs.npy \
-a npz/atac_wt_across_acc_fold1_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_acc_fold1_atac_profile/modisco.h5
modisco report -i modiscolite/atac_wt_across_acc_fold1_atac_profile/modisco.h5 \
-o modiscolite/atac_wt_across_acc_fold1_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 33 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_acc_fold2_atac_counts.h5 \
--seqs npz/atac_wt_across_acc_fold2_atac_counts_seqs \
--scores npz/atac_wt_across_acc_fold2_atac_counts_scores
mkdir -p modiscolite/atac_wt_across_acc_fold2_atac_counts
modisco motifs \
-s npz/atac_wt_across_acc_fold2_atac_counts_seqs.npy \
-a npz/atac_wt_across_acc_fold2_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_acc_fold2_atac_counts/modisco.h5
modisco report -i modiscolite/atac_wt_across_acc_fold2_atac_counts/modisco.h5 \
-o modiscolite/atac_wt_across_acc_fold2_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 34 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_acc_fold2_atac_profile.h5 \
--seqs npz/atac_wt_across_acc_fold2_atac_profile_seqs \
--scores npz/atac_wt_across_acc_fold2_atac_profile_scores
mkdir -p modiscolite/atac_wt_across_acc_fold2_atac_profile
modisco motifs \
-s npz/atac_wt_across_acc_fold2_atac_profile_seqs.npy \
-a npz/atac_wt_across_acc_fold2_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_acc_fold2_atac_profile/modisco.h5
modisco report -i modiscolite/atac_wt_across_acc_fold2_atac_profile/modisco.h5 \
-o modiscolite/atac_wt_across_acc_fold2_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 35 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_acc_fold3_atac_counts.h5 \
--seqs npz/atac_wt_across_acc_fold3_atac_counts_seqs \
--scores npz/atac_wt_across_acc_fold3_atac_counts_scores
mkdir -p modiscolite/atac_wt_across_acc_fold3_atac_counts
modisco motifs \
-s npz/atac_wt_across_acc_fold3_atac_counts_seqs.npy \
-a npz/atac_wt_across_acc_fold3_atac_counts_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_acc_fold3_atac_counts/modisco.h5
modisco report -i modiscolite/atac_wt_across_acc_fold3_atac_counts/modisco.h5 \
-o modiscolite/atac_wt_across_acc_fold3_atac_counts \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 36 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/shapToNumpy.py \
--h5 shap/atac_wt_across_acc_fold3_atac_profile.h5 \
--seqs npz/atac_wt_across_acc_fold3_atac_profile_seqs \
--scores npz/atac_wt_across_acc_fold3_atac_profile_scores
mkdir -p modiscolite/atac_wt_across_acc_fold3_atac_profile
modisco motifs \
-s npz/atac_wt_across_acc_fold3_atac_profile_seqs.npy \
-a npz/atac_wt_across_acc_fold3_atac_profile_scores.npy \
-n 10000 \
-w 1000 \
--output modiscolite/atac_wt_across_acc_fold3_atac_profile/modisco.h5
modisco report -i modiscolite/atac_wt_across_acc_fold3_atac_profile/modisco.h5 \
-o modiscolite/atac_wt_across_acc_fold3_atac_profile \
-m /n/projects/mw2098/publications/2024_weilert_acc/public/databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt

fi


