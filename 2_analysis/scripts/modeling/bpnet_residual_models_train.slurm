#!/usr/bin/bash
#SBATCH --job-name bpnet_training
#SBATCH --output=slurm_%j.log
#SBATCH --array=1-9%3
#SBATCH --mem=500gb
#SBATCH --time=72:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
source /home/mw2098/.bashrc
conda deactivate
conda activate bpreveal_404
ml bedtools
cd /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/
if [[ ${SLURM_ARRAY_TASK_ID} == 1 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_genome_fold1_residual.json
/home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_genome_fold1_combined.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 2 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_genome_fold2_residual.json
/home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_genome_fold2_combined.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 3 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_genome_fold3_residual.json
/home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_genome_fold3_combined.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 4 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_split_fold1_residual.json
/home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_split_fold1_combined.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 5 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_split_fold2_residual.json
/home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_split_fold2_combined.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 6 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_split_fold3_residual.json
/home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_split_fold3_combined.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 7 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_acc_fold1_residual.json
/home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_acc_fold1_combined.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 8 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_acc_fold2_residual.json
/home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_acc_fold2_combined.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 9 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_acc_fold3_residual.json
/home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/predictPeaks_atac_wt_across_acc_fold3_combined.json

fi


