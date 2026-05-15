#!/usr/bin/bash
#SBATCH --job-name bpnet_training
#SBATCH --output=slurm_%j.log
#SBATCH --array=1-6%3
#SBATCH --mem=200gb
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

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_atac_wt_across_split_fold1_atac.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 2 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_atac_wt_across_split_fold2_atac.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 3 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_atac_wt_across_split_fold3_atac.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 4 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_atac_wt_across_acc_fold1_atac.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 5 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_atac_wt_across_acc_fold2_atac.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 6 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_atac_wt_across_acc_fold3_atac.json

fi


