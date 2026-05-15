#!/usr/bin/bash
#SBATCH --job-name bpnet_training
#SBATCH --output=slurm_%j.log
#SBATCH --array=1-18%20
#SBATCH --mem=200gb
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=10
source /home/mw2098/.bashrc
conda deactivate
conda activate bpreveal_404
ml bedtools
cd /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/
if [[ ${SLURM_ARRAY_TASK_ID} == 1 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_genome_fold1_residual.h5                     --bw preds/atac_wt_across_genome_fold1_atac_cutsite_residual.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 2 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_genome_fold1_combined.h5                     --bw preds/atac_wt_across_genome_fold1_atac_cutsite_combined.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 3 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_genome_fold2_residual.h5                     --bw preds/atac_wt_across_genome_fold2_atac_cutsite_residual.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 4 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_genome_fold2_combined.h5                     --bw preds/atac_wt_across_genome_fold2_atac_cutsite_combined.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 5 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_genome_fold3_residual.h5                     --bw preds/atac_wt_across_genome_fold3_atac_cutsite_residual.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 6 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_genome_fold3_combined.h5                     --bw preds/atac_wt_across_genome_fold3_atac_cutsite_combined.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 7 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_split_fold1_residual.h5                     --bw preds/atac_wt_across_split_fold1_atac_cutsite_residual.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 8 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_split_fold1_combined.h5                     --bw preds/atac_wt_across_split_fold1_atac_cutsite_combined.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 9 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_split_fold2_residual.h5                     --bw preds/atac_wt_across_split_fold2_atac_cutsite_residual.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 10 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_split_fold2_combined.h5                     --bw preds/atac_wt_across_split_fold2_atac_cutsite_combined.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 11 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_split_fold3_residual.h5                     --bw preds/atac_wt_across_split_fold3_atac_cutsite_residual.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 12 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_split_fold3_combined.h5                     --bw preds/atac_wt_across_split_fold3_atac_cutsite_combined.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 13 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_acc_fold1_residual.h5                     --bw preds/atac_wt_across_acc_fold1_atac_cutsite_residual.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 14 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_acc_fold1_combined.h5                     --bw preds/atac_wt_across_acc_fold1_atac_cutsite_combined.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 15 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_acc_fold2_residual.h5                     --bw preds/atac_wt_across_acc_fold2_atac_cutsite_residual.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 16 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_acc_fold2_combined.h5                     --bw preds/atac_wt_across_acc_fold2_atac_cutsite_combined.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 17 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_acc_fold3_residual.h5                     --bw preds/atac_wt_across_acc_fold3_atac_cutsite_residual.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 18 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                     --h5 preds/atac_wt_across_acc_fold3_combined.h5                     --bw preds/atac_wt_across_acc_fold3_atac_cutsite_combined.bw                     --head-id 0                     --task-id 0 --mode profile --verbose

fi


