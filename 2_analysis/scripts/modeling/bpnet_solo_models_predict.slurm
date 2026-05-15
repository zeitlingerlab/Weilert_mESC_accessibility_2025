#!/usr/bin/bash
#SBATCH --job-name bpnet_training
#SBATCH --output=slurm_%j.log
#SBATCH --array=1-90%20
#SBATCH --mem=200gb
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=10
source /home/mw2098/.bashrc
conda deactivate
conda activate bpreveal_510
ml bedtools
cd /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/
if [[ ${SLURM_ARRAY_TASK_ID} == 1 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold1.h5                 --bw preds/bpnet_osknz_across_genome_fold1_oct4_pos.bw                 --head-id 0                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 2 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold1.h5                 --bw preds/bpnet_osknz_across_genome_fold1_oct4_neg.bw                 --head-id 0                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 3 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold1.h5                 --bw preds/bpnet_osknz_across_genome_fold1_sox2_pos.bw                 --head-id 1                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 4 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold1.h5                 --bw preds/bpnet_osknz_across_genome_fold1_sox2_neg.bw                 --head-id 1                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 5 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold1.h5                 --bw preds/bpnet_osknz_across_genome_fold1_nanog_pos.bw                 --head-id 2                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 6 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold1.h5                 --bw preds/bpnet_osknz_across_genome_fold1_nanog_neg.bw                 --head-id 2                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 7 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold1.h5                 --bw preds/bpnet_osknz_across_genome_fold1_klf4_pos.bw                 --head-id 3                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 8 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold1.h5                 --bw preds/bpnet_osknz_across_genome_fold1_klf4_neg.bw                 --head-id 3                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 9 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold1.h5                 --bw preds/bpnet_osknz_across_genome_fold1_zic3_pos.bw                 --head-id 4                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 10 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold1.h5                 --bw preds/bpnet_osknz_across_genome_fold1_zic3_neg.bw                 --head-id 4                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 11 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold2.h5                 --bw preds/bpnet_osknz_across_genome_fold2_oct4_pos.bw                 --head-id 0                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 12 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold2.h5                 --bw preds/bpnet_osknz_across_genome_fold2_oct4_neg.bw                 --head-id 0                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 13 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold2.h5                 --bw preds/bpnet_osknz_across_genome_fold2_sox2_pos.bw                 --head-id 1                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 14 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold2.h5                 --bw preds/bpnet_osknz_across_genome_fold2_sox2_neg.bw                 --head-id 1                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 15 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold2.h5                 --bw preds/bpnet_osknz_across_genome_fold2_nanog_pos.bw                 --head-id 2                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 16 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold2.h5                 --bw preds/bpnet_osknz_across_genome_fold2_nanog_neg.bw                 --head-id 2                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 17 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold2.h5                 --bw preds/bpnet_osknz_across_genome_fold2_klf4_pos.bw                 --head-id 3                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 18 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold2.h5                 --bw preds/bpnet_osknz_across_genome_fold2_klf4_neg.bw                 --head-id 3                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 19 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold2.h5                 --bw preds/bpnet_osknz_across_genome_fold2_zic3_pos.bw                 --head-id 4                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 20 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold2.h5                 --bw preds/bpnet_osknz_across_genome_fold2_zic3_neg.bw                 --head-id 4                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 21 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold3.h5                 --bw preds/bpnet_osknz_across_genome_fold3_oct4_pos.bw                 --head-id 0                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 22 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold3.h5                 --bw preds/bpnet_osknz_across_genome_fold3_oct4_neg.bw                 --head-id 0                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 23 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold3.h5                 --bw preds/bpnet_osknz_across_genome_fold3_sox2_pos.bw                 --head-id 1                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 24 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold3.h5                 --bw preds/bpnet_osknz_across_genome_fold3_sox2_neg.bw                 --head-id 1                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 25 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold3.h5                 --bw preds/bpnet_osknz_across_genome_fold3_nanog_pos.bw                 --head-id 2                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 26 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold3.h5                 --bw preds/bpnet_osknz_across_genome_fold3_nanog_neg.bw                 --head-id 2                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 27 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold3.h5                 --bw preds/bpnet_osknz_across_genome_fold3_klf4_pos.bw                 --head-id 3                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 28 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold3.h5                 --bw preds/bpnet_osknz_across_genome_fold3_klf4_neg.bw                 --head-id 3                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 29 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold3.h5                 --bw preds/bpnet_osknz_across_genome_fold3_zic3_pos.bw                 --head-id 4                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 30 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_genome_fold3.h5                 --bw preds/bpnet_osknz_across_genome_fold3_zic3_neg.bw                 --head-id 4                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 31 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold1.h5                 --bw preds/bpnet_osknz_across_split_fold1_oct4_pos.bw                 --head-id 0                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 32 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold1.h5                 --bw preds/bpnet_osknz_across_split_fold1_oct4_neg.bw                 --head-id 0                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 33 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold1.h5                 --bw preds/bpnet_osknz_across_split_fold1_sox2_pos.bw                 --head-id 1                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 34 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold1.h5                 --bw preds/bpnet_osknz_across_split_fold1_sox2_neg.bw                 --head-id 1                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 35 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold1.h5                 --bw preds/bpnet_osknz_across_split_fold1_nanog_pos.bw                 --head-id 2                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 36 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold1.h5                 --bw preds/bpnet_osknz_across_split_fold1_nanog_neg.bw                 --head-id 2                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 37 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold1.h5                 --bw preds/bpnet_osknz_across_split_fold1_klf4_pos.bw                 --head-id 3                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 38 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold1.h5                 --bw preds/bpnet_osknz_across_split_fold1_klf4_neg.bw                 --head-id 3                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 39 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold1.h5                 --bw preds/bpnet_osknz_across_split_fold1_zic3_pos.bw                 --head-id 4                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 40 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold1.h5                 --bw preds/bpnet_osknz_across_split_fold1_zic3_neg.bw                 --head-id 4                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 41 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold2.h5                 --bw preds/bpnet_osknz_across_split_fold2_oct4_pos.bw                 --head-id 0                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 42 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold2.h5                 --bw preds/bpnet_osknz_across_split_fold2_oct4_neg.bw                 --head-id 0                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 43 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold2.h5                 --bw preds/bpnet_osknz_across_split_fold2_sox2_pos.bw                 --head-id 1                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 44 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold2.h5                 --bw preds/bpnet_osknz_across_split_fold2_sox2_neg.bw                 --head-id 1                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 45 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold2.h5                 --bw preds/bpnet_osknz_across_split_fold2_nanog_pos.bw                 --head-id 2                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 46 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold2.h5                 --bw preds/bpnet_osknz_across_split_fold2_nanog_neg.bw                 --head-id 2                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 47 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold2.h5                 --bw preds/bpnet_osknz_across_split_fold2_klf4_pos.bw                 --head-id 3                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 48 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold2.h5                 --bw preds/bpnet_osknz_across_split_fold2_klf4_neg.bw                 --head-id 3                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 49 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold2.h5                 --bw preds/bpnet_osknz_across_split_fold2_zic3_pos.bw                 --head-id 4                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 50 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold2.h5                 --bw preds/bpnet_osknz_across_split_fold2_zic3_neg.bw                 --head-id 4                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 51 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold3.h5                 --bw preds/bpnet_osknz_across_split_fold3_oct4_pos.bw                 --head-id 0                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 52 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold3.h5                 --bw preds/bpnet_osknz_across_split_fold3_oct4_neg.bw                 --head-id 0                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 53 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold3.h5                 --bw preds/bpnet_osknz_across_split_fold3_sox2_pos.bw                 --head-id 1                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 54 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold3.h5                 --bw preds/bpnet_osknz_across_split_fold3_sox2_neg.bw                 --head-id 1                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 55 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold3.h5                 --bw preds/bpnet_osknz_across_split_fold3_nanog_pos.bw                 --head-id 2                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 56 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold3.h5                 --bw preds/bpnet_osknz_across_split_fold3_nanog_neg.bw                 --head-id 2                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 57 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold3.h5                 --bw preds/bpnet_osknz_across_split_fold3_klf4_pos.bw                 --head-id 3                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 58 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold3.h5                 --bw preds/bpnet_osknz_across_split_fold3_klf4_neg.bw                 --head-id 3                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 59 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold3.h5                 --bw preds/bpnet_osknz_across_split_fold3_zic3_pos.bw                 --head-id 4                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 60 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_split_fold3.h5                 --bw preds/bpnet_osknz_across_split_fold3_zic3_neg.bw                 --head-id 4                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 61 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold1.h5                 --bw preds/bpnet_osknz_across_acc_fold1_oct4_pos.bw                 --head-id 0                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 62 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold1.h5                 --bw preds/bpnet_osknz_across_acc_fold1_oct4_neg.bw                 --head-id 0                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 63 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold1.h5                 --bw preds/bpnet_osknz_across_acc_fold1_sox2_pos.bw                 --head-id 1                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 64 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold1.h5                 --bw preds/bpnet_osknz_across_acc_fold1_sox2_neg.bw                 --head-id 1                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 65 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold1.h5                 --bw preds/bpnet_osknz_across_acc_fold1_nanog_pos.bw                 --head-id 2                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 66 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold1.h5                 --bw preds/bpnet_osknz_across_acc_fold1_nanog_neg.bw                 --head-id 2                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 67 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold1.h5                 --bw preds/bpnet_osknz_across_acc_fold1_klf4_pos.bw                 --head-id 3                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 68 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold1.h5                 --bw preds/bpnet_osknz_across_acc_fold1_klf4_neg.bw                 --head-id 3                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 69 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold1.h5                 --bw preds/bpnet_osknz_across_acc_fold1_zic3_pos.bw                 --head-id 4                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 70 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold1.h5                 --bw preds/bpnet_osknz_across_acc_fold1_zic3_neg.bw                 --head-id 4                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 71 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold2.h5                 --bw preds/bpnet_osknz_across_acc_fold2_oct4_pos.bw                 --head-id 0                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 72 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold2.h5                 --bw preds/bpnet_osknz_across_acc_fold2_oct4_neg.bw                 --head-id 0                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 73 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold2.h5                 --bw preds/bpnet_osknz_across_acc_fold2_sox2_pos.bw                 --head-id 1                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 74 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold2.h5                 --bw preds/bpnet_osknz_across_acc_fold2_sox2_neg.bw                 --head-id 1                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 75 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold2.h5                 --bw preds/bpnet_osknz_across_acc_fold2_nanog_pos.bw                 --head-id 2                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 76 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold2.h5                 --bw preds/bpnet_osknz_across_acc_fold2_nanog_neg.bw                 --head-id 2                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 77 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold2.h5                 --bw preds/bpnet_osknz_across_acc_fold2_klf4_pos.bw                 --head-id 3                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 78 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold2.h5                 --bw preds/bpnet_osknz_across_acc_fold2_klf4_neg.bw                 --head-id 3                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 79 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold2.h5                 --bw preds/bpnet_osknz_across_acc_fold2_zic3_pos.bw                 --head-id 4                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 80 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold2.h5                 --bw preds/bpnet_osknz_across_acc_fold2_zic3_neg.bw                 --head-id 4                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 81 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold3.h5                 --bw preds/bpnet_osknz_across_acc_fold3_oct4_pos.bw                 --head-id 0                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 82 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold3.h5                 --bw preds/bpnet_osknz_across_acc_fold3_oct4_neg.bw                 --head-id 0                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 83 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold3.h5                 --bw preds/bpnet_osknz_across_acc_fold3_sox2_pos.bw                 --head-id 1                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 84 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold3.h5                 --bw preds/bpnet_osknz_across_acc_fold3_sox2_neg.bw                 --head-id 1                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 85 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold3.h5                 --bw preds/bpnet_osknz_across_acc_fold3_nanog_pos.bw                 --head-id 2                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 86 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold3.h5                 --bw preds/bpnet_osknz_across_acc_fold3_nanog_neg.bw                 --head-id 2                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 87 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold3.h5                 --bw preds/bpnet_osknz_across_acc_fold3_klf4_pos.bw                 --head-id 3                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 88 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold3.h5                 --bw preds/bpnet_osknz_across_acc_fold3_klf4_neg.bw                 --head-id 3                 --task-id 1 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 89 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold3.h5                 --bw preds/bpnet_osknz_across_acc_fold3_zic3_pos.bw                 --head-id 4                 --task-id 0 --mode profile --verbose

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 90 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/predictToBigwig.py                 --h5 preds/bpnet_osknz_across_acc_fold3.h5                 --bw preds/bpnet_osknz_across_acc_fold3_zic3_neg.bw                 --head-id 4                 --task-id 1 --mode profile --verbose

fi


