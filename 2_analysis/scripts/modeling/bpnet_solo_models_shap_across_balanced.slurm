#!/usr/bin/bash
#SBATCH --job-name bpnet_training
#SBATCH --output=slurm_%j.log
#SBATCH --array=1-45%3
#SBATCH --mem=200gb
#SBATCH --time=72:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
source /home/mw2098/.bashrc
conda deactivate
conda activate bpreveal_510
ml bedtools
cd /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/
if [[ ${SLURM_ARRAY_TASK_ID} == 1 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_genome_fold1_oct4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 2 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_genome_fold1_sox2_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 3 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_genome_fold1_nanog_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 4 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_genome_fold1_klf4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 5 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_genome_fold1_zic3_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 6 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_genome_fold2_oct4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 7 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_genome_fold2_sox2_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 8 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_genome_fold2_nanog_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 9 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_genome_fold2_klf4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 10 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_genome_fold2_zic3_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 11 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_genome_fold3_oct4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 12 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_genome_fold3_sox2_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 13 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_genome_fold3_nanog_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 14 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_genome_fold3_klf4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 15 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_genome_fold3_zic3_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 16 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_split_fold1_oct4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 17 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_split_fold1_sox2_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 18 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_split_fold1_nanog_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 19 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_split_fold1_klf4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 20 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_split_fold1_zic3_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 21 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_split_fold2_oct4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 22 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_split_fold2_sox2_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 23 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_split_fold2_nanog_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 24 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_split_fold2_klf4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 25 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_split_fold2_zic3_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 26 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_split_fold3_oct4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 27 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_split_fold3_sox2_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 28 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_split_fold3_nanog_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 29 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_split_fold3_klf4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 30 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_split_fold3_zic3_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 31 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_acc_fold1_oct4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 32 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_acc_fold1_sox2_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 33 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_acc_fold1_nanog_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 34 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_acc_fold1_klf4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 35 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_acc_fold1_zic3_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 36 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_acc_fold2_oct4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 37 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_acc_fold2_sox2_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 38 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_acc_fold2_nanog_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 39 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_acc_fold2_klf4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 40 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_acc_fold2_zic3_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 41 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_acc_fold3_oct4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 42 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_acc_fold3_sox2_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 43 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_acc_fold3_nanog_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 44 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_acc_fold3_klf4_across_balanced.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 45 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/interpretFlat.py json/shapFlat_bpnet_osknz_across_acc_fold3_zic3_across_balanced.json

fi


