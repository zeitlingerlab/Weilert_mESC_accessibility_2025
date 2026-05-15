#!/usr/bin/bash
#SBATCH --job-name bpnet_training
#SBATCH --output=slurm_%j.log
#SBATCH --array=1-54%20
#SBATCH --mem=200gb
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=10
source /home/mw2098/.bashrc
conda deactivate
conda activate bpreveal_404
ml bedtools
cd /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/
if [[ ${SLURM_ARRAY_TASK_ID} == 1 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold1_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold1_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold1_atac_cutsite_residual_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 2 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold1_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold1_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold1_atac_cutsite_combined_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 3 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold1_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold1_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold1_atac_cutsite_residual_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 4 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold1_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold1_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold1_atac_cutsite_combined_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 5 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold1_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold1_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold1_atac_cutsite_residual_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 6 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold1_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold1_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold1_atac_cutsite_combined_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 7 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold2_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold2_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold2_atac_cutsite_residual_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 8 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold2_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold2_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold2_atac_cutsite_combined_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 9 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold2_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold2_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold2_atac_cutsite_residual_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 10 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold2_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold2_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold2_atac_cutsite_combined_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 11 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold2_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold2_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold2_atac_cutsite_residual_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 12 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold2_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold2_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold2_atac_cutsite_combined_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 13 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold3_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold3_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold3_atac_cutsite_residual_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 14 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold3_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold3_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold3_atac_cutsite_combined_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 15 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold3_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold3_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold3_atac_cutsite_residual_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 16 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold3_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold3_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold3_atac_cutsite_combined_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 17 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold3_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold3_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold3_atac_cutsite_residual_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 18 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_genome_fold3_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_genome_fold3_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_genome_fold3_atac_cutsite_combined_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 19 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold1_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_split_fold1_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold1_atac_cutsite_residual_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 20 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold1_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_split_fold1_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold1_atac_cutsite_combined_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 21 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold1_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_split_fold1_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold1_atac_cutsite_residual_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 22 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold1_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_split_fold1_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold1_atac_cutsite_combined_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 23 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold1_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_split_fold1_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold1_atac_cutsite_residual_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 24 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold1_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_split_fold1_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold1_atac_cutsite_combined_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 25 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold2_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_split_fold2_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold2_atac_cutsite_residual_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 26 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold2_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_split_fold2_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold2_atac_cutsite_combined_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 27 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold2_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_split_fold2_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold2_atac_cutsite_residual_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 28 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold2_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_split_fold2_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold2_atac_cutsite_combined_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 29 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold2_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_split_fold2_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold2_atac_cutsite_residual_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 30 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold2_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_split_fold2_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold2_atac_cutsite_combined_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 31 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold3_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_split_fold3_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold3_atac_cutsite_residual_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 32 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold3_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_split_fold3_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold3_atac_cutsite_combined_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 33 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold3_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_split_fold3_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold3_atac_cutsite_residual_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 34 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold3_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_split_fold3_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold3_atac_cutsite_combined_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 35 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold3_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_split_fold3_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold3_atac_cutsite_residual_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 36 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_split_fold3_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_split_fold3_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_split_fold3_atac_cutsite_combined_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 37 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold1_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold1_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold1_atac_cutsite_residual_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 38 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold1_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold1_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold1_atac_cutsite_combined_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 39 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold1_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold1_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold1_atac_cutsite_residual_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 40 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold1_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold1_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold1_atac_cutsite_combined_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 41 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold1_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold1_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold1_atac_cutsite_residual_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 42 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold1_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold1_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold1_atac_cutsite_combined_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 43 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold2_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold2_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold2_atac_cutsite_residual_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 44 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold2_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold2_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold2_atac_cutsite_combined_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 45 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold2_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold2_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold2_atac_cutsite_residual_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 46 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold2_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold2_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold2_atac_cutsite_combined_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 47 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold2_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold2_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold2_atac_cutsite_residual_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 48 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold2_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold2_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold2_atac_cutsite_combined_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 49 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold3_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold3_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold3_atac_cutsite_residual_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 50 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold3_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold3_test.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold3_atac_cutsite_combined_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 51 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold3_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold3_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold3_atac_cutsite_residual_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 52 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold3_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold3_val.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold3_atac_cutsite_combined_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 53 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold3_atac_cutsite_residual.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold3_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold3_atac_cutsite_residual_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 54 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                             --reference bw/mesc_native_atac_cutsites_combined.bw                             --pred preds/atac_wt_across_acc_fold3_atac_cutsite_combined.bw                             --regions bed/bpreveal/atac_wt_across_acc_fold3_train.bed                             --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_atac_wt_across_acc_fold3_atac_cutsite_combined_train.json

fi


