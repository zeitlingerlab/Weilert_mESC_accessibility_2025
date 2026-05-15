#!/usr/bin/bash
#SBATCH --job-name bpnet_training
#SBATCH --output=slurm_%j.log
#SBATCH --array=1-270%20
#SBATCH --mem=200gb
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=10
source /home/mw2098/.bashrc
conda deactivate
conda activate bpreveal_510
ml bedtools
cd /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/
if [[ ${SLURM_ARRAY_TASK_ID} == 1 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold1_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_oct4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 2 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold1_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_oct4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 3 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold1_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_sox2_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 4 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold1_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_sox2_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 5 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold1_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_nanog_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 6 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold1_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_nanog_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 7 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold1_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_klf4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 8 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold1_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_klf4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 9 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold1_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_zic3_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 10 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold1_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_zic3_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 11 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold1_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_oct4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 12 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold1_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_oct4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 13 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold1_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_sox2_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 14 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold1_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_sox2_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 15 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold1_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_nanog_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 16 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold1_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_nanog_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 17 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold1_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_klf4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 18 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold1_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_klf4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 19 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold1_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_zic3_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 20 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold1_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_zic3_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 21 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold1_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_oct4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 22 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold1_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_oct4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 23 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold1_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_sox2_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 24 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold1_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_sox2_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 25 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold1_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_nanog_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 26 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold1_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_nanog_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 27 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold1_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_klf4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 28 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold1_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_klf4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 29 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold1_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_zic3_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 30 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold1_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold1_zic3_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 31 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold2_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_oct4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 32 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold2_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_oct4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 33 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold2_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_sox2_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 34 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold2_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_sox2_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 35 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold2_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_nanog_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 36 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold2_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_nanog_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 37 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold2_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_klf4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 38 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold2_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_klf4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 39 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold2_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_zic3_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 40 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold2_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_zic3_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 41 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold2_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_oct4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 42 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold2_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_oct4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 43 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold2_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_sox2_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 44 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold2_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_sox2_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 45 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold2_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_nanog_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 46 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold2_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_nanog_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 47 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold2_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_klf4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 48 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold2_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_klf4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 49 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold2_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_zic3_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 50 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold2_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_zic3_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 51 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold2_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_oct4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 52 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold2_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_oct4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 53 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold2_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_sox2_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 54 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold2_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_sox2_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 55 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold2_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_nanog_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 56 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold2_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_nanog_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 57 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold2_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_klf4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 58 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold2_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_klf4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 59 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold2_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_zic3_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 60 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold2_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold2_zic3_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 61 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold3_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_oct4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 62 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold3_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_oct4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 63 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold3_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_sox2_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 64 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold3_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_sox2_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 65 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold3_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_nanog_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 66 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold3_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_nanog_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 67 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold3_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_klf4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 68 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold3_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_klf4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 69 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold3_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_zic3_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 70 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold3_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_zic3_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 71 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold3_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_oct4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 72 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold3_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_oct4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 73 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold3_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_sox2_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 74 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold3_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_sox2_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 75 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold3_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_nanog_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 76 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold3_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_nanog_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 77 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold3_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_klf4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 78 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold3_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_klf4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 79 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold3_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_zic3_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 80 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold3_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_zic3_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 81 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold3_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_oct4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 82 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold3_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_oct4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 83 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold3_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_sox2_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 84 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold3_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_sox2_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 85 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold3_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_nanog_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 86 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold3_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_nanog_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 87 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold3_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_klf4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 88 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold3_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_klf4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 89 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_genome_fold3_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_zic3_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 90 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_genome_fold3_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_genome_fold3_zic3_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 91 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold1_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_oct4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 92 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold1_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_oct4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 93 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold1_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_sox2_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 94 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold1_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_sox2_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 95 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold1_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_nanog_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 96 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold1_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_nanog_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 97 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold1_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_klf4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 98 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold1_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_klf4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 99 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold1_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_zic3_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 100 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold1_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_zic3_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 101 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold1_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_oct4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 102 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold1_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_oct4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 103 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold1_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_sox2_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 104 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold1_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_sox2_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 105 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold1_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_nanog_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 106 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold1_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_nanog_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 107 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold1_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_klf4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 108 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold1_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_klf4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 109 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold1_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_zic3_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 110 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold1_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_zic3_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 111 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold1_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_oct4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 112 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold1_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_oct4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 113 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold1_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_sox2_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 114 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold1_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_sox2_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 115 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold1_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_nanog_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 116 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold1_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_nanog_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 117 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold1_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_klf4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 118 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold1_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_klf4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 119 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold1_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_zic3_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 120 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold1_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold1_zic3_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 121 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold2_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_oct4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 122 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold2_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_oct4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 123 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold2_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_sox2_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 124 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold2_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_sox2_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 125 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold2_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_nanog_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 126 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold2_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_nanog_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 127 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold2_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_klf4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 128 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold2_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_klf4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 129 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold2_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_zic3_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 130 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold2_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_zic3_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 131 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold2_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_oct4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 132 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold2_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_oct4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 133 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold2_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_sox2_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 134 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold2_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_sox2_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 135 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold2_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_nanog_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 136 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold2_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_nanog_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 137 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold2_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_klf4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 138 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold2_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_klf4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 139 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold2_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_zic3_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 140 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold2_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_zic3_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 141 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold2_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_oct4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 142 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold2_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_oct4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 143 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold2_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_sox2_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 144 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold2_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_sox2_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 145 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold2_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_nanog_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 146 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold2_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_nanog_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 147 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold2_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_klf4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 148 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold2_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_klf4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 149 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold2_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_zic3_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 150 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold2_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold2_zic3_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 151 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold3_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_oct4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 152 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold3_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_oct4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 153 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold3_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_sox2_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 154 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold3_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_sox2_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 155 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold3_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_nanog_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 156 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold3_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_nanog_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 157 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold3_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_klf4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 158 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold3_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_klf4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 159 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold3_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_zic3_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 160 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold3_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_zic3_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 161 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold3_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_oct4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 162 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold3_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_oct4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 163 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold3_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_sox2_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 164 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold3_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_sox2_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 165 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold3_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_nanog_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 166 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold3_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_nanog_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 167 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold3_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_klf4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 168 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold3_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_klf4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 169 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold3_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_zic3_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 170 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold3_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_zic3_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 171 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold3_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_oct4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 172 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold3_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_oct4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 173 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold3_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_sox2_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 174 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold3_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_sox2_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 175 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold3_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_nanog_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 176 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold3_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_nanog_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 177 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold3_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_klf4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 178 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold3_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_klf4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 179 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_split_fold3_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_zic3_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 180 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_split_fold3_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_split_fold3_zic3_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 181 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold1_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_oct4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 182 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold1_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_oct4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 183 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold1_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_sox2_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 184 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold1_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_sox2_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 185 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold1_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_nanog_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 186 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold1_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_nanog_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 187 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold1_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_klf4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 188 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold1_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_klf4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 189 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold1_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_zic3_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 190 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold1_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_zic3_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 191 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold1_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_oct4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 192 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold1_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_oct4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 193 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold1_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_sox2_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 194 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold1_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_sox2_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 195 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold1_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_nanog_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 196 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold1_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_nanog_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 197 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold1_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_klf4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 198 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold1_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_klf4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 199 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold1_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_zic3_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 200 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold1_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_zic3_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 201 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold1_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_oct4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 202 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold1_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_oct4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 203 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold1_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_sox2_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 204 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold1_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_sox2_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 205 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold1_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_nanog_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 206 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold1_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_nanog_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 207 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold1_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_klf4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 208 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold1_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_klf4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 209 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold1_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_zic3_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 210 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold1_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold1_zic3_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 211 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold2_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_oct4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 212 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold2_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_oct4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 213 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold2_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_sox2_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 214 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold2_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_sox2_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 215 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold2_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_nanog_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 216 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold2_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_nanog_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 217 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold2_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_klf4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 218 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold2_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_klf4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 219 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold2_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_zic3_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 220 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold2_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_zic3_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 221 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold2_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_oct4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 222 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold2_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_oct4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 223 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold2_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_sox2_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 224 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold2_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_sox2_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 225 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold2_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_nanog_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 226 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold2_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_nanog_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 227 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold2_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_klf4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 228 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold2_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_klf4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 229 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold2_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_zic3_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 230 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold2_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_zic3_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 231 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold2_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_oct4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 232 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold2_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_oct4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 233 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold2_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_sox2_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 234 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold2_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_sox2_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 235 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold2_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_nanog_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 236 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold2_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_nanog_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 237 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold2_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_klf4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 238 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold2_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_klf4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 239 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold2_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_zic3_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 240 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold2_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold2_zic3_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 241 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold3_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_oct4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 242 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold3_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_oct4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 243 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold3_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_sox2_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 244 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold3_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_sox2_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 245 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold3_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_nanog_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 246 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold3_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_nanog_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 247 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold3_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_klf4_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 248 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold3_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_klf4_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 249 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold3_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_zic3_pos_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 250 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold3_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_zic3_neg_peaks_test.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 251 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold3_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_oct4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 252 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold3_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_oct4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 253 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold3_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_sox2_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 254 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold3_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_sox2_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 255 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold3_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_nanog_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 256 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold3_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_nanog_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 257 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold3_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_klf4_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 258 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold3_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_klf4_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 259 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold3_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_zic3_pos_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 260 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold3_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_zic3_neg_peaks_val.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 261 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold3_oct4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_oct4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 262 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_oct4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold3_oct4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_oct4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 263 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold3_sox2_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_sox2_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 264 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_sox2_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold3_sox2_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_sox2_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 265 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold3_nanog_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_nanog_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 266 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_nanog_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold3_nanog_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_nanog_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 267 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold3_klf4_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_klf4_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 268 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_klf4_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold3_klf4_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_klf4_neg_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 269 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_positive.bw                         --pred preds/bpnet_osknz_across_acc_fold3_zic3_pos.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_zic3_pos_peaks_train.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 270 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/metrics.py                         --reference bw/mesc_zic3_nexus_combined_negative.bw                         --pred preds/bpnet_osknz_across_acc_fold3_zic3_neg.bw                         --regions bed/bpreveal/atac_wt_fold1_all.bed                         --threads 70 --json-output --apply-abs --skip-zeroes | tee json/metrics_bpnet_osknz_across_acc_fold3_zic3_neg_peaks_train.json

fi


