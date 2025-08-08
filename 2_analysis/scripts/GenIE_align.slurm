#!/bin/bash

#SBATCH --job-name GenIE_align
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=GenIE_align_log%j.log
#SBATCH --cpus-per-task 21
#SBATCH --array=1-32%50

source /home/hj2486/.bashrc
conda deactivate
ml zeitlinger

cd /n/projects/hj2486/analysis/mw2098/GenIEATAC/
if [[ ${SLURM_ARRAY_TASK_ID} == 1 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_AB/mm10_mesc_r1_Btbd11_WT_coop_AB -1 3_data/fa/demultiplexed/GenIEATAC_wt_SOI_forw_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_wt_SOI_forw_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_SOI_forw.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_SOI_forw.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_SOI_forw.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 2 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_A/mm10_mesc_r1_Btbd11_WT_coop_A -1 3_data/fa/demultiplexed/GenIEATAC_A_SOI_forw_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_A_SOI_forw_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_SOI_forw.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_SOI_forw.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_SOI_forw.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 3 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_B/mm10_mesc_r1_Btbd11_WT_coop_B -1 3_data/fa/demultiplexed/GenIEATAC_B_SOI_forw_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_B_SOI_forw_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_SOI_forw.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_SOI_forw.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_SOI_forw.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 4 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_null/mm10_mesc_r1_Btbd11_WT_coop_null -1 3_data/fa/demultiplexed/GenIEATAC_null_SOI_forw_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_null_SOI_forw_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_SOI_forw.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_SOI_forw.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_SOI_forw.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 5 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_AB/mm10_mesc_r1_Btbd11_WT_coop_AB -1 3_data/fa/demultiplexed/GenIEATAC_wt_SOI_revs_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_wt_SOI_revs_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_SOI_revs.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_SOI_revs.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_SOI_revs.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 6 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_A/mm10_mesc_r1_Btbd11_WT_coop_A -1 3_data/fa/demultiplexed/GenIEATAC_A_SOI_revs_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_A_SOI_revs_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_SOI_revs.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_SOI_revs.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_SOI_revs.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 7 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_B/mm10_mesc_r1_Btbd11_WT_coop_B -1 3_data/fa/demultiplexed/GenIEATAC_B_SOI_revs_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_B_SOI_revs_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_SOI_revs.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_SOI_revs.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_SOI_revs.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 8 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_null/mm10_mesc_r1_Btbd11_WT_coop_null -1 3_data/fa/demultiplexed/GenIEATAC_null_SOI_revs_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_null_SOI_revs_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_SOI_revs.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_SOI_revs.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_SOI_revs.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 9 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_AB/mm10_mesc_r1_Btbd11_WT_coop_AB -1 3_data/fa/demultiplexed/GenIEATAC_wt_Adj_forw_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_wt_Adj_forw_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_Adj_forw.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_Adj_forw.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_Adj_forw.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 10 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_A/mm10_mesc_r1_Btbd11_WT_coop_A -1 3_data/fa/demultiplexed/GenIEATAC_A_Adj_forw_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_A_Adj_forw_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_Adj_forw.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_Adj_forw.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_Adj_forw.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 11 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_B/mm10_mesc_r1_Btbd11_WT_coop_B -1 3_data/fa/demultiplexed/GenIEATAC_B_Adj_forw_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_B_Adj_forw_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_Adj_forw.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_Adj_forw.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_Adj_forw.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 12 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_null/mm10_mesc_r1_Btbd11_WT_coop_null -1 3_data/fa/demultiplexed/GenIEATAC_null_Adj_forw_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_null_Adj_forw_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_Adj_forw.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_Adj_forw.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_Adj_forw.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 13 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_AB/mm10_mesc_r1_Btbd11_WT_coop_AB -1 3_data/fa/demultiplexed/GenIEATAC_wt_Adj_revs_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_wt_Adj_revs_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_Adj_revs.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_Adj_revs.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_Adj_revs.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 14 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_A/mm10_mesc_r1_Btbd11_WT_coop_A -1 3_data/fa/demultiplexed/GenIEATAC_A_Adj_revs_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_A_Adj_revs_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_Adj_revs.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_Adj_revs.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_Adj_revs.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 15 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_B/mm10_mesc_r1_Btbd11_WT_coop_B -1 3_data/fa/demultiplexed/GenIEATAC_B_Adj_revs_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_B_Adj_revs_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_Adj_revs.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_Adj_revs.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_Adj_revs.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 16 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_null/mm10_mesc_r1_Btbd11_WT_coop_null -1 3_data/fa/demultiplexed/GenIEATAC_null_Adj_revs_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_null_Adj_revs_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_Adj_revs.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_Adj_revs.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_Adj_revs.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 17 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_AB/mm10_mesc_r1_Btbd11_WT_coop_AB -1 3_data/fa/demultiplexed/GenIEATAC_wt_null1_forw_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_wt_null1_forw_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_null1_forw.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_null1_forw.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_null1_forw.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 18 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_A/mm10_mesc_r1_Btbd11_WT_coop_A -1 3_data/fa/demultiplexed/GenIEATAC_A_null1_forw_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_A_null1_forw_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_null1_forw.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_null1_forw.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_null1_forw.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 19 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_B/mm10_mesc_r1_Btbd11_WT_coop_B -1 3_data/fa/demultiplexed/GenIEATAC_B_null1_forw_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_B_null1_forw_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_null1_forw.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_null1_forw.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_null1_forw.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 20 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_null/mm10_mesc_r1_Btbd11_WT_coop_null -1 3_data/fa/demultiplexed/GenIEATAC_null_null1_forw_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_null_null1_forw_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_null1_forw.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_null1_forw.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_null1_forw.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 21 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_AB/mm10_mesc_r1_Btbd11_WT_coop_AB -1 3_data/fa/demultiplexed/GenIEATAC_wt_null1_revs_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_wt_null1_revs_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_null1_revs.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_null1_revs.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_null1_revs.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 22 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_A/mm10_mesc_r1_Btbd11_WT_coop_A -1 3_data/fa/demultiplexed/GenIEATAC_A_null1_revs_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_A_null1_revs_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_null1_revs.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_null1_revs.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_null1_revs.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 23 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_B/mm10_mesc_r1_Btbd11_WT_coop_B -1 3_data/fa/demultiplexed/GenIEATAC_B_null1_revs_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_B_null1_revs_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_null1_revs.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_null1_revs.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_null1_revs.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 24 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_null/mm10_mesc_r1_Btbd11_WT_coop_null -1 3_data/fa/demultiplexed/GenIEATAC_null_null1_revs_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_null_null1_revs_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_null1_revs.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_null1_revs.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_null1_revs.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 25 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_AB/mm10_mesc_r1_Btbd11_WT_coop_AB -1 3_data/fa/demultiplexed/GenIEATAC_wt_null2_forw_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_wt_null2_forw_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_null2_forw.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_null2_forw.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_null2_forw.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 26 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_A/mm10_mesc_r1_Btbd11_WT_coop_A -1 3_data/fa/demultiplexed/GenIEATAC_A_null2_forw_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_A_null2_forw_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_null2_forw.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_null2_forw.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_null2_forw.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 27 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_B/mm10_mesc_r1_Btbd11_WT_coop_B -1 3_data/fa/demultiplexed/GenIEATAC_B_null2_forw_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_B_null2_forw_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_null2_forw.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_null2_forw.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_null2_forw.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 28 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_null/mm10_mesc_r1_Btbd11_WT_coop_null -1 3_data/fa/demultiplexed/GenIEATAC_null_null2_forw_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_null_null2_forw_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_null2_forw.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_null2_forw.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_null2_forw.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 29 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_AB/mm10_mesc_r1_Btbd11_WT_coop_AB -1 3_data/fa/demultiplexed/GenIEATAC_wt_null2_revs_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_wt_null2_revs_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_null2_revs.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_null2_revs.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_wt_null2_revs.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 30 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_A/mm10_mesc_r1_Btbd11_WT_coop_A -1 3_data/fa/demultiplexed/GenIEATAC_A_null2_revs_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_A_null2_revs_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_null2_revs.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_null2_revs.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_A_null2_revs.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 31 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_B/mm10_mesc_r1_Btbd11_WT_coop_B -1 3_data/fa/demultiplexed/GenIEATAC_B_null2_revs_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_B_null2_revs_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_null2_revs.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_null2_revs.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_B_null2_revs.bam -@ 20 
fi 


if [[ ${SLURM_ARRAY_TASK_ID} == 32 ]] ; then

    bowtie2 -x /l/Zeitlinger/ZeitlingerLab/pipeline/data/supplemental_data/indexes/bowtie2/mm10_mesc_r1_Btbd11_WT_coop_null/mm10_mesc_r1_Btbd11_WT_coop_null -1 3_data/fa/demultiplexed/GenIEATAC_null_null2_revs_R1_demultiplexed.fastq.gz -2 3_data/fa/demultiplexed/GenIEATAC_null_null2_revs_R2_demultiplexed.fastq.gz \
    -X 600 -p 20 --no-unal --no-mixed 2> /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_null2_revs.stats | \
    samtools view -F 4 -Sb -h | samtools sort --threads 20 -o /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_null2_revs.bam - 
    samtools index /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis/bam/GenIEATAC_null_null2_revs.bam -@ 20 
fi 

