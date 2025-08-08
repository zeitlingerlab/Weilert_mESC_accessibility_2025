#!/bin/bash
#SBATCH --job-name genomic_perturbations
#SBATCH --output=slurm_%j.log
#SBATCH --mem=500gb
#SBATCH --time=72:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=30
source /home/mw2098/.bashrc
conda deactivate
conda activate bpreveal_404
cd /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis
/home/mw2098/anaconda3/envs/bpreveal_404/bin/python scripts/py/tools/generate_genomic_perturbs.py \
--motifs_df_path tsv/genomic/motifs/instances_formatted_for_genomic_perturbations_0based.tsv.gz \
--motif_groups_df_path tsv/genomic/motifs/regions_formatted_for_genomic_perturbations_0based.tsv.gz \
--model_dir models/atac_0h_fold1_residual.model/ \
--output_prefix tsv/genomic/motifs/genomic_pertubs_all_motifs_atac_0h \
--tasks_in_model_order atac \
--region_index_column region_id \
--motif_window_start_column motif_window_start \
--motif_window_end_column motif_window_end \
--motif_name_column motif \
--motif_chrom_column chrom \
--fasta_file ../0_setup/fa/mm10.fa \
--input_seqlen 2032 --output_seqlen 1000 \
--comb_max island_count --comb_min island_count \
--nodes 24 --trials 16 \
--use_whole_window 
