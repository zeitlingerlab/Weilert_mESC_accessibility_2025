#!/bin/bash
#SBATCH --job-name marginalizations
#SBATCH --output=slurm_%j.log
#SBATCH --mem=100gb
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=10
source /home/mw2098/.bashrc
conda deactivate
conda activate bpreveal_404
cd /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis
/home/mw2098/anaconda3/envs/bpreveal_404/bin/python scripts/py/tools/generate_insilico_perturbs_across_affinities.py \
--model_dir models/atac_3h_fold1_residual.model/ \
--null_sequence_np_filepath npz/open_dinuc_shuffled_seqs_input_length_2032_seed_2356_trials_256_array.npz \
--tasks_separated_by_commas atac \
--output_tsv_path tsv/insilico/marginalizations/dinucleotides/marginalizations_atac_3h.tsv.gz \
--seqA_txt_filepath tsv/insilico/marginalizations/dinucleotides/dinucleotides.tsv.gz \
--motifA_name all \
--motifA_output_position 500 \
--seqB_txt_filepath tsv/mapped_motifs/filler.tsv.gz \
--motifB_name empty \
--motifB_output_position 530 \
--input_length 2032 \
--output_length 1000
