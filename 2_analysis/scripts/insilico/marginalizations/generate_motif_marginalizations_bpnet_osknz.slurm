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
--model_dir models/bpnet_osknz_fold1.model \
--null_sequence_np_filepath npz/random_seqs_seed_2356_input_2032_trials_256_group_mid_gc_0.45_to_0.525_cpg_0.25_to_0.5_array.npz \
--tasks_separated_by_commas oct4,sox2,klf4,nanog,zic3 \
--output_tsv_path tsv/insilico/marginalizations/motifs/marginalizations_bpnet_osknz.tsv.gz \
--seqA_txt_filepath tsv/insilico/marginalizations/motifs/all_sequences.tsv.gz \
--motifA_name all \
--motifA_output_position 500 \
--seqB_txt_filepath tsv/mapped_motifs/filler.tsv.gz \
--motifB_name empty \
--motifB_output_position 530 \
--input_length 2032 \
--output_length 1000
