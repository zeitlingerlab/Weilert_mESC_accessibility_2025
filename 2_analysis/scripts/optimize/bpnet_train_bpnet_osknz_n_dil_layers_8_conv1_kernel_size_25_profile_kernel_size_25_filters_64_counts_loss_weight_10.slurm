#!/usr/bin/bash
#SBATCH --job-name bpnet_training
#SBATCH --output=slurm_%j.log
#SBATCH --mem=200gb
#SBATCH --time=72:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=10
source /home/mw2098/.bashrc
conda deactivate
conda activate bpreveal_404
cd /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis
/home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/trainSoloModel.py json/optimize/trainSoloModel_bpnet_osknz_n_dil_layers_8_conv1_kernel_size_25_profile_kernel_size_25_filters_64_counts_loss_weight_10.json
/home/mw2098/anaconda3/envs/bpreveal_404/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal_404//src/makePredictionsBed.py json/optimize/predictPeaks_bpnet_osknz_n_dil_layers_8_conv1_kernel_size_25_profile_kernel_size_25_filters_64_counts_loss_weight_10.json
