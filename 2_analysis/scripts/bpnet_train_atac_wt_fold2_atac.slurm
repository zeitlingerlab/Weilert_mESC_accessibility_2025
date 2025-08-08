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
conda activate /home/mw2098/anaconda3/envs/acc_vs_bind
cd /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis
/home/mw2098/anaconda3/envs/acc_vs_bind/bin/python /n/projects/mw2098/publications/2024_weilert_acc/public/software/bpreveal//src/interpretFlat.py json/shapFlat_atac_wt_fold2_atac.json
