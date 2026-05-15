#!/usr/bin/bash
#SBATCH --array=1-8%15
#SBATCH --job-name=modiscolite
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%j.log
source /home/mw2098/.bashrc
conda deactivate
ml meme
conda activate bpreveal_404
cd /n/projects/mw2098/publications/2024_weilert_acc/code/2_analysis
if [[ ${SLURM_ARRAY_TASK_ID} == 1 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/motifScan.py json/prepareMotifSeqletCutoffs_bunina_atac_day0_across_combined_atac_counts.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 2 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/motifScan.py json/prepareMotifSeqletCutoffs_bunina_atac_day0_across_combined_atac_profile.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 3 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/motifScan.py json/prepareMotifSeqletCutoffs_bunina_atac_day4_across_combined_atac_counts.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 4 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/motifScan.py json/prepareMotifSeqletCutoffs_bunina_atac_day4_across_combined_atac_profile.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 5 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/motifScan.py json/prepareMotifSeqletCutoffs_bunina_atac_day8_across_combined_atac_counts.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 6 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/motifScan.py json/prepareMotifSeqletCutoffs_bunina_atac_day8_across_combined_atac_profile.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 7 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/motifScan.py json/prepareMotifSeqletCutoffs_bunina_atac_day12_across_combined_atac_counts.json

fi


if [[ ${SLURM_ARRAY_TASK_ID} == 8 ]] ; then

    /home/mw2098/anaconda3/envs/bpreveal_404/bin/python /home/mw2098/bin/bpreveal_404//src/motifScan.py json/prepareMotifSeqletCutoffs_bunina_atac_day12_across_combined_atac_profile.json

fi


