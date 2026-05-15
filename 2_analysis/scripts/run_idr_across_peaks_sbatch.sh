#!/bin/bash
#SBATCH --job-name=mw2098-2026-02-02
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=200GB
#SBATCH --output=mw2098-2026-02-02.out

#!/bin/bash
source ~/.bashrc; conda deactivate; conda activate idr
~/anaconda3/envs/idr/bin/idr --peak-list macs2/Bunina_mesc_atac_day0_peaks.narrowPeak --samples ../1_processing/peaks/Bunina_mesc_atac_day0_1_peaks.narrowPeak ../1_processing/peaks/Bunina_mesc_atac_day0_2_peaks.narrowPeak --idr-threshold 0.05 --input-file-type narrowPeak --output-file idr/Bunina_mesc_atac_day0_1_vs_2_idr.txt
~/anaconda3/envs/idr/bin/idr --peak-list macs2/Bunina_mesc_atac_day4_peaks.narrowPeak --samples ../1_processing/peaks/Bunina_mesc_atac_day4_1_peaks.narrowPeak ../1_processing/peaks/Bunina_mesc_atac_day4_2_peaks.narrowPeak --idr-threshold 0.05 --input-file-type narrowPeak --output-file idr/Bunina_mesc_atac_day4_1_vs_2_idr.txt
~/anaconda3/envs/idr/bin/idr --peak-list macs2/Bunina_mesc_atac_day8_peaks.narrowPeak --samples ../1_processing/peaks/Bunina_mesc_atac_day8_1_peaks.narrowPeak ../1_processing/peaks/Bunina_mesc_atac_day8_2_peaks.narrowPeak --idr-threshold 0.05 --input-file-type narrowPeak --output-file idr/Bunina_mesc_atac_day8_1_vs_2_idr.txt
~/anaconda3/envs/idr/bin/idr --peak-list macs2/Bunina_mesc_atac_day12_peaks.narrowPeak --samples ../1_processing/peaks/Bunina_mesc_atac_day12_1_peaks.narrowPeak ../1_processing/peaks/Bunina_mesc_atac_day12_2_peaks.narrowPeak --idr-threshold 0.05 --input-file-type narrowPeak --output-file idr/Bunina_mesc_atac_day12_1_vs_2_idr.txt
