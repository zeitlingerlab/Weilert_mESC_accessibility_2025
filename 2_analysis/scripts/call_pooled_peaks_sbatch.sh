#!/bin/bash
#SBATCH --job-name=mw2098-2026-02-02
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=200GB
#SBATCH --output=mw2098-2026-02-02.out

#!/bin/bash
ml zeitlinger
macs2 callpeak -f BAMPE -t bam/Bunina_mesc_atac_day0_combined.bam --outdir macs2 -n Bunina_mesc_atac_day0
macs2 callpeak -f BAMPE -t bam/Bunina_mesc_atac_day12_combined.bam --outdir macs2 -n Bunina_mesc_atac_day12
macs2 callpeak -f BAMPE -t bam/Bunina_mesc_atac_day4_combined.bam --outdir macs2 -n Bunina_mesc_atac_day4
macs2 callpeak -f BAMPE -t bam/Bunina_mesc_atac_day8_combined.bam --outdir macs2 -n Bunina_mesc_atac_day8
