#!/bin/bash
#SBATCH --job-name=mw2098-2024-11-05
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=200GB
#SBATCH --output=mw2098-2024-11-05.out

#!/bin/bash
ml zeitlinger
samtools merge --threads 12 bam/GSE174774_mesc_atac_0h_combined.bam ../1_processing/bam/GSE174774_mesc_atac_0h_1.bam ../1_processing/bam/GSE174774_mesc_atac_0h_2.bam
samtools merge --threads 12 bam/GSE174774_mesc_atac_12h_combined.bam ../1_processing/bam/GSE174774_mesc_atac_12h_1.bam ../1_processing/bam/GSE174774_mesc_atac_12h_2.bam
samtools merge --threads 12 bam/GSE174774_mesc_atac_15h_combined.bam ../1_processing/bam/GSE174774_mesc_atac_15h_1.bam ../1_processing/bam/GSE174774_mesc_atac_15h_2.bam
samtools merge --threads 12 bam/GSE174774_mesc_atac_3h_combined.bam ../1_processing/bam/GSE174774_mesc_atac_3h_1.bam ../1_processing/bam/GSE174774_mesc_atac_3h_2.bam
samtools merge --threads 12 bam/GSE174774_mesc_atac_6h_combined.bam ../1_processing/bam/GSE174774_mesc_atac_6h_1.bam ../1_processing/bam/GSE174774_mesc_atac_6h_2.bam
samtools merge --threads 12 bam/GSE174774_mesc_atac_9h_combined.bam ../1_processing/bam/GSE174774_mesc_atac_9h_1.bam ../1_processing/bam/GSE174774_mesc_atac_9h_2.bam
samtools merge --threads 12 bam/GSE174774_mesc_ttseq_0h_combined.bam ../1_processing/bam/GSE174774_mesc_ttseq_0h_1.bam ../1_processing/bam/GSE174774_mesc_ttseq_0h_2.bam
samtools merge --threads 12 bam/GSE174774_mesc_ttseq_12h_combined.bam ../1_processing/bam/GSE174774_mesc_ttseq_12h_1.bam ../1_processing/bam/GSE174774_mesc_ttseq_12h_2.bam
samtools merge --threads 12 bam/GSE174774_mesc_ttseq_15h_combined.bam ../1_processing/bam/GSE174774_mesc_ttseq_15h_1.bam ../1_processing/bam/GSE174774_mesc_ttseq_15h_2.bam
samtools merge --threads 12 bam/GSE174774_mesc_ttseq_3h_combined.bam ../1_processing/bam/GSE174774_mesc_ttseq_3h_1.bam ../1_processing/bam/GSE174774_mesc_ttseq_3h_2.bam
samtools merge --threads 12 bam/GSE174774_mesc_ttseq_6h_combined.bam ../1_processing/bam/GSE174774_mesc_ttseq_6h_1.bam ../1_processing/bam/GSE174774_mesc_ttseq_6h_2.bam
samtools merge --threads 12 bam/GSE174774_mesc_ttseq_9h_combined.bam ../1_processing/bam/GSE174774_mesc_ttseq_9h_1.bam ../1_processing/bam/GSE174774_mesc_ttseq_9h_2.bam
samtools merge --threads 12 bam/mesc_Akr1cl_scenario_context_mut_A_atac_combined.bam ../1_processing/bam/mesc_Akr1cl_scenario_context_mut_A_atac_1.bam ../1_processing/bam/mesc_Akr1cl_scenario_context_mut_A_atac_2.bam
samtools merge --threads 12 bam/mesc_Akr1cl_scenario_context_mut_B_atac_combined.bam ../1_processing/bam/mesc_Akr1cl_scenario_context_mut_B_atac_1.bam ../1_processing/bam/mesc_Akr1cl_scenario_context_mut_B_atac_2.bam
samtools merge --threads 12 bam/mesc_Btbd11_scenario_WT_coop_WT_atac_combined.bam ../1_processing/bam/mesc_Btbd11_scenario_WT_coop_WT_atac_1.bam ../1_processing/bam/mesc_Btbd11_scenario_WT_coop_WT_atac_2.bam
samtools merge --threads 12 bam/mesc_Btbd11_scenario_WT_coop_mut_A_atac_combined.bam ../1_processing/bam/mesc_Btbd11_scenario_WT_coop_mut_A_atac_1.bam ../1_processing/bam/mesc_Btbd11_scenario_WT_coop_mut_A_atac_2.bam
samtools merge --threads 12 bam/mesc_Btbd11_scenario_WT_coop_mut_B_atac_combined.bam ../1_processing/bam/mesc_Btbd11_scenario_WT_coop_mut_B_atac_1.bam ../1_processing/bam/mesc_Btbd11_scenario_WT_coop_mut_B_atac_2.bam
samtools merge --threads 12 bam/mesc_Btbd11_scenario_WT_coop_mut_null_atac_combined.bam ../1_processing/bam/mesc_Btbd11_scenario_WT_coop_mut_null_atac_1.bam ../1_processing/bam/mesc_Btbd11_scenario_WT_coop_mut_null_atac_2.bam
samtools merge --threads 12 bam/mesc_Btbd11_scenario_enh_coop_WT_atac_combined.bam ../1_processing/bam/mesc_Btbd11_scenario_enh_coop_WT_atac_1.bam ../1_processing/bam/mesc_Btbd11_scenario_enh_coop_WT_atac_2.bam
samtools merge --threads 12 bam/mesc_Btbd11_scenario_enh_coop_mut_A_atac_combined.bam ../1_processing/bam/mesc_Btbd11_scenario_enh_coop_mut_A_atac_1.bam ../1_processing/bam/mesc_Btbd11_scenario_enh_coop_mut_A_atac_2.bam
samtools merge --threads 12 bam/mesc_Btbd11_scenario_enh_coop_mut_B_atac_combined.bam ../1_processing/bam/mesc_Btbd11_scenario_enh_coop_mut_B_atac_1.bam ../1_processing/bam/mesc_Btbd11_scenario_enh_coop_mut_B_atac_2.bam
samtools merge --threads 12 bam/mesc_Btbd11_scenario_enh_coop_mut_null_atac_combined.bam ../1_processing/bam/mesc_Btbd11_scenario_enh_coop_mut_null_atac_1.bam ../1_processing/bam/mesc_Btbd11_scenario_enh_coop_mut_null_atac_2.bam
samtools merge --threads 12 bam/mesc_ctcf_nexus_combined.bam ../1_processing/bam/mesc_ctcf_nexus_2.bam ../1_processing/bam/mesc_ctcf_nexus_5.bam ../1_processing/bam/mesc_ctcf_nexus_6.bam
samtools merge --threads 12 bam/mesc_esrrb_nexus_combined.bam ../1_processing/bam/mesc_esrrb_nexus_1.bam ../1_processing/bam/mesc_esrrb_nexus_2.bam ../1_processing/bam/mesc_esrrb_nexus_4.bam ../1_processing/bam/mesc_esrrb_nexus_5.bam ../1_processing/bam/mesc_esrrb_nexus_6.bam ../1_processing/bam/mesc_esrrb_nexus_7.bam ../1_processing/bam/mesc_esrrb_nexus_8.bam
samtools merge --threads 12 bam/mesc_klf4_nexus_combined.bam ../1_processing/bam/mesc_klf4_nexus_1.bam ../1_processing/bam/mesc_klf4_nexus_2.bam ../1_processing/bam/mesc_klf4_nexus_3.bam ../1_processing/bam/mesc_klf4_nexus_4.bam ../1_processing/bam/mesc_klf4_nexus_5.bam
samtools merge --threads 12 bam/mesc_nanog_nexus_combined.bam ../1_processing/bam/mesc_nanog_nexus_1.bam ../1_processing/bam/mesc_nanog_nexus_2.bam ../1_processing/bam/mesc_nanog_nexus_3.bam ../1_processing/bam/mesc_nanog_nexus_4.bam ../1_processing/bam/mesc_nanog_nexus_5.bam
samtools merge --threads 12 bam/mesc_native_atac_combined.bam ../1_processing/bam/mesc_native_atac_1.bam ../1_processing/bam/mesc_native_atac_2.bam ../1_processing/bam/mesc_native_atac_3.bam
samtools merge --threads 12 bam/mesc_oct4_nexus_combined.bam ../1_processing/bam/mesc_oct4_nexus_1.bam ../1_processing/bam/mesc_oct4_nexus_2.bam ../1_processing/bam/mesc_oct4_nexus_3.bam ../1_processing/bam/mesc_oct4_nexus_4.bam ../1_processing/bam/mesc_oct4_nexus_5.bam ../1_processing/bam/mesc_oct4_nexus_6.bam
samtools merge --threads 12 bam/mesc_sox2_nexus_combined.bam ../1_processing/bam/mesc_sox2_nexus_1.bam ../1_processing/bam/mesc_sox2_nexus_2.bam ../1_processing/bam/mesc_sox2_nexus_4.bam ../1_processing/bam/mesc_sox2_nexus_6.bam
samtools merge --threads 12 bam/mesc_sox2crispr1_b11_atac_combined.bam ../1_processing/bam/mesc_sox2crispr1_b11_atac_1.bam ../1_processing/bam/mesc_sox2crispr1_b11_atac_2.bam
samtools merge --threads 12 bam/mesc_wt_r1_atac_combined.bam ../1_processing/bam/mesc_wt_r1_atac_1.bam ../1_processing/bam/mesc_wt_r1_atac_2.bam
samtools merge --threads 12 bam/mesc_zic3_nexus_combined.bam ../1_processing/bam/mesc_zic3_nexus_1.bam ../1_processing/bam/mesc_zic3_nexus_2.bam
