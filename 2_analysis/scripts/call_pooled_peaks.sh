#!/bin/bash
ml zeitlinger
macs2 callpeak -f BAMPE -t bam/GSE174774_mesc_atac_0h_combined.bam --outdir macs2 -n GSE174774_mesc_atac_0h
macs2 callpeak -f BAMPE -t bam/GSE174774_mesc_atac_12h_combined.bam --outdir macs2 -n GSE174774_mesc_atac_12h
macs2 callpeak -f BAMPE -t bam/GSE174774_mesc_atac_15h_combined.bam --outdir macs2 -n GSE174774_mesc_atac_15h
macs2 callpeak -f BAMPE -t bam/GSE174774_mesc_atac_3h_combined.bam --outdir macs2 -n GSE174774_mesc_atac_3h
macs2 callpeak -f BAMPE -t bam/GSE174774_mesc_atac_6h_combined.bam --outdir macs2 -n GSE174774_mesc_atac_6h
macs2 callpeak -f BAMPE -t bam/GSE174774_mesc_atac_9h_combined.bam --outdir macs2 -n GSE174774_mesc_atac_9h
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/GSE174774_mesc_ttseq_0h_combined.bam --outdir macs2 -n GSE174774_mesc_ttseq_0h
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/GSE174774_mesc_ttseq_12h_combined.bam --outdir macs2 -n GSE174774_mesc_ttseq_12h
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/GSE174774_mesc_ttseq_15h_combined.bam --outdir macs2 -n GSE174774_mesc_ttseq_15h
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/GSE174774_mesc_ttseq_3h_combined.bam --outdir macs2 -n GSE174774_mesc_ttseq_3h
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/GSE174774_mesc_ttseq_6h_combined.bam --outdir macs2 -n GSE174774_mesc_ttseq_6h
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/GSE174774_mesc_ttseq_9h_combined.bam --outdir macs2 -n GSE174774_mesc_ttseq_9h
macs2 callpeak -f BAMPE -t bam/mesc_Akr1cl_scenario_context_mut_A_atac_combined.bam --outdir macs2 -n mesc_Akr1cl_scenario_context_mut_A_atac
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/mesc_Akr1cl_scenario_context_mut_A_klf4_nexus_combined.bam --outdir macs2 -n mesc_Akr1cl_scenario_context_mut_A_klf4_nexus
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/mesc_Akr1cl_scenario_context_mut_A_sox2_nexus_combined.bam --outdir macs2 -n mesc_Akr1cl_scenario_context_mut_A_sox2_nexus
macs2 callpeak -f BAMPE -t bam/mesc_Akr1cl_scenario_context_mut_B_atac_combined.bam --outdir macs2 -n mesc_Akr1cl_scenario_context_mut_B_atac
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/mesc_Akr1cl_scenario_context_mut_B_klf4_nexus_combined.bam --outdir macs2 -n mesc_Akr1cl_scenario_context_mut_B_klf4_nexus
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/mesc_Akr1cl_scenario_context_mut_B_sox2_nexus_combined.bam --outdir macs2 -n mesc_Akr1cl_scenario_context_mut_B_sox2_nexus
macs2 callpeak -f BAMPE -t bam/mesc_Btbd11_scenario_WT_coop_WT_atac_combined.bam --outdir macs2 -n mesc_Btbd11_scenario_WT_coop_WT_atac
macs2 callpeak -f BAMPE -t bam/mesc_Btbd11_scenario_WT_coop_mut_A_atac_combined.bam --outdir macs2 -n mesc_Btbd11_scenario_WT_coop_mut_A_atac
macs2 callpeak -f BAMPE -t bam/mesc_Btbd11_scenario_WT_coop_mut_B_atac_combined.bam --outdir macs2 -n mesc_Btbd11_scenario_WT_coop_mut_B_atac
macs2 callpeak -f BAMPE -t bam/mesc_Btbd11_scenario_WT_coop_mut_null_atac_combined.bam --outdir macs2 -n mesc_Btbd11_scenario_WT_coop_mut_null_atac
macs2 callpeak -f BAMPE -t bam/mesc_Btbd11_scenario_enh_coop_WT_atac_combined.bam --outdir macs2 -n mesc_Btbd11_scenario_enh_coop_WT_atac
macs2 callpeak -f BAMPE -t bam/mesc_Btbd11_scenario_enh_coop_mut_A_atac_combined.bam --outdir macs2 -n mesc_Btbd11_scenario_enh_coop_mut_A_atac
macs2 callpeak -f BAMPE -t bam/mesc_Btbd11_scenario_enh_coop_mut_B_atac_combined.bam --outdir macs2 -n mesc_Btbd11_scenario_enh_coop_mut_B_atac
macs2 callpeak -f BAMPE -t bam/mesc_Btbd11_scenario_enh_coop_mut_null_atac_combined.bam --outdir macs2 -n mesc_Btbd11_scenario_enh_coop_mut_null_atac
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/mesc_R1_klf4_nexus_combined.bam --outdir macs2 -n mesc_R1_klf4_nexus
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/mesc_R1_sox2_nexus_combined.bam --outdir macs2 -n mesc_R1_sox2_nexus
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/mesc_ctcf_nexus_combined.bam --outdir macs2 -n mesc_ctcf_nexus
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/mesc_esrrb_nexus_combined.bam --outdir macs2 -n mesc_esrrb_nexus
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/mesc_klf4_nexus_combined.bam --outdir macs2 -n mesc_klf4_nexus
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/mesc_nanog_nexus_combined.bam --outdir macs2 -n mesc_nanog_nexus
macs2 callpeak -f BAMPE -t bam/mesc_native_atac_combined.bam --outdir macs2 -n mesc_native_atac
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/mesc_oct4_nexus_combined.bam --outdir macs2 -n mesc_oct4_nexus
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/mesc_sox2_nexus_combined.bam --outdir macs2 -n mesc_sox2_nexus
macs2 callpeak -f BAM --keep-dup all --nomodel --shift -75 --extsize 150 -t bam/mesc_zic3_nexus_combined.bam --outdir macs2 -n mesc_zic3_nexus
