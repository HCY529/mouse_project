#排序
samtools sort -n ../SRR833527_1_val_1_bismark_bt2_pe.deduplicated.bam



#甲基化信息提取
bismark_methylation_extractor --comprehensive  SRR833527_1_val_1_bismark_bt2_pe.deduplicated.bam
