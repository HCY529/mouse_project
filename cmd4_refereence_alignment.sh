#参考基因组
bismark_genome_preparation /home/huchenyuan/mouse_project/data/reference


#比对bismark
bismark -1 ../trim_galore_data/SRR833528_1_val_1.fq -2 ../trim_galore_data/SRR833528_2_val_2.fq  /home/huchenyuan/mouse_project/data/reference/
