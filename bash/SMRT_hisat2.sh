
############################## Crtl_0h_R2 #######################################################################
hisat2 -t -p 15 --qc-filter -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/Crtl_0h_R2_S7_R1_001.fastq -2 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/Crtl_0h_R2_S7_R2_001.fastq | samtools view -bS - >hisat2_out/Crtl_0h_R2.bam
samtools sort -@ 10 hisat2_out/Crtl_0h_R2.bam -o hisat2_out/Crtl_0h_R2.srt.bam
samtools index hisat2_out/Crtl_0h_R2.srt.bam

############################## Crtl_0h_R3 #######################################################################
hisat2 -t -p 15 --qc-filter -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/Crtl_0h_R3_S1_R1_001.fastq -2 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/Crtl_0h_R3_S1_R2_001.fastq | samtools view -bS - >hisat2_out/Crtl_0h_R3.bam
samtools sort -@ 10 hisat2_out/Crtl_0h_R3.bam -o hisat2_out/Crtl_0h_R3.srt.bam
samtools index hisat2_out/Crtl_0h_R3.srt.bam

############################## Ctrl_CpG_2h_R2 #######################################################################
hisat2 -t -p 15 --qc-filter -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/Ctrl_CpG_2h_R2_S9_R1_001.fastq -2 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/Ctrl_CpG_2h_R2_S9_R2_001.fastq | samtools view -bS - >hisat2_out/Ctrl_CpG_2h_R2.bam
samtools sort -@ 10 hisat2_out/Ctrl_CpG_2h_R2.bam -o hisat2_out/Ctrl_CpG_2h_R2.srt.bam
samtools index hisat2_out/Ctrl_CpG_2h_R2.srt.bam

############################## Ctrl_CpG_2h_R3 #######################################################################
hisat2 -t -p 15 --qc-filter -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/Ctrl_CpG_2h_R3_S3_R1_001.fastq -2 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/Ctrl_CpG_2h_R3_S3_R2_001.fastq | samtools view -bS - >hisat2_out/Ctrl_CpG_2h_R3.bam
samtools sort -@ 10 hisat2_out/Ctrl_CpG_2h_R3.bam -o hisat2_out/Ctrl_CpG_2h_R3.srt.bam
samtools index hisat2_out/Ctrl_CpG_2h_R3.srt.bam

############################## Ctrl_CpG_6h_R2 #######################################################################
hisat2 -t -p 15 --qc-filter -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/Ctrl_CpG_6h_R2_S11_R1_001.fastq -2 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/Ctrl_CpG_6h_R2_S11_R2_001.fastq | samtools view -bS - >hisat2_out/Ctrl_CpG_6h_R2.bam
samtools sort -@ 10 hisat2_out/Ctrl_CpG_6h_R2.bam -o hisat2_out/Ctrl_CpG_6h_R2.srt.bam
samtools index hisat2_out/Ctrl_CpG_6h_R2.srt.bam

############################## Ctrl_CpG_6h_R3 #######################################################################
hisat2 -t -p 15 --qc-filter -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/Ctrl_CpG_6h_R3_S5_R1_001.fastq -2 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/Ctrl_CpG_6h_R3_S5_R2_001.fastq | samtools view -bS - >hisat2_out/Ctrl_CpG_6h_R3.bam
samtools sort -@ 10 hisat2_out/Ctrl_CpG_6h_R3.bam -o hisat2_out/Ctrl_CpG_6h_R3.srt.bam
samtools index hisat2_out/Ctrl_CpG_6h_R3.srt.bam

############################## SMRT_0h_R2 #######################################################################
hisat2 -t -p 15 --qc-filter -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/SMRT_0h_R2_S8_R1_001.fastq -2 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/SMRT_0h_R2_S8_R2_001.fastq | samtools view -bS - >hisat2_out/SMRT_0h_R2.bam
samtools sort -@ 10 hisat2_out/SMRT_0h_R2.bam -o hisat2_out/SMRT_0h_R2.srt.bam
samtools index hisat2_out/SMRT_0h_R2.srt.bam

############################## SMRT_0h_R3 #######################################################################
hisat2 -t -p 15 --qc-filter -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/SMRT_0h_R3_S2_R1_001.fastq -2 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/SMRT_0h_R3_S2_R2_001.fastq | samtools view -bS - >hisat2_out/SMRT_0h_R3.bam
samtools sort -@ 10 hisat2_out/SMRT_0h_R3.bam -o hisat2_out/SMRT_0h_R3.srt.bam
samtools index hisat2_out/SMRT_0h_R3.srt.bam

############################## SMRT_CpG_2h_R2 #######################################################################
hisat2 -t -p 15 --qc-filter -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/SMRT_CpG_2h_R2_S10_R1_001.fastq -2 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/SMRT_CpG_2h_R2_S10_R2_001.fastq | samtools view -bS - >hisat2_out/SMRT_CpG_2h_R2.bam
samtools sort -@ 10 hisat2_out/SMRT_CpG_2h_R2.bam -o hisat2_out/SMRT_CpG_2h_R2.srt.bam
samtools index hisat2_out/SMRT_CpG_2h_R2.srt.bam

############################## SMRT_CpG_2h_R3 #######################################################################
hisat2 -t -p 15 --qc-filter -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/SMRT_CpG_2h_R3_S4_R1_001.fastq -2 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/SMRT_CpG_2h_R3_S4_R2_001.fastq | samtools view -bS - >hisat2_out/SMRT_CpG_2h_R3.bam
samtools sort -@ 10 hisat2_out/SMRT_CpG_2h_R3.bam -o hisat2_out/SMRT_CpG_2h_R3.srt.bam
samtools index hisat2_out/SMRT_CpG_2h_R3.srt.bam

############################## SMRT_CpG_6h_R2 #######################################################################
hisat2 -t -p 15 --qc-filter -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/SMRT_CpG_6h_R2_S12_R1_001.fastq -2 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/SMRT_CpG_6h_R2_S12_R2_001.fastq | samtools view -bS - >hisat2_out/SMRT_CpG_6h_R2.bam
samtools sort -@ 10 hisat2_out/SMRT_CpG_6h_R2.bam -o hisat2_out/SMRT_CpG_6h_R2.srt.bam
samtools index hisat2_out/SMRT_CpG_6h_R2.srt.bam

############################## SMRT_CpG_6h_R3 #######################################################################
hisat2 -t -p 15 --qc-filter -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/SMRT_CpG_6h_R3_S6_R1_001.fastq -2 /media/Hard_disk1/RAW_data/NCOR1/ILS_seq_data/190202_NB551648_0003_AHT3H2BGX9/fastq_files/SMRT_CpG_6h_R3_S6_R2_001.fastq | samtools view -bS - >hisat2_out/SMRT_CpG_6h_R3.bam
samtools sort -@ 10 hisat2_out/SMRT_CpG_6h_R3.bam -o hisat2_out/SMRT_CpG_6h_R3.srt.bam
samtools index hisat2_out/SMRT_CpG_6h_R3.srt.bam
