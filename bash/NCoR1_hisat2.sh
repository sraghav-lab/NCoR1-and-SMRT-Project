
############################## Emp_CpG_R1 #######################################################################
hisat2 -t -p 20  -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1  Emp_CpG_R1_1.fastq.gz -2  Emp_CpG_R1_2.fastq.gz | samtools view -bS - >hisat2_out/Emp_CpG_R1.bam
samtools sort -@ 10 hisat2_out/Emp_CpG_R1.bam -o hisat2_out/Emp_CpG_R1.srt.bam
samtools index hisat2_out/Emp_CpG_R1.srt.bam

############################## Emp_CpG_R2 #######################################################################
hisat2 -t -p 20  -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1  Emp_CpG_R2_1.fastq.gz -2  Emp_CpG_R2_2.fastq.gz | samtools view -bS - >hisat2_out/Emp_CpG_R2.bam
samtools sort -@ 10 hisat2_out/Emp_CpG_R2.bam -o hisat2_out/Emp_CpG_R2.srt.bam
samtools index hisat2_out/Emp_CpG_R2.srt.bam

############################## Emp_Uns_R1 #######################################################################
hisat2 -t -p 20  -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1  Emp_Uns_R1_1.fastq.gz -2  Emp_Uns_R1_2.fastq.gz | samtools view -bS - >hisat2_out/Emp_Uns_R1.bam
samtools sort -@ 10 hisat2_out/Emp_Uns_R1.bam -o hisat2_out/Emp_Uns_R1.srt.bam
samtools index hisat2_out/Emp_Uns_R1.srt.bam

############################## Emp_Uns_R2 #######################################################################
hisat2 -t -p 20  -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1  Emp_Uns_R2_1.fastq.gz -2  Emp_Uns_R2_2.fastq.gz | samtools view -bS - >hisat2_out/Emp_Uns_R2.bam
samtools sort -@ 10 hisat2_out/Emp_Uns_R2.bam -o hisat2_out/Emp_Uns_R2.srt.bam
samtools index hisat2_out/Emp_Uns_R2.srt.bam

############################## NCoR1_CpG_R1 #######################################################################
hisat2 -t -p 20  -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1  NCoR1_CpG_R1_1.fastq.gz -2  NCoR1_CpG_R1_2.fastq.gz | samtools view -bS - >hisat2_out/NCoR1_CpG_R1.bam
samtools sort -@ 10 hisat2_out/NCoR1_CpG_R1.bam -o hisat2_out/NCoR1_CpG_R1.srt.bam
samtools index hisat2_out/NCoR1_CpG_R1.srt.bam

############################## NCoR1_CpG_R2 #######################################################################
hisat2 -t -p 20  -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1  NCoR1_CpG_R2_1.fastq.gz -2  NCoR1_CpG_R2_2.fastq.gz | samtools view -bS - >hisat2_out/NCoR1_CpG_R2.bam
samtools sort -@ 10 hisat2_out/NCoR1_CpG_R2.bam -o hisat2_out/NCoR1_CpG_R2.srt.bam
samtools index hisat2_out/NCoR1_CpG_R2.srt.bam

############################## NCoR1_Uns_R1 #######################################################################
hisat2 -t -p 20  -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1  NCoR1_Uns_R1_1.fastq.gz -2  NCoR1_Uns_R1_2.fastq.gz | samtools view -bS - >hisat2_out/NCoR1_Uns_R1.bam
samtools sort -@ 10 hisat2_out/NCoR1_Uns_R1.bam -o hisat2_out/NCoR1_Uns_R1.srt.bam
samtools index hisat2_out/NCoR1_Uns_R1.srt.bam

############################## NCoR1_Uns_R2 #######################################################################
hisat2 -t -p 20  -x /Toolbox/aln_index/mm10_hisat2_Refseq/mm10 -1  NCoR1_Uns_R2_1.fastq.gz -2  NCoR1_Uns_R2_2.fastq.gz | samtools view -bS - >hisat2_out/NCoR1_Uns_R2.bam
samtools sort -@ 10 hisat2_out/NCoR1_Uns_R2.bam -o hisat2_out/NCoR1_Uns_R2.srt.bam
samtools index hisat2_out/NCoR1_Uns_R2.srt.bam
