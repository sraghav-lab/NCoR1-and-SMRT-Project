# Program:featureCounts v1.6.2
# Command:
featureCounts --verbose -p -B -T 20 
          -a mm10_UCSC_genes.gtf 
          -o NCoR1_SMRT_CPG_all_condition_count.txt 
          Crtl_0h_R2.srt.bam 
          Crtl_0h_R3.srt.bam 
          Ctrl_CpG_2h_R2.srt.bam 
          Ctrl_CpG_2h_R3.srt.bam 
          Ctrl_CpG_6h_R2.srt.bam 
          Ctrl_CpG_6h_R3.srt.bam 
          SMRT_0h_R2.srt.bam 
          SMRT_0h_R3.srt.bam 
          SMRT_CpG_2h_R2.srt.bam 
          SMRT_CpG_2h_R3.srt.bam 
          SMRT_CpG_6h_R2.srt.bam 
          SMRT_CpG_6h_R3.srt.bam 
          Emp_CpG_R2.srt.bam 
          Emp_Uns_R2.srt.bam 
          NCoR1_CpG_R2.srt.bam 
          NCoR1_Uns_R2.srt.bam
