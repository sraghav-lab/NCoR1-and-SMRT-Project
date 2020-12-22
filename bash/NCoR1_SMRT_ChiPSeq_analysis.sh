#!/usr/bin/bash

#normalize NCoR1 Chipseq (0hr and 6hr CpG to SMRT depth)
#normalize Input Chipseq to lowest no. o sequencing depth)

#sequencing depth
#Input = 53156118 
#NCoR1 0hr = 35464501
#NCoR1 6hr CpG = 35465372
#SMRT 0hr = 9061286
#SMRT 6hr CpG = 9742021


#PICCARD
#java -jar /apps/picard.jar DownsampleSam INPUT=/home/imgsb/results/Gyan/NCOR1/Analysis_result_of_CHIP_seq_using_homer/NCOR1_mm10_alignment/Uniquley_mapped_reads/sorted_NCoR1_Ch_Input_R1.bam OUTPUT=./Input_9M.bam PROBABILITY=0.170
#java -jar /apps/picard.jar DownsampleSam INPUT=/home/imgsb/results/Gyan/NCOR1/Analysis_result_of_CHIP_seq_using_homer/normalized_ChIPseq_bam_file_to_input/NCoR1_CpG_R1.bam OUTPUT=./NCoR1_CpG_9M.bam PROBABILITY=0.255
#java -jar /apps/picard.jar DownsampleSam INPUT=/home/imgsb/results/Gyan/NCOR1/Analysis_result_of_CHIP_seq_using_homer/normalized_ChIPseq_bam_file_to_input/NCoR1_Uns_R1.bam OUTPUT=./NCoR1_Uns_9M.bam PROBABILITY=0.255
#java -jar /apps/picard.jar DownsampleSam INPUT=/home/imgsb/results/Gyan/NCOR1/Analysis_result_of_NCoR1_ChIP_sep2017/bowtie_out/SMRT_CpG.dup.filtered.srt.q10.bam OUTPUT=./SMRT_CpG_9M.bam PROBABILITY=0.930


#Run make tag directory
mkdir tag_dir
cat bam_file.txt| parallel --verbose 'makeTagDirectory tag_dir{} {}' ::: *9M.bam
wait
cat bam_file.txt |parallel --verbose 'makeUCSCfile tag_dir{} -o auto -bigWig /home/imgsb/results/Gyan/mm10/mm10_size' ::: *9M.bam
wait
head -4 bam_file.txt |parallel --verbose 'findPeaks tag_dir/{} -style factor -i tag_dir/Input_9M.bam -o auto' 
wait
head -4 bam_file.txt |parallel --verbose 'pos2bed.pl tag_dir/{}/peaks.txt > tag_dir/{}/peaks.bed' 
wait
head -4 bam_file.txt |parallel --verbose 'cat tag_dir/{}/peaks.bed | grep -v random > tag_dir/{}/{}.peaks.bed'
wait
head -4 bam_file.txt |parallel --verbose 'sort-bed tag_dir/{}/{}.peaks.bed > tag_dir/{}/{}.peaks.srt.bed' 
wait
head -4 bam_file.txt |parallel --verbose 'cp tag_dir/{}/{}.peaks.bed .' 
wait
head -4 bam_file.txt |parallel --verbose 'cp tag_dir/{}/{}.ucsc.bigWig .'
