library(ChIPQC)
## Load sample data
samples <- data.frame(SampleID=c("SMRT_0hr","SMRT_6hr_CpG"),
                      Replicate=c(1,1),
                      Tissue = c("DCs","DCs"),
                      Factor = c("SMRT","SMRT"),
                      Condition =  c("0hr","6hr CpG"),
                      bamReads = c("/home/imgsb/Gyan/NCOR1/Analysis_result_of_NCoR1_ChIP_sep2017/bowtie_out/SMRT_0hr.dup.filtered.srt.q10.bam",
                                  "/home/imgsb/Gyan/NCOR1/NCoR1_SMRT_analysis/SMRT_CpG_9M.bam"),
                      ControlID = c("Input",
                                    "Input"), 
                      bamControl = c("/home/imgsb/Gyan/NCOR1/NCoR1_SMRT_analysis/Input_9M.bam",
                                    "/home/imgsb/Gyan/NCOR1/NCoR1_SMRT_analysis/Input_9M.bam"),
                      Peaks = c("/home/imgsb/Gyan/NCOR1/NCoR1_SMRT_analysis/tag_dir/peak_file/SMRT_0hr_9M.peaks_wdot_input.noSV.bed",
                                "/home/imgsb/Gyan/NCOR1/NCoR1_SMRT_analysis/tag_dir/peak_file/SMRT_CpG_9M.peaks_wdot_input.noSV.bed"),
                      PeakCaller = c("bed", "bed"))

## Create ChIPQC object
chipObj <- ChIPQC(samples, annotation="mm10",blacklist = "/home/imgsb/Gyan/mm10/mm10_blacklisted_region.bed",chromosomes = NULL) 
## Create ChIPQC report
#dir.create("results/ChIPQC")
ChIPQCreport(chipObj, reportName="ChIP QC report: SMRT", reportFolder="results/ChIPQC/")

samples <- data.frame(SampleID=c("NCoR1_0hr","NCoR1_6hr_CpG"),
                      Replicate=c(1,1),
                      Tissue = c("DCs","DCs"),
                      Factor = c("NCoR1","NCoR1"),
                      Condition =  c("0hr","6hr CpG"),
                      bamReads = c("/home/imgsb/Gyan/NCOR1/NCoR1_SMRT_analysis/NCoR1_Uns_9M.bam",
                                   "/home/imgsb/Gyan/NCOR1/NCoR1_SMRT_analysis/NCoR1_CpG_9M.bam"),
                      ControlID = c("Input",
                                    "Input"), 
                      bamControl = c("/home/imgsb/Gyan/NCOR1/NCoR1_SMRT_analysis/Input_9M.bam",
                                     "/home/imgsb/Gyan/NCOR1/NCoR1_SMRT_analysis/Input_9M.bam"),
                      Peaks = c("/home/imgsb/Gyan/NCOR1/NCoR1_SMRT_analysis/tag_dir/peak_file/NCoR1_Uns_9M.bam_peaks.sorted.bed",
                                "/home/imgsb/Gyan/NCOR1/NCoR1_SMRT_analysis/tag_dir/peak_file/NCoR1_CpG_9M.bam_peaks.sorted.bed"),
                      PeakCaller = c("bed", "bed"))

## Create ChIPQC object
chipObj <- ChIPQC(samples, annotation="mm10",blacklist = "/home/imgsb/Gyan/mm10/mm10_blacklisted_region.bed",chromosomes = NULL) 
## Create ChIPQC report
dir.create("results/NCoR1_ChIPQC")
ChIPQCreport(chipObj, reportName="ChIP QC report: SMRT", reportFolder="results/NCoR1_ChIPQC/")