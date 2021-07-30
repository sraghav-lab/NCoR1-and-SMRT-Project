# Feature Count

# list Alignment file (BAM) from Batch 1 RNASeq data
batch1.files = list.files(path = "/home/imgsb/Gyan/NCOR1/Analysis_result_of_SMRT_RNAseq_data/hisat2_out/",
                          pattern = "*.srt.bam$",full.names = TRUE)
#Remove 2hr data
batch1.files = batch1.files[-c(3,4,9,10)]

# list Alignment file (BAM) from Batch 2 RNASeq data
batch2.files = list.files(path = "/media/Hard_disk1/Analysis_backup/SMRT_RNASeq_July_2021/hisat2_out/",
                          pattern = "*.bam$",full.names = TRUE)


CountData.batch1 = featureCounts(files = batch1.files,
                          annot.ext = "/home/imgsb/Gyan/mm10/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf",
                          isGTFAnnotationFile = TRUE,
                          GTF.attrType ="gene_name",
                          minMQS = 30,verbose = TRUE,nthreads = 20,requireBothEndsMapped =TRUE,
                          isPairedEnd = TRUE,primaryOnly = TRUE)


CountData.batch2 = featureCounts(files = batch2.files,
                          annot.ext = "/home/imgsb/Gyan/mm10/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf",
                          isGTFAnnotationFile = TRUE,
                          GTF.attrType ="gene_name",
                          minMQS = 30,verbose = TRUE,nthreads = 20,requireBothEndsMapped =TRUE,
                          isPairedEnd = TRUE,primaryOnly = TRUE)

CountData = merge(CountData.batch1$counts,CountData.batch2$counts,by=0) %>% column_to_rownames("Row.names")
colnames(CountData) = gsub("Crtl","Ctrl",colnames(CountData))
#######################################################################################################################
Stimulation = gsub("Ctrl.|SMRT.|.R[0-9].srt.bam|.R[0-9].bam","",colnames(CountData))
genotype = gsub(".0h.R[0-9]|.CpG.6h.R[0-9]|.srt.bam|.bam","",colnames(CountData))
Batch = c(rep(1,8),rep(2,12))
MetaData <- data.frame(Sample=colnames(CountData),
                       Stimulation=factor(Stimulation),
                       genotype = factor(genotype),
                       batch=factor(Batch))
MetaData$Group = paste0(MetaData$Stimulation,"_",MetaData$genotype)

dds <- DESeqDataSetFromMatrix(countData=CountData[,c(9:20)], colData=MetaData[c(9:20),], design=~Group)
dds <- DESeq(dds)


Sample.vst <- vst(dds, blind=FALSE)
plt <- plotPCA(Sample.vst, intgroup=c("Group","batch"),returnData=TRUE)
percentVar <- round(100 * attr(plt, "percentVar"))
plt$genotype = MetaData[c(1:8),]$genotype
########################################################################################################
pdf("Figures/NCoR1_SMRT_sample_PCA.pdf",width = 8,height = 5)
ggplot(plt, aes(PC1, PC2, shape=genotype)) + 
  geom_point(size=8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
  #scale_color_manual(values = c("grey","#058983","#DEB132"))+
  gg_theme+
  annotate("rect", xmin=c(-30,16), xmax=c(-18,33), ymin=c(-22,-22) , ymax=c(22,22), alpha=0,size=1, color=c("blue","red"))+
  scale_color_manual(values = c("grey","#DEB132"),labels = c("Feb 2019", "Jul 2021"))+
  scale_shape_discrete(labels = c("Control", "SMRT KD")) 
dev.off()

########################################################################################################
SMRT.volcano.plot.list.batch.1_2 = list()
SMRT_DE.list.batch.1_2 = list()

# 0hr
SU_vs_EU.1_2 <- results(dds, contrast=c("Group","0h_SMRT","0h_Ctrl")) %>% as.data.frame() %>% rownames_to_column("Gene")
SMRT_DE.list.batch.1_2[[1]] = SU_vs_EU.1_2 %>% filter(log2FoldChange >=1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.batch.1_2)[1]= "SU_vs_EU_up"
SMRT_DE.list.batch.1_2[[2]] = SU_vs_EU.1_2 %>% filter(log2FoldChange <=-1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.batch.1_2)[2]= "SU_vs_EU_down"
#SU_vs_EU %>% filter(abs(log2FoldChange) >=1 & padj <= 0.05) %>%
#  write.table(., file="results/SMRT_0h_KD_vs_ctrl_0h_DE.tsv", sep = "\t",quote = FALSE,row.names = FALSE)
SMRT.volcano.plot.list.batch.1_2[[1]] = volcano.plot(SU_vs_EU.1_2,"SMRT KD vs Control (0hr)")
################################################################################################################################

#############################################################################################################################
# 6hr
SC_vs_EC.6.1_2 <- results(dds,contrast=c("Group","CpG.6h_SMRT","CpG.6h_Ctrl")) %>% as.data.frame() %>% rownames_to_column("Gene")
SMRT_DE.list.batch.1_2[[3]] = SC_vs_EC.6.1_2 %>% filter(log2FoldChange >=1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.batch.1_2)[3]= "SC_vs_EC.6_up"
SMRT_DE.list.batch.1_2[[4]] = SC_vs_EC.6.1_2 %>% filter(log2FoldChange <=-1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.batch.1_2)[4]= "SC_vs_EC.6_down"
#SC_vs_EC.6.1_2 %>% filter(abs(log2FoldChange) >=1 & padj <= 0.05) %>%
#                   write.table(., file="results/SMRT_6h_KD_vs_ctrl_6h_DE.tsv", sep = "\t",quote = FALSE,row.names = FALSE)
SMRT.volcano.plot.list.batch.1_2[[3]] = volcano.plot(SC_vs_EC.6.1_2,"SMRT KD vs Control (6hr)")

