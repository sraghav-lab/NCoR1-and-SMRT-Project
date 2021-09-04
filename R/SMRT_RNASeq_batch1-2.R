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

saveRDS(CountData,file = "results/SMRT_batch1-2_RNASeq_RawCount.rds")
#######################################################################################################################
# Prepare Metadata
Stimulation = gsub("Ctrl.|SMRT.|.R[0-9].srt.bam|.R[0-9].bam","",colnames(CountData))
genotype = gsub(".0h.R[0-9]|.CpG.6h.R[0-9]|.srt.bam|.bam","",colnames(CountData))
Batch = c(rep(1,8),rep(2,12))
MetaData <- data.frame(Sample=colnames(CountData),
                       Stimulation=factor(Stimulation),
                       genotype = factor(genotype),
                       batch=factor(Batch))
MetaData$Group = paste0(MetaData$Stimulation,"_",MetaData$genotype)
MetaData$Replicate = Replicate = c(rep(c("Rep 1","Rep 2"),4),rep(c("Rep 3","Rep 4" ,"Rep 5"),4))

dds.5 <- DESeqDataSetFromMatrix(countData=CountData, colData=MetaData, design=~Group+batch)
dds.5 <- DESeq(dds.5)


SMRT_Sample.vst.Rep5 <- vst(dds.5, blind=FALSE)
saveRDS(SMRT_Sample.vst.Rep5,file = "results/SMRT_batch1-2_RNASeq_vst.rds")
plt <- plotPCA(Sample.vst, intgroup=c("Group","batch"),returnData=TRUE,ntop =1000)
percentVar <- round(100 * attr(plt, "percentVar"))
plt$genotype = MetaData$genotype
########################################################################################################
pdf("Figures/SMRT_sample_PCA.pdf",width = 8,height = 5)
ggplot(plt, aes(PC1, PC2, color = batch,shape=genotype)) + 
       geom_point(size=8) +
       xlab(paste0("PC1: ",percentVar[1],"% variance")) +
       ylab(paste0("PC2: ",percentVar[2],"% variance")) +
       #scale_color_manual(values = c("grey","#058983","#DEB132"))+
       gg_theme+
       annotate("rect", xmin=c(-37,16), xmax=c(-23,42), ymin=c(-30,-30) , ymax=c(28,28), alpha=0,size=1, color=c("blue","red"))+
       scale_color_manual(values = c("grey","#DEB132"),labels = c("Old Batch", "New Batch"))+
       scale_shape_discrete(labels = c("Control", "SMRT KD")) 
dev.off()

########################################################################################################
SMRT.volcano.plot.list.5Rep = list()
SMRT_DE.list.5Rep = list()
set.seed(100)
# 0hr
SU_vs_EU.5Rep <- results(dds.5, contrast=c("Group","0h_SMRT","0h_Ctrl")) %>% as.data.frame() %>% rownames_to_column("Gene")
SMRT_DE.list.5Rep[[1]] = SU_vs_EU.5Rep %>% filter(log2FoldChange >=1 & padj <= 0.01) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.5Rep)[1]= "SU_vs_EU_up"
SMRT_DE.list.5Rep[[2]] = SU_vs_EU.5Rep %>% filter(log2FoldChange <=-1 & padj <= 0.01) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.5Rep)[2]= "SU_vs_EU_down"
SMRT.volcano.plot.list.5Rep[[1]] = volcano.plot(SU_vs_EU.5Rep,"SMRT KD vs Control (0hr)")
write.xlsx(x=SU_vs_EU.5Rep %>% filter(abs(log2FoldChange) >=1 & padj <= 0.01),
           file="results/SMRT_DEGs_DESeq2.xlsx",sheetName ="SMRT KD vs Control(0hr)", append=T)
# 6hr
SC_vs_EC.6.5Rep <- results(dds.5,contrast=c("Group","CpG.6h_SMRT","CpG.6h_Ctrl")) %>% as.data.frame() %>% rownames_to_column("Gene")
SMRT_DE.list.5Rep[[3]] = SC_vs_EC.6.5Rep %>% filter(log2FoldChange >=1 & padj <= 0.01) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.5Rep)[3]= "SC_vs_EC.6_up"
SMRT_DE.list.5Rep[[4]] = SC_vs_EC.6.5Rep %>% filter(log2FoldChange <=-1 & padj <= 0.01) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.5Rep)[4]= "SC_vs_EC.6_down"
SMRT.volcano.plot.list.5Rep[[2]] = volcano.plot(SC_vs_EC.6.5Rep,"SMRT KD vs Control (6hr CpG)")
write.xlsx(x=SC_vs_EC.6.5Rep %>% filter(abs(log2FoldChange) >=1 & padj <= 0.01),
           file="results/SMRT_DEGs_DESeq2.xlsx",sheetName ="SMRT KD vs Control(6hr CpG)", append=T)

# 0hr
SMRT_DE.list.5Rep[[5]] = SU_vs_EU.5Rep %>% filter(log2FoldChange >=1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.5Rep)[5]= "SU_vs_EU_up"
SMRT_DE.list.5Rep[[6]] = SU_vs_EU.5Rep %>% filter(log2FoldChange <=-1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.5Rep)[6]= "SU_vs_EU_down"

# 6hr
SMRT_DE.list.5Rep[[7]] = SC_vs_EC.6.5Rep %>% filter(log2FoldChange >=1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.5Rep)[7]= "SC_vs_EC.6_up"
SMRT_DE.list.5Rep[[8]] = SC_vs_EC.6.5Rep %>% filter(log2FoldChange <=-1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.5Rep)[8]= "SC_vs_EC.6_down"

pdf("Figures/SMRT_KD_vs_control_volcano.pdf",height = 5,width = 10)
do.call(gridExtra::grid.arrange,c(SMRT.volcano.plot.list.5Rep,ncol=2))
dev.off()



################################################################################################################################
#############################################################################################################################

dds.3 <- DESeqDataSetFromMatrix(countData=CountData[,c(9:20)], colData=MetaData[c(9:20),], design=~Group)
dds.3 <- DESeq(dds.3)

SMRT.volcano.plot.list.3Rep = list()
SMRT_DE.list.3Rep = list()

# 0hr
SU_vs_EU.3Rep <- results(dds.3, contrast=c("Group","0h_SMRT","0h_Ctrl")) %>% as.data.frame() %>% rownames_to_column("Gene")
SMRT_DE.list.3Rep[[1]] = SU_vs_EU.3Rep %>% filter(log2FoldChange >=1 & padj <= 0.01) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.3Rep)[1]= "SU_vs_EU_up"
SMRT_DE.list.3Rep[[2]] = SU_vs_EU.3Rep %>% filter(log2FoldChange <=-1 & padj <= 0.01) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.3Rep)[2]= "SU_vs_EU_down"

# 6hr
SC_vs_EC.6.3Rep <- results(dds.3,contrast=c("Group","CpG.6h_SMRT","CpG.6h_Ctrl")) %>% as.data.frame() %>% rownames_to_column("Gene")
SMRT_DE.list.3Rep[[3]] = SC_vs_EC.6.3Rep %>% filter(log2FoldChange >=1 & padj <= 0.01) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.3Rep)[3]= "SC_vs_EC.6_up"
SMRT_DE.list.3Rep[[4]] = SC_vs_EC.6.3Rep %>% filter(log2FoldChange <=-1 & padj <= 0.01) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.3Rep)[4]= "SC_vs_EC.6_down"

# 0hr
SMRT_DE.list.3Rep[[5]] = SU_vs_EU.3Rep %>% filter(log2FoldChange >=1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.3Rep)[5]= "SU_vs_EU_up"
SMRT_DE.list.3Rep[[6]] = SU_vs_EU.3Rep %>% filter(log2FoldChange <=-1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.3Rep)[6]= "SU_vs_EU_down"

# 6hr
SMRT_DE.list.3Rep[[7]] = SC_vs_EC.6.3Rep %>% filter(log2FoldChange >=1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.3Rep)[7]= "SC_vs_EC.6_up"
SMRT_DE.list.3Rep[[8]] = SC_vs_EC.6.3Rep %>% filter(log2FoldChange <=-1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list.3Rep)[8]= "SC_vs_EC.6_down"


####################################################################################################################
SMRT_Sample_Dist = assay(SMRT_Sample.vst.Rep5) %>% 
                        as.data.frame() %>% 
                        filter(rownames(.) %in% unique(unlist(SMRT_DE.list.5Rep[c(1,2,3,4)],use.names = FALSE))) %>%
                        t() %>%
                        dist() 
SMRT_Sample_Dist_mat<- as.matrix(SMRT_Sample_Dist)
rownames(SMRT_Sample_Dist_mat) <- colnames(SMRT_Sample_Dist_mat) <- with(colData(dds.5),paste0(Group, " ", Replicate,":",batch))
my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 50)
SMRT_Sample_Dist.map = Heatmap(SMRT_Sample_Dist_mat,col = my_palette,name = "Euclidean Distance",
                               column_title="Euclidean distance",
                               cell_fun = function(j, i, x, y, width, height, fill) {
                                       grid.text(sprintf("%.0f", SMRT_Sample_Dist_mat[i, j]), x, y, gp = gpar(fontsize = 10))
                               })

# All genes 
SMRT_Sample_Cor.all.genes = assay(SMRT_Sample.vst.Rep5) %>% 
                                        as.data.frame() %>%
                                        t() %>% scale() %>% t() %>% as.data.frame() %>%
                                        filter(Ctrl.0h.R4.bam != "NaN") %>%
                                        cor(method = "pearson") 

rownames(SMRT_Sample_Cor.all.genes) <- colnames(SMRT_Sample_Cor.all.genes) <- with(colData(dds.5),paste0(Group, " ", Replicate,":",batch))
my_palette <- colorRampPalette(c( "blue", "white","red"))(100)

SMRT_Sample_all.genes.Cor.map = Heatmap(SMRT_Sample_Cor.all.genes,col = my_palette,
                                        name = "Pearson\nCorrelation",
                                        column_title="Total Genes(24421)",
                                        cell_fun = function(j, i, x, y, width, height, fill) {
                                        grid.text(sprintf("%.2f", SMRT_Sample_Cor.all.genes[i, j]), 
                                                 x, y, gp = gpar(fontsize = 10))})

# All DE genes 
SMRT_Sample_Cor.DE.genes = assay(SMRT_Sample.vst.Rep5) %>% 
                                as.data.frame() %>% #head()
                                filter(rownames(.) %in% unique(unlist(SMRT_DE.list.5Rep[c(1,2,3,4)],use.names = FALSE))) %>% 
                                t() %>% scale() %>% t() %>% 
                                cor(method = "pearson") 

rownames(SMRT_Sample_Cor.DE.genes) <- colnames(SMRT_Sample_Cor.DE.genes) <- with(colData(dds.5),paste0(Group, " ", Replicate,":",batch))
my_palette <- colorRampPalette(c( "blue", "white","red"))(100)
SMRT_Sample_Cor.DE.genes.map = Heatmap(SMRT_Sample_Cor.DE.genes,col = my_palette,
                                        name = "Pearson\nCorrelation",
                                        column_title="Differentially expressed genes(2804)",
                                        cell_fun = function(j, i, x, y, width, height, fill) {
                                                grid.text(sprintf("%.2f", SMRT_Sample_Cor.DE.genes[i, j]), x, y, gp = gpar(fontsize = 10))
                                        })
# Top 500 genes
SMRT_Sample.vst.Rep5.Var = assay(SMRT_Sample.vst.Rep5) %>% 
                                as.data.frame() %>% #rownames_to_column("Genes") %>%
                                mutate(Var = rowVars(as.matrix(.))) %>% 
                                arrange(desc(Var)) %>% 
                                top_n(1000)
SMRT_Sample_Cor = SMRT_Sample.vst.Rep5.Var %>% dplyr::select(-Var) %>% 
                        t() %>% scale() %>% t() %>% 
                        cor(method = "pearson") 

rownames(SMRT_Sample_Cor) <- colnames(SMRT_Sample_Cor) <- with(colData(dds.5),paste0(Group, " ", Replicate,":",batch))
my_palette <- colorRampPalette(c( "blue", "white","red"))(100)
SMRT_Sample_Cor.top500.map = Heatmap(SMRT_Sample_Cor,col = my_palette,
                              name = "Pearson\nCorrelation",
                              column_title="Top 500 variable genes",
                              cell_fun = function(j, i, x, y, w, h, fill) {
                                      #grid.text(sprintf("%.2f", SMRT_Sample_Cor[i, j]), x, y, gp = gpar(fontsize = 10))
                                      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                              })
pdf("Figures/SMRT_sample_Correlation_plot.pdf",height = 10,width = 10)
corrplot(SMRT_Sample_Cor,method ="square", type="upper", order="hclust", col=my_palette,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45,#Text label color and rotation
         number.cex = 0.7,
         number.font = 1)
dev.off()

pdf("Figures/SMRT_sample_Correlation_plot.pdf",height = 10,width = 30)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 3)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(SMRT_Sample_all.genes.Cor.map, column_title = "", show_heatmap_legend = TRUE, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
draw(SMRT_Sample_Cor.DE.genes.map , column_title = "", show_heatmap_legend = TRUE, newpage = FALSE)
upViewport()


pushViewport(viewport(layout.pos.row=1, layout.pos.col=3))
draw(SMRT_Sample_Cor.top500.map , column_title = "", show_heatmap_legend = TRUE, newpage = FALSE)
upViewport()

upViewport()
dev.off()


