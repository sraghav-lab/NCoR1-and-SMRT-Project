library(DESeq2)

#count data
setwd("/home/imgsb/Gyan/NCOR1/Analysis_result_of_NCoR1_CpG_RNAseq_data/")
count <- read.csv("/home/imgsb/Gyan/NCOR1/Analysis_result_of_NCoR1_CpG_RNAseq_data/NCoR1_SMRT_CPG_all_condition_count.txt",sep = "\t",header = T)
head(count)
rownames(count) <- count$Geneid
#count1 <- count[,c(2:ncol(count))]
count1 <- count[,c(2:13,16,17,14,15,20,21,18,19)]
#Remove rows if count is <  in 90% of sample
# rem <- function(x){
#   x <- as.matrix(x)
#   x <- t(apply(x,1,as.numeric))
#   r <- as.numeric(apply(x,1,function(i) sum(i == 0) ))
#   remove <- which(r > dim(x)[2]*0.1)
#   return(remove)
# }
# remove <- rem(count1)
# countdata <- count1[-remove,]
#######################################################################################################################
countdata <- count1[rowSums(count1)>1,]
head(countdata)
condition <- factor(c(rep(c(rep("0h",2),rep("2h",2),rep("6h",2)),2),rep(c(rep("0h",2),rep("6h",2)),2)))
genotype <- factor(c(rep("WT",6),rep("S_KD",6),rep("WT","4"),rep("N_KD",4)))
Rep <- factor(c(rep(c(rep("R1",1),rep("R2",1)),10)))
batch <- factor(c(rep("1",12),rep("2",8)))
coldata <- data.frame(row.names=colnames(countdata),genotype,condition,batch)
coldata$group<- factor(paste0(coldata$condition,"_",coldata$genotype))
#design <- ~ condition + genotype + condition:genotype
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~group)
dds <- DESeq(dds)

######################################################################################################
#normalized count
setwd("/home/imgsb/Gyan/NCOR1/NCoR1_SMRT_analysis/")
norm_count <- as.data.frame(counts(dds, normalized=TRUE))[,]
write.table(norm_count, file="NCoR1_SMRT_all_condition_norm_count.txt", sep = "\t",quote = FALSE)

######################################################################################################
rld<- rlogTransformation(dds, blind=TRUE)
write.table(as.data.frame(assay(rld)),file="/home/imgsb/Gyan/NCOR1/NCoR1_SMRT_analysis/NCoR1_SMRT_all_condition_rld.txt",sep="\t",quote = FALSE)

#######################################################################################################
vst<-varianceStabilizingTransformation(dds, blind=TRUE)
write.table(as.data.frame(assay(vst)),"/home/imgsb/Gyan/NCOR1/NCoR1_SMRT_analysis/NCoR1_SMRT_all_condition_vst.txt",sep = "\t",quote = FALSE)
######################## PCA #########################################################################
col <- read.csv("/home/imgsb/Gyan/NCOR1/NCoR1_SMRT_analysis/NCoR1_SMRT_all_condition_rld.txt",header = T,sep = "\t",row.names = 1)
head(col)
pca_data=prcomp(t(col[,c(1:ncol(col))]), center=TRUE, scale=TRUE)
pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
plot(pca_data, type = "l")
df_pca_data=data.frame(PC1 = pca_data$x[,1],
                       PC2 = pca_data$x[,2], sample = colnames(col[,c(1:ncol(col))]))
png("../RNAseq_chipseq_integration/expression/NCoR1_PCA_plot.png",width = 12,height = 10,units = "cm",res=500 )
pca_plt <- ggplot(df_pca_data, aes(PC1,PC2,colour =sample,label=sample))+
  geom_point(size=3)+
  geom_text_repel(
    data          = df_pca_data,
    size          = 7,
    #angle         = 45,
    #fill = df_pca_data$sample,
    #arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "first"),
    #nudge_y       = 0.5,
    segment.size  = 0,
    segment.color = "black",
    direction     = "y"
  ) +
  #geom_text(aes(label = sample), size = 5, position=position_jitter(width=6, height=5))+#position = position_nudge(y = 4),check_overlap = T)+
  labs(x=paste0("PC1 (",pca_data_perc[1],"%",")"), y=paste0("PC2 (",pca_data_perc[2],"%",")"))+
  theme_bw() +
  theme(axis.text.x=element_text(size=rel(2), angle=90))+
  theme(axis.text.y=element_text(size=rel(2), angle=90),legend.position = "none")
pca_plt
dev.off()

#####################################################################################################
### PCA for NCoR1 and SMRT Data #####################################################################
library(ggplot2)
library(ggrepel)
col <- read.csv("/home/imgsb/Gyan/NCOR1/NCoR1_SMRT_analysis/WGCNA_with_more_DE_genes/NCoR1_SMRT_all_condition_DE_genes_1.5fold_list_rld.txt",header = T,sep = "\t",row.names = 1)
head(col)
col <- col[,c(1,2,13,14,17,18,7,8,5,6,15,16,19,20,11,12)]
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0) ))
  remove <- which(r > dim(x)[2]*0.2)
  return(remove)
}
remove <- rem(col)
col1 <- col[-remove,]
pca_data=prcomp(t(col[,c(1:ncol(col))]), center=TRUE, scale=FALSE)
pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
plot(pca_data, type = "l")
df_pca_data=data.frame(PC1 = pca_data$x[,1],
                       PC2 = pca_data$x[,2], sample = colnames(col[,c(1:ncol(col))]))
#png("../RNAseq_chipseq_integration/expression/NCoR1_PCA_plot.png",width = 12,height = 10,units = "cm",res=500 )
pca_plt <- ggplot(df_pca_data, aes(PC1,PC2,colour =batch))+
  geom_point(aes(shape=genotype,size=5))+
  geom_mark_ellipse(aes(fill = time_point)) +
  #geom_text_repel(
    #data          = df_pca_data,
    #size          = 7,
    #angle         = 45,
    #fill = df_pca_data$sample,
    #arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "first"),
    #nudge_y       = 0.5,
    #segment.size  = 0,
    #segment.color = "black",
    #direction     = "y"
 # ) +
  #geom_text(aes(label = sample), size = 5, position=position_jitter(width=6, height=5))+#position = position_nudge(y = 4),check_overlap = T)+
  labs(x=paste0("PC1 (",pca_data_perc[1],"%",")"), y=paste0("PC2 (",pca_data_perc[2],"%",")"))+
  theme_bw() +
  theme(axis.text.x=element_text(size=rel(2), angle=90))+
  theme(axis.text.y=element_text(size=rel(2), angle=90))
pca_plt
dev.off()
####### 3D PCA  ######################################################################################
library(scatterplot3d)
library()
pca_data1 <- as.data.frame(pca_data$x)
pca_data1$Sample <- rownames(pca_data1)
pca_data1$Sample <- factor(pca_data1$Sample)
pca_data1$sample <- factor(c(rep("Ctrl 0hr",4),rep("NCoR1 KD 0hr",2),rep("SMRT KD 0hr",2),rep("Ctrl 6hr",4),rep("NCoR1 KD 6hr",2),rep("SMRT KD 6hr",2)))
#pca_data1$sample <- factor(pca_data1$sample,levels = legend)
#colors <- c("blue","blue","blue","blue","black","black", "#058983", "#058983","#DEB132","#DEB132","#DEB132","#DEB132","#977411","#977411","red","red")
colors <- c("blue","#6E0070", "#058983","#DEB132","#977411","red")
colors <- colors[as.numeric(pca_data1$sample)]
layout(matrix(c(1,2), nrow = 1), widths = c(0.5, 0.5))
par(mar = c(5, 4, 4, 2) + 0.1)
###################
s3d <- scatterplot3d(pca_data1[,1:3], pch = 19,type="h",color=colors,
                     scale.y=1,lty.hplot=2,
                     xlab=paste("PC1, ", round(pca_data_perc[1], 2), "%"), 
                     ylab=paste("PC2, ", round(pca_data_perc[2], 2), "%"), 
                     zlab=paste("PC3, ", round(pca_data_perc[3], 2), "%"), 
                     grid=TRUE, box=FALSE)
#text(s3d$xyz.convert(pca_data1[, 1:3]), labels = rownames(pca_data1),cex= 0.7, col = "black",pos=2.5)
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
addgrids3d(pca_data1[,1:3], grid = c("xy", "xz", "yz"))
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("topleft",legend = levels(pca_data1$sample),
       col =  c("blue","#6E0070", "#058983","#DEB132","#977411","red"),
       pch = 16,xpd=FALSE, inset=c(-.0, -.0), cex=1.2)
######################################################################################################




vsd <- vst(dds, blind=TRUE)
test <- as.data.frame(assay(vsd))
plt <- plotPCA(vsd, intgroup=c("grp","cnd","batch"))
plt

mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat, batch = col_data$batch,design = count.dat.mod)
assay(vsd) <- mat
plotPCA(vsd,intgroup=c("grp","cnd","batch"))



##################################################################################################################
