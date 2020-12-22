library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)

# read raw count data from featureCount
count <- read.csv("data/SMRT_all_condition_count.txt",sep = "\t",header = T,row.names = 1)
head(count)
#######################################################################################################################
condition <- factor(rep(c(rep("0h",2),rep("2h",2),rep("6h",2)),2))
genotype <- factor(c(rep("WT",6),rep("KO",6)),levels = c("WT","KO"))
coldata <- data.frame(row.names=colnames(count),condition,genotype)
coldata$group<- paste0(coldata$condition,coldata$genotype)
#design <- ~ condition + genotype + condition:genotype
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~ group)
dds <- DESeq(dds)
resultsNames(dds)

################################################################################################################################
SMRT.volcano.plot.list = list()
SMRT_DE.list = list()

# 0hr
SU_vs_EU <- results(dds, contrast=c("group","0hKO","0hWT")) %>% as.data.frame() %>% rownames_to_column("Gene")
SMRT_DE.list[[1]] = SU_vs_EU %>% filter(log2FoldChange >=1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list)[1]= "SU_vs_EU_up"
SMRT_DE.list[[2]] = SU_vs_EU %>% filter(log2FoldChange <=-1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list)[2]= "SU_vs_EU_down"
SU_vs_EU %>% filter(abs(log2FoldChange) >=1 & padj <= 0.05) %>%
write.table(., file="results/SMRT_0h_KD_vs_ctrl_0h_DE.tsv", sep = "\t",quote = FALSE,row.names = FALSE)
SMRT.volcano.plot.list[[1]] = volcano.plot(SU_vs_EU,"SMRT KD vs Control (0hr)")
################################################################################################################################
# 2hr
SC_vs_EC.2 <- results(dds,contrast=c("group","2hKO","2hWT")) %>% as.data.frame() %>% rownames_to_column("Gene")
SMRT_DE.list[[3]] = SC_vs_EC.2 %>% filter(log2FoldChange >=1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list)[3]= "SC_vs_EC.2_up"
SMRT_DE.list[[4]] = SC_vs_EC.2 %>% filter(log2FoldChange <=-1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list)[4]= "SC_vs_EC.2_down"
SC_vs_EC.2 %>% filter(abs(log2FoldChange) >=1 & padj <= 0.05) %>%
write.table(., file="results/SMRT_2h_KD_vs_ctrl_2h_DE.tsv", sep = "\t",quote = FALSE,row.names = FALSE)
SMRT.volcano.plot.list[[2]] = volcano.plot(SC_vs_EC.2,"SMRT KD vs Control (2hr)")
################################################################################################################################
################################################################################################################################
# 6hr
SC_vs_EC.6 <- results(dds,contrast=c("group","6hKO","6hWT")) %>% as.data.frame() %>% rownames_to_column("Gene")
SMRT_DE.list[[5]] = SC_vs_EC.6 %>% filter(log2FoldChange >=1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list)[5]= "SC_vs_EC.6_up"
SMRT_DE.list[[6]] = SC_vs_EC.6 %>% filter(log2FoldChange <=-1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(SMRT_DE.list)[6]= "SC_vs_EC.6_down"
SC_vs_EC.6 %>% filter(abs(log2FoldChange) >=1 & padj <= 0.05) %>%
  write.table(., file="results/SMRT_6h_KD_vs_ctrl_6h_DE.tsv", sep = "\t",quote = FALSE,row.names = FALSE)
SMRT.volcano.plot.list[[3]] = volcano.plot(SC_vs_EC.6,"SMRT KD vs Control (6hr)")
################################################################################################################################

pdf("Figures/SMRT_KD_vs_control_volcano.pdf",height = 5,width = 15)
do.call(gridExtra::grid.arrange,c(SMRT.volcano.plot.list,ncol=3))
dev.off()
################################################################################################################################

rld  <- rlogTransformation(dds, blind=TRUE)
write.table(as.data.frame(assay(rld)),file="results/SMRT_all_condition_rld.txt",sep="\t",quote = FALSE)

######################## PCA #########################################################################
col <- read.csv("results/SMRT_all_condition_rld.txt",header = T,sep = "\t") %>% filter(rowSums(.[,]) >1)
pca_data=prcomp(t(col[,c(1:ncol(col))]), center=TRUE, scale=TRUE)
pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
plot(pca_data, type = "l")
df_pca_data=data.frame(PC1 = pca_data$x[,1],
                       PC2 = pca_data$x[,2], sample = colnames(col[,c(1:ncol(col))]))
pdf("Figures/SMRT_sample_PCA_plot.pdf",width = 8,height = 6)
pca_plt <- ggplot(df_pca_data, aes(PC1,PC2,colour =sample,label=sample))+
  geom_point(size=3)+
  geom_text_repel(
    data          = df_pca_data,
    size          = 4,
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

######################################################################################################
#                             Scatter plot                                                           #
######################################################################################################
volcano.plot = function(smrt,title){
  smrt <- smrt[which(smrt$log2FoldChange != "NA" & smrt$padj != "NA"),]
  smrt <- smrt %>% mutate(reg = case_when(
                    smrt$log2FoldChange >= 1 & smrt$padj < 0.05 ~ "UP",
                    smrt$log2FoldChange <= -1 & smrt$padj < 0.05 ~ "DOWN",
                    abs(smrt$log2FoldChange) < 1 & smrt$padj >= 0.05 ~ "no_change",
                    abs(smrt$log2FoldChange) < 1 & smrt$padj <= 0.05 ~ "no_change",
                    abs(smrt$log2FoldChange) > 1 & smrt$padj >0.05 ~ "no_change")) %>%
  mutate(reg = factor(reg, 
                      levels = c("UP", "no_change","DOWN")))
  label <- c("Il10","Il27","Ido1","Ido2","Cd274","Ctla4","Il6","Il12b",
           "Lag3", "Pparg", "Ifnb1", "Myd88", "Akt3", "Lag3","Cd83","Il12a","Ncor2")

  up = dim(smrt[which(smrt$reg =="UP"),])[1]
  down = dim(smrt[which(smrt$reg =="DOWN"),])[1]
  data= subset(smrt, Gene %in% label)
  data= data[which(abs(data$log2FoldChange) >=1),]
  data= data[which(abs(data$padj) <=0.05),]
  plt_smrt <- ggplot(smrt,aes(x=log2FoldChange,y=-log10(padj),label = Gene))+
                geom_point(aes(color=reg))+
                #xlim(-10,10)+
                #ylim(0,50)+
                scale_color_manual(name = "Differential \n regulation",
                                   values = c("DOWN" = "#058983",
                                              "no_change" = "grey",
                                              "UP" = "#DEB132"),
                                   labels = c("UP", "No Change", "DOWN"))+
                theme_bw()+
                xlab("log2 (Fold Change)")+ylab("-log10(adj p-value)")+
                theme(axis.text.x=element_text(size=16,color="black")
                      , axis.title.x=element_text(size=16)
                      , axis.text.y=element_text(size=16,color="black")
                      , axis.title.y=element_text(size=16)
                      ,legend.justification=c(0,1)
                      ,legend.box.background = element_rect(colour = "black")
                      ,plot.title = element_text(hjust = 0.5))+
                guides(colour = guide_legend(override.aes = list(size=6)))+
                geom_text_repel(
                  data          = data,
                  size          = 5,
                  direction    = "x",
                  angle        = 0,
                  vjust        = 0,
                  nudge_y       = 5 + data$log2FoldChange ,
                  segment.size  = 0.5,
                  segment.color = "black"
                ) +
                geom_vline(xintercept=c(-1,1), linetype="dashed",size=0.5)+
                geom_hline(yintercept=c(1.3), linetype="dashed",size=0.5)
  plt_smrt <- plt_smrt+annotate("text", x = -7, y = 200, label = down, color="#058983",size=8)
  plt_smrt <- plt_smrt+annotate("text", x = 7, y = 200,  label = up, color="#DEB132",size=8) +ggtitle(title) 
  return(plt_smrt)
}
################################################################################################################

