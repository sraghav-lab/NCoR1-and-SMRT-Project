#count data
count <- read.csv("data/NCoR1_SMRT_CPG_all_condition_count.txt",sep = "\t",header = T,row.names = 1)
count <- count[,c(15,16,13,14,19,20,17,18)]


count_dat = merge(count,CountData,by=0) %>% column_to_rownames("Row.names")
saveRDS(count_dat, file = "results/NCoR1_SMRT_RNASeq_RawCount.rds")

group <- c(rep("NCoR1_Ctrl",4),rep("NCoR1_KD",4),rep("SMRT_Ctrl",4),rep("SMRT_KD",4),rep("SMRT_Ctrl",6),rep("SMRT_KD",6))
cond <- c(rep(c(rep("0hr",2),rep("6hr_CpG",2)),4),rep(c(rep("0hr",3),rep("6hr_CpG",3)),2))
batch <- c(rep(1,8),rep(2,8),rep(3,12))
grp_cnd <- factor(paste0(group,"_",cond,"_",batch))
col_data <- data.frame(Sample=colnames(count_dat),
                       grp=factor(group),
                       #ind=factor(ind),
                       cnd=factor(cond),
                       group=factor(grp_cnd),
                       batch=factor(batch))

dds.all <- DESeqDataSetFromMatrix(countData=count_dat, colData=col_data, design=~group)
dds.all <- DESeq(dds.all)


vst.all <- vst(dds.all, blind=FALSE)
saveRDS(vst.all, file = "results/NCoR1_SMRT_RNASeq_vst.rds")
plt <- plotPCA(vst.all,intgroup=c("grp","cnd","batch"),returnData=TRUE)
percentVar <- round(100 * attr(plt, "percentVar"))
colnames(plt)[c(4,5,6)] = c("Group","Stimulation","Batch")
plt$grp = col_data$group
plt$Genotype = gsub("_.*","",plt$Group)
########################################################################################################
pdf("Figures/NCoR1_SMRT_sample_PCA1.pdf",width = 8,height = 5)
ggplot(plt, aes(PC1, PC2, color=Group,shape=Genotype)) + 
  geom_point(size=8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c("grey","#DEB132","grey","#058983"))+
  gg_theme+
  annotate("rect", xmin=c(-44,14), xmax=c(-22,55), ymin=c(-25,-25) , ymax=c(20,20), alpha=0,size=1, color=c("blue","red"))
dev.off()
######################################################################################################
#normalized count
NCoR1_SMRT_norm_count <- as.data.frame(counts(dds.all, normalized=TRUE))[,] %>% rownames_to_column("Gene")
saveRDS(NCoR1_SMRT_norm_count, file = "results/NCoR1_SMRT_RNASeq_normalizeCount.rds")

######################################################################################################
# Average of log of variance stabilized normalized count 
colnames(vst.all)
NCoR1_SMRT_vst_rep_merged <- as.data.frame(assay(vst.all)) %>% 
                                rownames_to_column("Gene") %>% 
                                melt() %>% mutate(variable=gsub("_R[0-9]|.R[0-9].srt.bam|.R[0-9].bam","",variable)) %>% 
                                dplyr::group_by(Gene,variable) %>% 
                                dplyr::summarise(Avg =mean(value)) %>% 
                                dcast(Gene~variable,value.var = "Avg")
write.table(NCoR1_SMRT_vst_rep_merged,file="results/NCoR1_SMRT_vst_rep_merged.tsv",sep="\t",quote = FALSE)


NCoR1_SMRT_rld_rep_merged.scaled <- scale(NCoR1_SMRT_rld_rep_merged,center = TRUE) %>% as.data.frame()
#######################################################################################################
# Since there were variability between SMRT Control and NCoR1 control sample due to batch effect we did 
# differential expression analyses for NCoR1 sample with its matched control

#count data
count <- read.csv("data/NCoR1_SMRT_CPG_all_condition_count.txt",sep = "\t",header = T,row.names = 1)
count <- count[,c(15,16,13,14,19,20,17,18)]

#######################################################################################################################
coldata <- data.frame(row.names=colnames(count),Conditon = gsub("_R[0-9]","",colnames(count)))
dds.control <- DESeqDataSetFromMatrix(countData=count, colData=coldata, design=~Conditon)
dds.control <- DESeq(dds.control)
#######################################################################################################################

NCoR1.volcano.plot.list = list()
NCoR1_DE.list = list()
Emp_DE.list = list()

# Control 6hr CpG vs Uns (0hr)
EC_vs_EU <- results(dds.control, contrast=c("Conditon","Emp_CpG_6h","Emp_Uns")) %>% as.data.frame() %>%  rownames_to_column('Gene')
Emp_DE.list[[1]] = EC_vs_EU %>% filter(log2FoldChange >=1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(Emp_DE.list)[1]= "EC_vs_EU_up"
Emp_DE.list[[2]] = EC_vs_EU %>% filter(log2FoldChange <=-1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(Emp_DE.list)[2]= "EC_vs_EU_down"
write.table(EC_vs_EU, file="results/Emp_6hrCpG_vs_0hr_DE.tsv", sep = "\t",quote = FALSE,row.names = FALSE)
volcano.plot(EC_vs_EU,"Control 6hr CpG vs 0hr")



# 0hr
NU_vs_EU <- results(dds.control, contrast=c("Conditon","NCoR1_Uns","Emp_Uns")) %>% as.data.frame() %>% rownames_to_column("Gene")
NCoR1_DE.list[[1]] = NU_vs_EU %>% filter(log2FoldChange >=1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(NCoR1_DE.list)[1]= "NU_vs_EU_up"
NCoR1_DE.list[[2]] = NU_vs_EU %>% filter(log2FoldChange <=-1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(NCoR1_DE.list)[2]= "NU_vs_EU_down"
write.table(NU_vs_EU, file="results/NCoR1_KD_vs_Emp_0hr_DE.tsv", sep = "\t",quote = FALSE,row.names = FALSE)
NCoR1.volcano.plot.list[[1]] = volcano.plot(NU_vs_EU,"NCoR1 KD vs Control (0hr)")



# 6hr
NC_vs_EC <- results(dds.control, contrast=c("Conditon","NCoR1_CpG_6h","Emp_CpG_6h")) %>% as.data.frame() %>% rownames_to_column("Gene")
NCoR1_DE.list[[3]] = NC_vs_EC %>% filter(log2FoldChange >=1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(NCoR1_DE.list)[3]= "NC_vs_EC_up"
NCoR1_DE.list[[4]] = NC_vs_EC %>% filter(log2FoldChange <=-1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(NCoR1_DE.list)[4]= "NC_vs_EC_down"
write.table(NC_vs_EC, file="results/NCoR1_KD_vs_Emp_6hr_DE.tsv", sep = "\t",quote = FALSE,row.names = FALSE)
NCoR1.volcano.plot.list[[2]] = volcano.plot(NC_vs_EC,"NCoR1 KD vs Control (6hr)")

#####################################################################################################################################
genes <- c("Il12b","Il10","Il6","Socs3","Nr4a1")
pdf("Figures/Gene_expression_bar_plot.pdf",height = 3.5,width = 3.5)
for (i in genes) {
  print(i) 
  df2 <- NCoR1_SMRT_norm_count %>% subset(NCoR1_SMRT_norm_count$Gene==i)
  print(df2)
  melted <- melt(df2)
  melted = merge(melted,col_data,by.x="variable",by.y="Sample") %>% dplyr::select(2,3,4,5)
  tgc <- summarySE(melted, measurevar="value", groupvars=c("grp","cnd"))
  tgc$condition <- c(rep("NCoR1",4),rep("SMRT",4))
  tgc$condition = factor(tgc$condition,levels = c("SMRT","NCoR1"))
  print(ggplot(tgc,aes(x=condition,y=value,fill=grp)) +
          geom_bar(stat="identity",position=position_dodge()) +
          facet_wrap(~cnd)+
          geom_errorbar( aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.9))+
          scale_fill_manual(values=c("grey","#DEB132","grey","#058983"))+
          scale_color_hue()+
          theme_bw()+
          theme(axis.text.x=element_text(size=10,angle =45,colour = "black",face= "bold",vjust = 1,hjust = 1),
                axis.text.y=element_text(size=10,colour = "black",face= "bold"),
                axis.title.y=element_text(size=10,colour = "black",face= "bold"),
                plot.title = element_text(size = 20, face = "italic",hjust = 0.5),
                legend.text = element_text(size=10,colour = "black"),
                legend.position = "none",
                strip.text = element_text(size = 15,face="bold")) + 
          labs(x = "",y="Normalized count",label =FALSE,title = "") +
          ggtitle(i))
  #print(p)
  
  
}
dev.off()


##########################################################################################
# Merge Replicate
merge_rep= function(col1){
  N <- ncol(col1)
  name <- colnames(col1)
  obj <- vector("list",ncol(col1)/2)
  k=1
  
  for(i in 1:N) {
    
    if(i%%2 ==1 && i <= N)    {
      
      #print(i)
      ID <-rowMeans(col1[,c(i,i+1)])
      obj[[k]] <- ID
      nam <- colnames(col1)[i]
      nam <- str_replace(nam,"_R[0-9]|.R[0-9].srt.bam|.R[0-9].bam","")
      names(obj)[k] <- nam
      names(obj[[k]]) <- rownames(col1)
      #print(k)
      k=k+1
    }
  }
  mat_merged <- as.data.frame(t(do.call(rbind, obj)))
  colnames(mat_merged) = names(obj)
  return(mat_merged)
}


