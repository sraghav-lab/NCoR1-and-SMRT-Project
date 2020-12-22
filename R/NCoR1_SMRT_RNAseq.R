#count data
count <- read.csv("data/NCoR1_SMRT_CPG_all_condition_count.txt",sep = "\t",header = T)
head(count1)
rownames(count) <- count$Geneid
count1 <- count[,c(2:ncol(count))]
count1 <- count[,c(2:13,16,17,14,15,20,21,18,19)]
#######################################################################################################################
# 0hr and 6hr sample
count_dat <- count1[,c(1,2,13,14,5,6,15,16,7,8,11,12,17,18,19,20)]
head(count_dat)
group <- c(rep("Ctrl",8),rep("SMRT_KD",4),rep("NCoR1_KD",4))
#ind <- c(1,1,1,2,2,2,3,3,4,4,1,1,1,2,2,2,1,1,2,2)
cond <- c(rep("0hr",4),rep("6hr_CpG",4),rep(c("0hr","0hr","6hr_CpG","6hr_CpG"),2))
batch <- c(1,1,2,2,1,1,2,2,1,1,1,1,2,2,2,2)
grp_cnd <- factor(paste0(group,"_",cond))
col_data <- data.frame(Sample=colnames(count_dat),
                       grp=factor(group),
                       #ind=factor(ind),
                       cnd=factor(cond),
                       group=factor(grp_cnd),
                       batch=factor(batch))

dds <- DESeqDataSetFromMatrix(countData=count_dat, colData=col_data, design=~group+batch)
dds <- DESeq(dds)


vst <- vst(dds, blind=FALSE)
plt <- plotPCA(vst, intgroup=c("grp","cnd","batch"),returnData=TRUE)
percentVar <- round(100 * attr(plt, "percentVar"))
colnames(plt)[c(4,5,6)] = c("Group","Stimulation","Batch")
########################################################################################################
pdf("Figures/NCoR1_SMRT_sample_PCA.pdf",width = 8,height = 5)
ggplot(plt, aes(PC1, PC2, shape=Batch, color=Group)) + geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_mark_ellipse(aes(fill = Stimulation))+
  gg_theme
dev.off()
######################################################################################################
#normalized count
norm_count <- as.data.frame(counts(dds, normalized=TRUE))[,] %>% rownames_to_column("Gene")
write.table(norm_count, file="results/NCoR1_SMRT_all_condition_norm_count.tsv", sep = "\t",quote = FALSE)

######################################################################################################
rld<- rlogTransformation(dds, blind=TRUE)
rld = as.data.frame(assay(rld))
write.table(rld,file="results/NCoR1_SMRT_all_condition_rld.tsv",sep="\t",quote = FALSE)

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
      nam <- str_replace(nam,"_R[0-9]","")
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
# Average of log of variance stabilized normalized count 
NCoR1_SMRT_rld_rep_merged <- merge_rep(rld)
write.table(NCoR1_SMRT_rld_rep_merged,file="results/NCoR1_SMRT_rld_rep_merged.tsv",sep="\t",quote = FALSE)

#######################################################################################################
write.table(as.data.frame(assay(vst)),"results/NCoR1_SMRT_all_condition_vst.tsv",sep = "\t",quote = FALSE)
#######################################################################################################

# Since there were variability between SMRT Control and NCoR1 control sample due to batch effect we did 
# differential expression analyses for NCoR1 sample with its matched control

#count data
count <- read.csv("data/NCoR1_SMRT_CPG_all_condition_count.txt",sep = "\t",header = T,row.names = 1)
count <- count[,c(15,16,13,14,19,20,17,18)]

#######################################################################################################################
coldata <- data.frame(row.names=colnames(count),Conditon = gsub("_R[0-9]","",colnames(count)))
dds <- DESeqDataSetFromMatrix(countData=count, colData=coldata, design=~Conditon)
dds <- DESeq(dds)
#######################################################################################################################

NCoR1.volcano.plot.list = list()
NCoR1_DE.list = list()

# 0hr
NU_vs_EU <- results(dds, contrast=c("Conditon","NCoR1_Uns","Emp_Uns")) %>% as.data.frame() %>% rownames_to_column("Gene")
NCoR1_DE.list[[1]] = NU_vs_EU %>% filter(log2FoldChange >=1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(NCoR1_DE.list)[1]= "NU_vs_EU_up"
NCoR1_DE.list[[2]] = NU_vs_EU %>% filter(log2FoldChange <=-1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(NCoR1_DE.list)[2]= "NU_vs_EU_down"
write.table(NU_vs_EU, file="results/NCoR1_KD_vs_Emp_0hr_DE.tsv", sep = "\t",quote = FALSE,row.names = FALSE)
NCoR1.volcano.plot.list[[1]] = volcano.plot(NU_vs_EU,"NCoR1 KD vs Control (0hr)")



# 6hr
NC_vs_EC <- results(dds, contrast=c("Conditon","NCoR1_CpG_6h","Emp_CpG_6h")) %>% as.data.frame() %>% rownames_to_column("Gene")
NCoR1_DE.list[[3]] = NC_vs_EC %>% filter(log2FoldChange >=1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(NCoR1_DE.list)[3]= "NC_vs_EC_up"
NCoR1_DE.list[[4]] = NC_vs_EC %>% filter(log2FoldChange <=-1 & padj <= 0.05) %>% dplyr::select("Gene") %>% pull(.,"Gene")
names(NCoR1_DE.list)[4]= "NC_vs_EC_down"
write.table(NC_vs_EC, file="results/NCoR1_KD_vs_Emp_6hr_DE.tsv", sep = "\t",quote = FALSE,row.names = FALSE)
NCoR1.volcano.plot.list[[2]] = volcano.plot(NC_vs_EC,"NCoR1 KD vs Control (6hr)")



