# K-means clustering of Total Diff expressed genes  (log2 FC >=1) in NCoR1 and SMRT KD (0hr and 6hr)
df <- read.csv("data/NCoR1_SMRT_0_6_DE_genes.txt",sep = "\t",header = 1,row.names = 1)
head(df)

km<- kmeans(df,6)
m.kmeans<- cbind(df, km$cluster)  
dim(m.kmeans)
o<- order(m.kmeans[,5],decreasing = TRUE) # order the last column
m.kmeans<- m.kmeans[o,] # order the matrix according to the order of the last column
m.kmeans <- m.kmeans[,c(1,3,2,4,5)]
m.kmeans <- m.kmeans %>%
  mutate(Clus = case_when(
    m.kmeans$Cluster == 5 ~1,
    m.kmeans$Cluster == 3 ~2,
    m.kmeans$Cluster == 4 ~3,
    m.kmeans$Cluster == 1 ~4,
    m.kmeans$Cluster == 2 ~5,
    m.kmeans$Cluster == 6 ~6,
  ))
write.table(m.kmeans[,c(1:4,6)],file = "results/NCoR1_SMRT_0_6_DE_genes_cluster_batch1.txt",sep = "\t",quote = FALSE)
############################################################################################################
NCoR1_SMRT_m.kmeans = read.csv("results/NCoR1_SMRT_0_6_DE_genes_cluster_batch1.txt.txt",sep = "\t",header = T,row.names = 1)
NCoR1_SMRT_m.kmeans.cluster.list = split(rownames(NCoR1_SMRT_m.kmeans),NCoR1_SMRT_m.kmeans$Clus)
split <- factor(NCoR1_SMRT_m.kmeans[,5])
p2 <- Heatmap(NCoR1_SMRT_m.kmeans[,c(1,2,3,4)],km = 1,
              name="log2(FC)",
              col= colorRamp2(c(-2,0,2),c("#058983","white","#DEB132")),
              heatmap_legend_param=list(at=c(-3,0,3),color_bar="continuous", 
                                        legend_direction="horizontal", legend_width=unit(5,"cm"),
                                        title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
              #Split heatmap rows by gene family
              split = split,          
              #rect_gp = gpar(col = "black"),# box colour
              #Row annotation configurations
              cluster_rows=TRUE,
              cluster_row_slices = FALSE,
              show_row_dend=TRUE,
              #row_title="Transcript", #overridden by 'split' it seems
              row_title_side="left",
              row_title_gp=gpar(fontsize=8),
              show_row_names=FALSE,
              row_names_gp =gpar(fontsize=6, fontface="bold"),
              row_names_side="right",
              #row_title_rot=30,
              
              #Column annotation configuratiions
              cluster_columns=FALSE,
              show_column_dend=TRUE,
              #column_title="Samples",
              column_title_side="top",
              column_title_gp=gpar(fontsize=10, fontface="bold"),
              #column_title_rot=45,
              column_names_rot = 45,
              column_names_gp = gpar(fontsize = 10),
              show_column_names=TRUE,
              
              #Dendrogram configurations: columns
              clustering_distance_columns="euclidean",
              clustering_method_columns="single",
              column_dend_height=unit(10,"mm"),
              
              #Dendrogram configurations: rows
              clustering_distance_rows="euclidean",
              clustering_method_rows="single",
              row_dend_width=unit(10,"mm"))

pdf("Figures/NCoR1_SMRT_DE_genes_Kmeans_cluster.pdf",height = 7,width = 3)
draw(p2,heatmap_legend_side="top")
dev.off()

###############################################################################################
# Expression box plot for NCOR1 and SMRT cluster genes.
library(reshape2)
exp_file <- read.csv("results/NCoR1_SMRT_rld_rep_merged.tsv",sep="\t",header = T) %>% rownames_to_column("Gene")
head(exp_file)
pdf("Figures/NCoR1_SMRT_0_6hr_DE_genes_rld.pdf",height = 5,width = 7)
for (i in sort(unique(NCoR1_SMRT_m.kmeans$Clus))){
  print(i)
  
  files <- NCoR1_SMRT_m.kmeans %>% dplyr::filter(Clus == i)
  files$Gene <- rownames(NCoR1_SMRT_m.kmeans[which(NCoR1_SMRT_m.kmeans$Clus==i),])
  files <- exp_file %>% dplyr::filter(Gene %in% files$Gene)
  #name <- gsub("_rld.txt","",files[i])
  #rld <- read.csv(files,sep = "\t",header = T)
  rld <- files[,c(2,8,5,10,4,9,7,11)]
  rld <- t(scale(t(rld),center = TRUE))
  rld_melt <- melt(rld,variable.names = "sample")
  colnames(rld_melt) <- c("Var1","sample","value")
  rld_melt$Time <- gsub(".*_","",rld_melt$sample)
  rld_melt$Time <- gsub("Uns","0h",rld_melt$Time)
  rld_melt$Condition <- gsub("_.*","",rld_melt$sample)
  rld_melt$Condition <- gsub("Emp","Ctrl",rld_melt$Condition)
  rld_melt$Condition <- gsub("Crtl","Ctrl",rld_melt$Condition)
  rld_melt$Time <- factor(rld_melt$Time,levels = c("0h","6h"), labels = c("0hr","6hr CpG"))
  rld_melt$Condition <- factor(rld_melt$Condition,levels = c("Ctrl","SMRT","NCoR1"))
  
  #Compare means
  my_comparisons <- list( c("Ctrl", "SMRT"), c("Ctrl", "NCoR1"),c("SMRT","NCoR1"))
  #png(paste0(path,"NCoR1_SMRT_0_6_DE_genes_cluster",i,"_rld.png"),width = 16,height = 12,units = "cm",res=400)
  print(ggplot(rld_melt,aes(x=Condition,y=value))+
          stat_boxplot(geom="errorbar", width=.5)+
          #geom_violin()+
          facet_grid(.~Time) +
          ylim(-2.5,5)+
          ylab("Scaled rlog (DESeq2)")+
          scale_color_manual(values = c("grey","#058983","#DEB132"))+
          gg_theme_xrot+
          ggtitle(i)+
          geom_boxplot(position="dodge",aes(color=Condition))+
          stat_summary(fun=median, geom="line", aes(group=1), lwd=1,linetype=2)+
          stat_compare_means(comparisons = my_comparisons,tip.length=0.01,label = "p.format",
                             y_position = 7,method = "wilcox",paired = FALSE,))
  #dev.off()
}
dev.off()
##########################################################################################################
# Pathaway enrichment analysis of each cluster genes from IPA
path = "data/NCoR1_SMRT_0_6_DE_cluster_IPA/"
cluster.list = list.files("data/NCoR1_SMRT_0_6_DE_cluster_IPA/",pattern = "*_IPA.txt")
pdf("Figures/NCoR1_SMRT_0_6_DE_cluster_IPA.pdf",height = 3,width = 6)
cluster_pathway = read.csv(paste0(path,cluster.list[6]),header = T,sep = "\t",skip = 2) 
# cluster1 = c(1,3,4,6,8)
# cluster2 = c(5,10,11,23,27,28)
# cluster3 = c(3,5,8,10)
# cluster4 = c(5,6,8,15)
# cluster5 = c(2,4,8)
# cluster6 = c(1,2,3,6)
print(cluster_pathway[ c(1,2,3,6),] %>%
ggplot(aes(x=reorder(Ingenuity.Canonical.Pathways,X.log.p.value.), y=X.log.p.value.)) +
  geom_bar(stat='identity',fill = "#787878") +
  coord_flip() +
  gg_theme+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))+
  scale_y_continuous(expand = c(0,0))) 
dev.off()
##########################################################################################################