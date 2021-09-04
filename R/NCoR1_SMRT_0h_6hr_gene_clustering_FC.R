df = as.data.frame(unique(sort(c(unlist(SMRT_DE.list.5Rep[c(1,2,3,4)]),
                                 unlist(NCoR1_DE.list[c(1,2,3,4)],use.names = FALSE))))) %>% 
        set_colnames("Gene") %>% 
        merge(SU_vs_EU.5Rep,by="Gene",all.x ="TRUE") %>% 
        dplyr::select(1,3) %>% 
        set_colnames(c("Gene","SU_vs_EU")) %>% 
        merge(SC_vs_EC.6.5Rep,by="Gene",all.x ="TRUE") %>% 
        dplyr::select(1,2,4) %>% 
        set_colnames(c("Gene","SU_vs_EU","SC_vs_EC")) %>% 
        merge(NU_vs_EU,by="Gene",all.x ="TRUE") %>% 
        dplyr::select(1,2,3,5) %>% 
        set_colnames(c("Gene","SU_vs_EU","SC_vs_EC","NU_vs_NU")) %>% 
        merge(NC_vs_EC,by="Gene",all.x ="TRUE") %>% 
        dplyr::select(1,2,3,4,6) %>% 
        set_colnames(c("Gene","SU_vs_EU","SC_vs_EC","NU_vs_EU","NC_vs_EC")) %>% 
        column_to_rownames("Gene") %>% 
        filter(!(is.na(NC_vs_EC))) %>% 
        filter(!(is.na(SU_vs_EU)))
  
# K-means clustering of Total Diff expressed genes  (log2 FC >=1) in NCoR1 and SMRT KD (0hr and 6hr)
###############################################################################################
NCoR1_SMRT_kmeans.batch1 = read.csv("results/NCoR1_SMRT_0_6_DE_genes_cluster_batch1.txt",sep = "\t",header = T)

NCoR1_SMRT_kmeans.Rep5 = merge(df,NCoR1_SMRT_kmeans.batch1,by.x = 0,by.y="Gene") %>% 
                                dplyr::select(c(1,6:10)) 
NCoR1_SMRT_kmeans.Rep5.list = split(NCoR1_SMRT_kmeans.Rep5$Row.names,NCoR1_SMRT_kmeans.Rep5$Clus)

write.table(NCoR1_SMRT_kmeans.Rep5,file = "results/NCoR1_SMRT_0_6_DE_genes_cluster_batch2.txt",sep = "\t",quote = FALSE)


split.Rep5 = NCoR1_SMRT_kmeans.Rep5$Clus
p2 <- Heatmap(NCoR1_SMRT_kmeans.Rep5[,c(-5)],km = 1,
              name="log2(FC)",
              col= colorRamp2(c(-3,0,3),c("#058983","white","#DEB132")),
              heatmap_legend_param=list(at=c(-3,0,3),color_bar="continuous", 
                                        legend_direction="horizontal", legend_width=unit(5,"cm"),
                                        title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
              #Split heatmap rows by gene 
              split = split.Rep5,    
              border = TRUE,
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

# Expression box plot for NCOR1 and SMRT cluster genes.
exp_file <- read.csv("results/NCoR1_SMRT_vst_rep_merged.tsv",sep="\t",header = T) #%>% rownames_to_column("Gene")
head(exp_file)
pdf("Figures/NCoR1_SMRT_0_6hr_DE_genes_rld.pdf",height = 4.7,width = 6)
for (i in sort(unique(NCoR1_SMRT_kmeans.Rep5$Clus))){
  print(i)
  
  files <- NCoR1_SMRT_kmeans.Rep5 %>% dplyr::filter(Clus == i)
  files$Gene <- rownames(NCoR1_SMRT_kmeans.Rep5[which(NCoR1_SMRT_kmeans.Rep5$Clus==i),])
  files <- exp_file %>% dplyr::filter(Gene %in% files$Gene)
  #name <- gsub("_rld.txt","",files[i])
  #rld <- read.csv(files,sep = "\t",header = T)
  rld <- files
  rld <- t(scale(t(rld[,-1]),center = TRUE))
  rld_melt <- melt(rld,variable.names = "sample")
  colnames(rld_melt) <- c("Var1","sample","value")
  rld_melt$Time <- gsub("Ctrl.|Emp_|NCoR1_|SMRT.","",rld_melt$sample)
  rld_melt$Time <- gsub("Uns","0h",rld_melt$Time)
  rld_melt$Time <- gsub("_",".",rld_melt$Time)
  rld_melt$Condition <- gsub("_.*|\\..*","",rld_melt$sample)
  rld_melt$Condition <- gsub("Emp","Ctrl",rld_melt$Condition)
  #rld_melt$Condition <- gsub("Crtl","Ctrl",rld_melt$Condition)
  rld_melt$Time <- factor(rld_melt$Time,levels = c("0h","CpG.6h"), labels = c("0hr","6hr CpG"))
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
NCoR1_SMRT_m.kmeans.cluster.list.geneId = lapply(NCoR1_SMRT_m.kmeans.cluster.list,function(i) bitr(i, fromType = "SYMBOL",
                                                        toType = "ENTREZID",
                                                        OrgDb = org.Mm.eg.db)$ENTREZID)
compKEGG <- compareCluster(geneClusters   = NCoR1_SMRT_m.kmeans.cluster.list.geneId,
                           fun           = "enrichGO",
                           pvalueCutoff=0.05,
                           OrgDb=org.Mm.eg.db)

compKEGG = as.data.frame(compKEGG)
entrez = strsplit(compKEGG$geneID,"/")
gene = lapply(entrez,function(i) bitr(i, fromType = "ENTREZID",
                                      toType = "SYMBOL",
                                      OrgDb = org.Mm.eg.db)$SYMBOL)
pathway_genes = unique(unlist(gene))
gene = lapply(gene, function(i) paste0(i, sep="/", collapse=""))
gene = do.call(rbind,gene)
compKEGG$gene_symbol = gene

write.table(compKEGG,file = "results/NCoR1_SMRT_diff_peaks_KEGG.tsv",sep = "\t",quote = FALSE,row.names = FALSE)

kegg_pathway = c("mmu04062","mmu04659","mmu04151","mmu04630","mmu04010","mmu04060","mmu04668","mmu04660",
                 "mmu05235","mmu04658","mmu04620","mmu04668","mmu05323","mmu04064",
                 "mmu04657","mmu04150","mmu04612","mmu05168")

#dotplot(compKEGG, showCategory = 5, title = "KEGG Pathway Enrichment Analysis") 
compKEGG = compKEGG %>% separate(GeneRatio, c("A","B"), "/")
compKEGG$A = as.numeric(compKEGG$A)
compKEGG$B = as.numeric(compKEGG$B)
compKEGG = transform(compKEGG, GeneRatio=A / B)
compKEGG$Cluster = factor(compKEGG$Cluster,levels = c("SMRT_Uns_CpG","SMRT_CpG","SMRT_Uns","NCoR1_SMRT_CpG","NCoR1_SMRT_Uns_CpG","NCoR1_CpG","NCoR1_Uns"))
pdf("Figures/NCoR1_SMRT_diff_peaks_KEGG.pdf",height = 6,width = 12)
compKEGG %>% dplyr::filter(ID %in% kegg_pathway) %>%  #filter(Cluster  %in% c("NCoR1_CpG","SMRT_CpG")) %>% # View()
  ggplot(aes(x=Cluster,y=Description,colour=-log(p.adjust),size=GeneRatio))+
  geom_point()+
  gg_theme+scale_y_discrete(position = "right")+
  # coord_flip()
  theme(axis.text.x = element_text(colour = "black",angle = 45, vjust=0.8,hjust = 0.05,size=15),
        axis.text.y = element_text(size=15)) +
  scale_color_gradient(low = "blue", high = "red",) +
  scale_x_discrete(position = "top")+
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()





