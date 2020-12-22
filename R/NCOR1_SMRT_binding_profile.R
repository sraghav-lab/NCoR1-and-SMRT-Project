####################################################################################################################
# Annotate NCoR1 and SMRT peak bed files
files = list.files(path = "data/",pattern = "*peak.bed",full.names = TRUE)
names(files) = c("NCoR1 Uns","NCoR1 CpG 6hr","SMRT Uns","SMRT CpG 6hr")
peakAnnoList <- lapply(files, annotatePeak,
                       TxDb=txdb,
                       tssRegion=c(-1000, 1000), 
                       annoDb = "org.Mm.eg.db",
                       genomicAnnotationPriority = c("Promoter", "Exon", "Intron", "Intergenic"))
NCoR1.SMRT.annotation.df = lapply(peakAnnoList, function(i) as.data.frame(i)) %>%
                                       do.call(rbind, .) %>% 
                                       mutate(Clusters= gsub("\\..*","",rownames(.),""))
NCoR1.SMRT.annotation.df.list = split(NCoR1.SMRT.annotation.df$SYMBOL,NCoR1.SMRT.annotation.df$Clusters)
NCoR1.SMRT.annotation.df.list = lapply(NCoR1.SMRT.annotation.df.list, unique)
pdf("Figures/NCoR1_SMRT_peak_distribution.pdf",height = 3.5,width = 9)
plotAnnoBar(peakAnnoList)+gg_theme
dev.off()
####################################################################################################################
#Annotate NCoR1 and SMRT differential peaks
files =  list.files(path = "data/","*_unique_peaks.bed$",full.names = TRUE)
files = files[c(4,1,3,2,7,5,6)]
file_name =  list.files(path = "data","*_unique_peaks.bed$")
file_name = gsub("_unique_peaks.bed","",file_name)
file_name = file_name[c(4,1,3,2,7,5,6)]
names(files) = file_name

peakAnnoList <- lapply(files, annotatePeak, 
                       TxDb=txdb,
                       tssRegion=c(-1000, 1000),
                       annoDb = "org.Mm.eg.db",
                       genomicAnnotationPriority = c("Promoter", "Exon", "Intron", "Intergenic"))

diff_NCoR1.SMRT.annotation.df = lapply(peakAnnoList, function(i) as.data.frame(i)) %>%
                                    do.call(rbind, .) %>% 
                                    mutate(Clusters= gsub("\\..*","",rownames(.),""))
#####################################################################################################################
#Extract number of features 
annotation = function(x){
    as.data.frame(x) %>% 
    dplyr::select(annotation) %>%
    separate(annotation, c("A","B"), " \\(",extra = "drop",fill = "right") %>% 
    dplyr::select(A) %>% table()  %>% as.data.frame() 
}
#####################################################################################################################
# Number of NCoR1 and SMRT Peak Distributed relative to TSS 
peak_annotation = lapply(peakAnnoList, annotation)
peak_annotation = do.call(cbind,peak_annotation)[c(1,3,4,5),c(1,2,4,6,8,10,12,14)]
colnames(peak_annotation) = gsub("\\.Freq","",colnames(peak_annotation))
colnames(peak_annotation)[1] = "Features"
peak_annotation_melt = melt(peak_annotation)
peak_annotation_melt$variable = factor(peak_annotation_melt$variable,
                                       levels = c("SMRT_Uns_CpG","SMRT_CpG","SMRT_Uns","NCoR1_SMRT_CpG",
                                                  "NCoR1_SMRT_Uns_CpG","NCoR1_CpG","NCoR1_Uns"))

pdf("Figures/NCoR1_SMRT_diff.peaks_distribition.pdf",height = 7,width = 12)
print(peak_annotation %>% melt() %>%
  ggplot(aes(x=variable,y=value,fill=Features))+
  geom_bar(stat = "identity")+
  facet_wrap(~variable,scales = "free",ncol = 7)+
  #coord_flip() +
  gg_theme+
  theme(axis.text.x = element_blank())) 
dev.off()
#####################################################################################################################
# KEGG pathway enrichment analysis for NCoR1 and SMRT differential peaks
genes = split(diff_NCoR1.SMRT.annotation.df$geneId,diff_NCoR1.SMRT.annotation.df$Clusters)
genes = lapply(genes, unique)
compKEGG <- compareCluster(geneClusters   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff=0.05,
                           org="mmu")

compKEGG = as.data.frame(compKEGG)
entrez = strsplit(compKEGG$geneID,"/")
gene = lapply(entrez,function(i) bitr(i, fromType = "ENTREZID",
                                      toType = "SYMBOL",
                                      OrgDb = org.Mm.eg.db)$SYMBOL)
pathway_genes = unique(unlist(gene))
gene = lapply(gene, function(i) paste0(i, sep="/", collapse="") )
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
compKEGG %>% dplyr::filter(ID %in% kegg_pathway) %>%  #!(Cluster  %in% c("NCoR1_CpG","SMRT_CpG"))) %>% # View()
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
################################################################################################################
#Jak-Stat signalling pathway gene expression heatmap 
jak_stat <- c("Il10","Stat5b","Stat3","Socs3","Pim1","Il4ra","Pik3cb","Pik3cd","Sos1","Il15ra","Bcl2l1","Il21r","Lifr")
df <- read.csv("data/NCoR1_SMRT_all_condition_DE_genes_1.5fold_list_FC.txt",sep = "\t",header = 1,row.names = 1)
head(df)
jak_stat_df = df %>% dplyr::filter(rownames(df) %in% jak_stat)
jak_stat_map = Heatmap(jak_stat_df,cluster_columns = FALSE,name = "Log2(Fold Change)",
        col= colorRamp2(c(-2,0,2),c("#058983","white","#DEB132")),
        rect_gp = gpar(col= "black"),
        heatmap_legend_param=list(at=c(-3,0,3),color_bar="continuous", 
                                  legend_direction="horizontal", legend_width=unit(5,"cm"),
                                  title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
        column_names_rot = 45, column_names_side = "top",)
pdf("Figures/Jak_stat_signalling_genes_exp.pdf",height = 6,width = 3)
draw(jak_stat_map,heatmap_legend_side="bottom")
dev.off()
################################################################################################################
motifs = read.csv("data/Enriched_motifs_NCoR1_SMRT_diff_peaks.txt",sep = "\t",header = T)
motifs_melt = melt(motifs)
motifs_melt$Category =factor(motifs_melt$Category)
motifs_melt$variable =factor(motifs_melt$variable,levels = resetwd("/home/imgsb/Gyan/NCOR1/NCoR1_SMRT_analysis/tag_dir/peak_file/diff_peaks")v(c("PU.1","RUNX","NFkB","Jun.Fos.Fra1","Stat3","IRF4","IRF8")))
my_pal <- colorRampPalette(brewer.pal(11, "OrRd"))(25)
pdf("Figures/NCoR1_SMRT_Enriched_Motifs_.pdf",height = 7,width = 7)
ggplot(motifs_melt,aes(x=variable,y=Category))+
      geom_tile(aes(fill = value),colour = "black")+
      theme_bw()+
      theme(axis.text.x = element_text(colour = "black",angle = 45, vjust=0.8,hjust = 0.05,size=15),
            axis.text.y = element_text(colour = "black",size=15),
            legend.position="right") +
      scale_x_discrete(position = "bottom")+
      scale_y_discrete(position = "right")+
      coord_flip()+
      scale_fill_gradientn(name="Percent",colors = my_pal, na.value = 'white')

dev.off()
#################################################################################################################