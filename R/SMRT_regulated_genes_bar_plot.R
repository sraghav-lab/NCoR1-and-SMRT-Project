####################################################################################################################
# Overlpa of SMRT DE genes with SMRT ChIP-seq target genes
gs.RNASeq = 25000 
SMRT_bound_DE = newGOM(SMRT_DE.list.5Rep[c(1,2,3,4)],
                                NCoR1.SMRT.annotation.df.list[c(3,4)],
                                genome.size=gs.RNASeq)
SMRT_DE.bound.pval = getMatrix(SMRT_bound_DE,name="pval") %>% 
                        as.data.frame() 


SMRT_DE.bound.NestedList = getNestedList(SMRT_bound_DE,name="intersection") 

SMRT_DE.bound = getMatrix(SMRT_bound_DE,name="intersection") %>% 
                            as.data.frame() %>% 
                            #rownames_to_column("Comp") %>% 
                            melt() %>% .[c(3,4),]
SMRT_DE.unbound = data.frame(variable= c("SMRT CpG 6hr","SMRT CpG 6hr"),
                                   value=c(lengths(SMRT_DE.list.5Rep[3])-SMRT_DE.bound[1,2],
                                           lengths(SMRT_DE.list.5Rep[4])-SMRT_DE.bound[2,2]),
                                   row.names = NULL)
SMRT_DE.regulated.genes =  rbind(SMRT_DE.bound,SMRT_DE.unbound)
SMRT_DE.regulated.genes$Bound <- factor(c("SMRT_Bound","SMRT_Bound","SMRT_Unbound","SMRT_Unbound"),levels = c("SMRT_Unbound","SMRT_Bound"))
SMRT_DE.regulated.genes$Reg <- factor(c("UP","DOWN","UP","DOWN"),levels = c("UP","DOWN"))

pdf("Figures/SMRT_6hr_Bound.genes.barplot.pdf",height = 4,width = 5.5)
ggplot(data=SMRT_DE.regulated.genes, aes(x=Reg, y=value, fill=Bound)) +
  ylim(0,1000)+
  geom_bar(stat='identity', position='stack')+#facet_grid( ~V3)+
  scale_fill_manual(name = "",values = c("#058983","#DEB132"))+
  geom_text(data=SMRT_DE.regulated.genes,aes(label=abs(value),y=value), vjust=1, hjust=0.5,color = "black",position=position_stack(),size=5)+
  gg_theme_xrot+
  labs(x = "",y="No. of Genes",label =FALSE,title = "") +
  guides(fill=guide_legend(title="",label.theme = element_text(angle = 0,size = 15))) +
  scale_y_continuous(expand = c(0,0),breaks = c(200,400,600,800,1000))
dev.off()



# Write genes that are up and bound by SMRT
SMRT_bound_DE.list = list()
SMRT_bound_DE.list[[1]] = SMRT_DE.bound.NestedList$`SMRT CpG 6hr`$SC_vs_EC.6_up
SMRT_bound_DE.list[[2]] = SMRT_DE.bound.NestedList$`SMRT CpG 6hr`$SC_vs_EC.6_down
names(SMRT_bound_DE.list) = c("SMRT_Bound_up_genes","SMRT_Bound_down_genes")

SMRT_DE.bound.NestedList$`SMRT CpG 6hr`$SC_vs_EC.6_up %>% 
  as.data.frame() %>% 
  set_colnames("SMRT_Bound_up_genes") %>% 
  write.table(.,file = "results/SMRT_bound_DE_up_genes_6hr_CpG.list.tsv",sep = "\t",quote = FALSE,row.names = FALSE)

SMRT_DE.bound.NestedList$`SMRT CpG 6hr`$SC_vs_EC.6_down %>% 
  as.data.frame() %>% dim()
  set_colnames("SMRT_Bound_down_genes") %>% 
  write.table(.,file = "results/SMRT_bound_DE_down_genes_6hr_CpG.list.tsv",sep = "\t",quote = FALSE,row.names = FALSE)
########################################################################################################################################################
#No. of  NCoR1 and SMRT DE genes at 0hr and 6hr
dat = data.frame(Comparison=names(c(lengths(SMRT_DE.list.5Rep[c(1,2,3,4)]),lengths(NCoR1_DE.list))),
                 Number = c(lengths(SMRT_DE.list.5Rep[c(1,2,3,4)]),lengths(NCoR1_DE.list)))
dat$Reg = factor(c(rep(c("Up","Down"),4)))
dat$comp = gsub("_up|_down","",dat$Comparison)
dat$comp <- factor(dat$comp,levels = c("SU_vs_EU","NU_vs_EU","SC_vs_EC.6","NC_vs_EC"),
                            labels =c("SMRT_Uns","NCoR1_Uns","SMRT_CpG","NCoR1_CpG"))
dat= dat %>% mutate(Number = ifelse(Reg =="Down",dat$Number*-1,dat$Number))
p2 <- ggplot(dat , aes(x=comp,y=Number,fill=Reg))+ 
  geom_bar(width=0.8,stat="identity")+
  geom_text(data=dat,aes(label=abs(Number),y=Number), vjust=1, hjust=0.5,color = "black",position=position_stack(),size=5)+
  scale_fill_manual(name = "",values = c("#DEB132","#058983"))+
  theme_bw(base_size = 20) +
  gg_theme_xrot+
  labs(x = "",y="No. of Genes",label =FALSE,title = "") + 
  scale_y_continuous(expand = c(0,0),limits = c(-1000,1300))
pdf("Figures/NCoR1_SMRT_DE_genes_number.pdf",height = 6,width = 5)
p2
dev.off()
##########################################################################################################################
# List of direct targets of SMRT
NCoR1.SMRT.annotation.df %>% 
  dplyr::filter(Clusters=="SMRT Uns") %>% 
  dplyr::filter(SYMBOL %in% SMRT_DE.list.5Rep[[1]]) %>% #dim() 
  write.table(.,file = "results/SMRT_bound_DE_up_genes_0hr.tsv",sep = "\t",quote = FALSE,row.names = FALSE)

NCoR1.SMRT.annotation.df %>% 
  filter(Clusters=="SMRT Uns") %>% 
  filter(SYMBOL %in% SMRT_DE.list.5Rep[[2]]) %>% dim()
  write.table(.,file = "results/SMRT_bound_DE_down_genes_0hr.tsv",sep = "\t",quote = FALSE,row.names = FALSE)

NCoR1.SMRT.annotation.df %>% 
  filter(Clusters=="SMRT CpG 6hr") %>% 
  filter(SYMBOL %in% SMRT_DE.list.5Rep[[3]]) %>% dim()
  write.table(.,file = "results/SMRT_bound_DE_up_genes_6hr.tsv",sep = "\t",quote = FALSE,row.names = FALSE)

NCoR1.SMRT.annotation.df %>% 
  filter(Clusters=="SMRT CpG 6hr") %>% 
  filter(SYMBOL %in% SMRT_DE.list.5Rep[[4]]) %>% #dim()
  write.table(.,file = "results/SMRT_bound_DE_down_genes_6hr.tsv",sep = "\t",quote = FALSE,row.names = FALSE)

##########################################################################################################################
# List of direct targets of NCoR1
NCoR1.SMRT.annotation.df %>% 
  dplyr::filter(Clusters=="NCoR1 Uns") %>% 
  dplyr::filter(SYMBOL %in% NCoR1_DE.list[[1]]) %>% #dim() 
  write.table(.,file = "results/NCoR1_bound_DE_up_genes_0hr.tsv",sep = "\t",quote = FALSE,row.names = FALSE)

NCoR1.SMRT.annotation.df %>% 
  filter(Clusters=="NCoR1 Uns") %>% 
  filter(SYMBOL %in% NCoR1_DE.list[[2]]) %>% 
write.table(.,file = "results/NCoR1_bound_DE_down_genes_0hr.tsv",sep = "\t",quote = FALSE,row.names = FALSE)

NCoR1.SMRT.annotation.df %>% 
  filter(Clusters=="NCoR1 CpG 6hr") %>% 
  filter(SYMBOL %in% NCoR1_DE.list[[3]]) %>% 
write.table(.,file = "results/NCoR1_bound_DE_up_genes_6hr.tsv",sep = "\t",quote = FALSE,row.names = FALSE)

NCoR1.SMRT.annotation.df %>% 
  filter(Clusters=="NCoR1 CpG 6hr") %>% 
  filter(SYMBOL %in% NCoR1_DE.list[[4]]) %>% #dim()
  write.table(.,file = "results/NCoR1_bound_DE_down_genes_6hr.tsv",sep = "\t",quote = FALSE,row.names = FALSE)

##########################################################################################################################

gene = lapply( SMRT_DE.list.5Rep[c(3,4)],function(i) bitr(i, fromType = "SYMBOL",
                                        toType = "ENTREZID",
                                        OrgDb = org.Mm.eg.db)$ENTREZID)
genes = lapply(genes, unique)
compKEGG <- compareCluster(geneClusters   = gene,
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
  