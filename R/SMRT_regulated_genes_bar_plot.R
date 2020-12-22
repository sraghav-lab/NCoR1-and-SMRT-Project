####################################################################################################################
# Overlpa of SMRT DE genes with SMRT ChIP-seq target genes
gs.RNASeq = 25000 
SMRT_bound_DE = newGOM(SMRT_DE.list[c(1,2,5,6)],
                                NCoR1.SMRT.annotation.df.list[c(3,4)],
                                genome.size=gs.RNASeq)


SMRT_DE.bound = getMatrix(SMRT_bound_DE,name="pval") %>% 
                            as.data.frame() %>% 
                           # rownames_to_column("Comp") %>% 
                            melt() %>% .[c(3,4),]
SMRT_DE.unbound = data.frame(variable= c("SMRT CpG 6hr","SMRT CpG 6hr"),
                                   value=c(lengths(SMRT_DE.list[5])-SMRT_bound_DE.genes[1,2],
                                           lengths(SMRT_DE.list[6])-SMRT_bound_DE.genes[2,2]),
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

########################################################################################################################################################
#No. of  NCoR1 and SMRT DE genes at 0hr and 6hr
dat = data.frame(Comparison=names(c(lengths(SMRT_DE.list[c(1,2,5,6)]),lengths(NCoR1_DE.list))),
                 Number = c(lengths(SMRT_DE.list[c(1,2,5,6)]),lengths(NCoR1_DE.list)))
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
  scale_y_continuous(expand = c(0,0),limits = c(-1000,1300))+
pdf("Figures/NCoR1_SMRT_DE_genes_number.pdf",height = 6,width = 5)
p2
dev.off()
##########################################################################################################################
# List of direct targets of SMRT
NCoR1.SMRT.annotation.df %>% 
  filter(Clusters=="SMRT Uns") %>% 
  filter(SYMBOL %in% SMRT_DE.list[[1]]) %>% 
  write.table(.,file = "results/SMRT_bound_DE_up_genes_0hr.tsv",sep = "\t",quote = FALSE,row.names = FALSE)

NCoR1.SMRT.annotation.df %>% 
  filter(Clusters=="SMRT Uns") %>% 
  filter(SYMBOL %in% SMRT_DE.list[[2]]) %>% 
  write.table(.,file = "results/SMRT_bound_DE_down_genes_0hr.tsv",sep = "\t",quote = FALSE,row.names = FALSE)

NCoR1.SMRT.annotation.df %>% 
  filter(Clusters=="SMRT CpG 6hr") %>% 
  filter(SYMBOL %in% SMRT_DE.list[[5]]) %>% 
  write.table(.,file = "results/SMRT_bound_DE_up_genes_6hr.tsv",sep = "\t",quote = FALSE,row.names = FALSE)

NCoR1.SMRT.annotation.df %>% 
  filter(Clusters=="SMRT CpG 6hr") %>% 
  filter(SYMBOL %in% SMRT_DE.list[[6]]) %>% 
  write.table(.,file = "results/SMRT_bound_DE_down_genes_6hr.tsv",sep = "\t",quote = FALSE,row.names = FALSE)
