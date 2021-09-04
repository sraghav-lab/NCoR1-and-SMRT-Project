# Overlap NCoR1 and SMRT bound gene with upregulated and downregulated gene upon 6hr CpG stimulation 

ncor_smrt.diff= read.csv("../NCoR1_SMRT_analysis/tag_dir/peak_file/diff_peaks/NCoR1_SMRT_0h_6h_diff_peaks.bed",header = T,sep="\t")
ncor_smrt.diff= read.csv("../NCoR1_SMRT_analysis/tag_dir/peak_file/diff_peaks/NCoR1_SMRT_Uns_CpG_unique_peaks.bed",header = F,sep="\t")
View(ncor_smrt.diff)

ncor_smrt.diff1 = ncor_smrt.diff %>% filter(NC_vs_NU_up >=2 & SC_vs_SU_up >=2) 
ncor_smrt.diff2 = ncor_smrt.diff %>% filter(NC_vs_NU_down ==0 & NC_vs_SC_down ==0 & 
                                                    NU_vs_SU_down ==0 & SC_vs_SU_down == 0 &
                                                    NC_vs_NU_up  ==0 & NC_vs_SC_up >=2 &
                                                    NU_vs_SU_up >=2 & SC_vs_SU_up ==0) 


# Upregulated promoter proximal binidng 
up_pp = diff_NCoR1.SMRT.annotation.df %>% 
        filter(SYMBOL %in% Emp_DE.list[[1]]) %>% 
        filter(distanceToTSS <= 1000 & distanceToTSS >=-1000) %>% 
        dplyr::select(17,19) %>% unique() %>%
        dplyr::select(Clusters) %>% group_by(Clusters) %>% 
        dplyr::summarise(up_pp=n())

# Downregulated promoter proximal binding
down_pp = diff_NCoR1.SMRT.annotation.df %>% 
        filter(SYMBOL %in% Emp_DE.list[[2]]) %>% 
        filter(distanceToTSS <= 1000 & distanceToTSS >=-1000) %>% 
        dplyr::select(17,19) %>% unique() %>%
        dplyr::select(Clusters) %>% group_by(Clusters) %>% 
        dplyr::summarise(down_pp=n()) 

# Upregulated far distal binidng 
up_fd =  diff_NCoR1.SMRT.annotation.df %>% 
        filter(SYMBOL %in% Emp_DE.list[[1]]) %>% 
        filter(distanceToTSS > 1000 | distanceToTSS < -1000) %>% 
        dplyr::select(17,19) %>% unique() %>%
        dplyr::select(Clusters) %>% group_by(Clusters) %>% 
        dplyr::summarise(up_fd=n())

# Downregulated far distal binding
down_fd = diff_NCoR1.SMRT.annotation.df %>% 
        filter(SYMBOL %in% Emp_DE.list[[2]]) %>% 
        filter(distanceToTSS > 1000 | distanceToTSS < -1000) %>% 
        dplyr::select(17,19) %>% unique() %>%
        dplyr::select(Clusters) %>% group_by(Clusters) %>% 
        dplyr::summarise(down_fd=n())

ncor_smrt_binding = Reduce(merge, list(up_pp,down_pp,up_fd,down_fd))  %>% 
        melt()# %>% #filter(variable %in% c("up_pp","up_fd"))
#group_by(Clusters) %>% 
#mutate(percent=100*value/sum(value)) 
ncor_smrt_binding$Status = gsub("_.*","",ncor_smrt_binding$variable)
ncor_smrt_binding$distToTSS = gsub(".*_","",ncor_smrt_binding$variable)
ncor_smrt_binding = ddply(ncor_smrt_binding, .(Clusters,Status), transform, percent = value/sum(value) * 100)
ncor_smrt_binding = ddply(ncor_smrt_binding, .(Clusters), transform, pos = (cumsum(percent) - 0.7 * percent))
ncor_smrt_binding$label = paste0(sprintf("%.0f", ncor_smrt_binding$percent), "%")
ncor_smrt_binding$Clusters = factor(ncor_smrt_binding$Clusters,levels = c("NCoR1_Uns","NCoR1_CpG","NCoR1_SMRT_CpG","NCoR1_SMRT_Uns_CpG","SMRT_Uns","SMRT_CpG","SMRT_Uns_CpG"))
png("Figures/NCoR1_SMRT_binding_Control_DE_genes.png",height = 40,width = 8.5,units = "cm",res=400)
ggplot(ncor_smrt_binding, aes(x=Status, y=value,fill = distToTSS)) +
        geom_bar(stat = "identity") +
        geom_text(aes(y=value,label=paste0("\n",round(percent,2),"%","(",value,")")),position="stack",size=4)+
        facet_wrap(~Clusters,scales = "free",ncol=1,nrow = 7)+
        gg_theme+
        theme(legend.position = "top",legend.direction = "vertical")+
        scale_fill_manual(labels = c("Far Distal", "Promoter-proximal(± 1kb to TSS)"),values=c("#33B5FF","#FFD433"))+
        ylab("Number of Binding sites")
dev.off()
########################################################################################################################
NCoR1_SMRT_bound_DE.boxplot =list()
up = diff_NCoR1.SMRT.annotation.df %>% 
        filter(SYMBOL %in% Emp_DE.list[[1]] |SYMBOL %in% Emp_DE.list[[2]]) %>% 
        filter(SYMBOL %in% NCoR1_DE.list[[1]] | SYMBOL %in% NCoR1_DE.list[[2]]  |
                       SYMBOL %in% SMRT_DE.list[[1]] | SYMBOL %in% SMRT_DE.list[[2]]) %>% 
        filter(Clusters == "NCoR1_Uns") %>% 
        dplyr::select("SYMBOL") %>% pull(.,"SYMBOL") %>% unique() 
#names(NCoR1_Uns.down) = rep("down",length(NCoR1_Uns.down))

NCoR1_SMRT_bound_DE.boxplot[[1]]= NCoR1_SMRT_rld_rep_merged %>% 
        filter(rownames(.) %in% up) %>%
        mutate(group= case_when(rownames(.) %in% Emp_DE.list[[1]] ~ "up",
                                rownames(.) %in% Emp_DE.list[[2]] ~ "down")) %>%
        dplyr::select(3,6,8,10,11) %>%
        melt() %>%
        mutate(Group= case_when(variable %in% c("Ctrl_CpG_6h","SMRT_CpG_6h") ~ "SMRT KD",
                                variable %in% c("Emp_CpG_6h","NCoR1_CpG_6h") ~ "NCoR1 KD")) %>%
        mutate(variable=factor(variable, levels = c("Emp_CpG_6h","NCoR1_CpG_6h","Ctrl_CpG_6h","SMRT_CpG_6h"),
                               labels = c("Ctrl 6h","NCoR1 KD 6h","Ctrl 6h","SMRT KD 6h"))) %>%
        ggplot(.,aes(x=Group,y=value,color=variable))+
        geom_boxplot(width=0.4)+
        facet_wrap(~group)+
        scale_color_manual(name = "Group",values = c("grey", "#DEB132","#058983"))+
        gg_theme +
        ylab("rlog (DESeq2)") +xlab("")

up = diff_NCoR1.SMRT.annotation.df %>% 
        filter(SYMBOL %in% Emp_DE.list[[1]] |SYMBOL %in% Emp_DE.list[[2]]) %>%
        filter(SYMBOL %in% NCoR1_DE.list[[3]] | SYMBOL %in% NCoR1_DE.list[[4]]  |
                       SYMBOL %in% SMRT_DE.list[[5]] | SYMBOL %in% SMRT_DE.list[[6]]) %>% 
        filter(Clusters == "NCoR1_CpG") %>%
        dplyr::select("SYMBOL") %>% pull(.,"SYMBOL")%>% unique()

#names(NCoR1_Uns.down) = rep("down",length(NCoR1_Uns.down))

NCoR1_SMRT_bound_DE.boxplot[[2]]= NCoR1_SMRT_rld_rep_merged %>% 
        filter(rownames(.) %in% up) %>%
        mutate(group= case_when(rownames(.) %in% Emp_DE.list[[1]] ~ "up",
                                rownames(.) %in% Emp_DE.list[[2]] ~ "down")) %>%
        dplyr::select(3,6,8,10,11) %>%
        melt() %>%
        mutate(Group= case_when(variable %in% c("Ctrl_CpG_6h","SMRT_CpG_6h") ~ "SMRT KD",
                                variable %in% c("Emp_CpG_6h","NCoR1_CpG_6h") ~ "NCoR1 KD")) %>%
        mutate(variable=factor(variable, levels = c("Emp_CpG_6h","NCoR1_CpG_6h","Ctrl_CpG_6h","SMRT_CpG_6h"),
                               labels = c("Ctrl 6h","NCoR1 KD 6h","Ctrl 6h","SMRT KD 6h"))) %>%
        ggplot(.,aes(x=Group,y=value,color=variable))+
        geom_boxplot(width=0.4)+
        facet_wrap(~group)+
        scale_color_manual(name = "Group",values = c("grey", "#DEB132","#058983"))+
        gg_theme +
        ylab("rlog (DESeq2)") +xlab("")


up = diff_NCoR1.SMRT.annotation.df %>% 
        filter(SYMBOL %in% Emp_DE.list[[1]] |SYMBOL %in% Emp_DE.list[[2]]) %>%
        filter(SYMBOL %in% NCoR1_DE.list[[3]] | SYMBOL %in% NCoR1_DE.list[[4]]  |
                       SYMBOL %in% SMRT_DE.list[[5]] | SYMBOL %in% SMRT_DE.list[[6]]) %>% 
        filter(Clusters == "NCoR1_SMRT_CpG") %>%
        dplyr::select("SYMBOL") %>% pull(.,"SYMBOL")%>% unique()


NCoR1_SMRT_bound_DE.boxplot[[3]]= NCoR1_SMRT_rld_rep_merged %>% 
        filter(rownames(.) %in% up) %>%
        mutate(group= case_when(rownames(.) %in% Emp_DE.list[[1]] ~ "up",
                                rownames(.) %in% Emp_DE.list[[2]] ~ "down")) %>%
        dplyr::select(3,6,8,10,11) %>%
        melt() %>%
        mutate(Group= case_when(variable %in% c("Ctrl_CpG_6h","SMRT_CpG_6h") ~ "SMRT KD",
                                variable %in% c("Emp_CpG_6h","NCoR1_CpG_6h") ~ "NCoR1 KD")) %>%
        mutate(variable=factor(variable, levels = c("Emp_CpG_6h","NCoR1_CpG_6h","Ctrl_CpG_6h","SMRT_CpG_6h"),
                               labels = c("Ctrl 6h","NCoR1 KD 6h","Ctrl 6h","SMRT KD 6h"))) %>%
        ggplot(.,aes(x=Group,y=value,color=variable))+
        geom_boxplot(width=0.4)+
        facet_wrap(~group)+
        scale_color_manual(name = "Group",values = c("grey", "#DEB132","#058983"))+
        gg_theme +
        ylab("rlog (DESeq2)") +xlab("")

up = diff_NCoR1.SMRT.annotation.df %>% 
        filter(SYMBOL %in% Emp_DE.list[[1]] |SYMBOL %in% Emp_DE.list[[2]]) %>%
        filter(SYMBOL %in% NCoR1_DE.list[[3]] | SYMBOL %in% NCoR1_DE.list[[4]]  |
                       SYMBOL %in% SMRT_DE.list[[5]] | SYMBOL %in% SMRT_DE.list[[6]]) %>% 
        filter(Clusters == "NCoR1_SMRT_Uns_CpG") %>%
        dplyr::select("SYMBOL") %>% pull(.,"SYMBOL")%>% unique()


NCoR1_SMRT_bound_DE.boxplot[[4]]= NCoR1_SMRT_rld_rep_merged %>% 
        filter(rownames(.) %in% up) %>%
        mutate(group= case_when(rownames(.) %in% Emp_DE.list[[1]] ~ "up",
                                rownames(.) %in% Emp_DE.list[[2]] ~ "down")) %>%
        dplyr::select(3,6,8,10,11) %>%
        melt() %>%
        mutate(Group= case_when(variable %in% c("Ctrl_CpG_6h","SMRT_CpG_6h") ~ "SMRT KD",
                                variable %in% c("Emp_CpG_6h","NCoR1_CpG_6h") ~ "NCoR1 KD")) %>%
        mutate(variable=factor(variable, levels = c("Emp_CpG_6h","NCoR1_CpG_6h","Ctrl_CpG_6h","SMRT_CpG_6h"),
                               labels = c("Ctrl 6h","NCoR1 KD 6h","Ctrl 6h","SMRT KD 6h"))) %>%
        ggplot(.,aes(x=Group,y=value,color=variable))+
        geom_boxplot(width=0.4)+
        facet_wrap(~group)+
        scale_color_manual(name = "Group",values = c("grey", "#DEB132","#058983"))+
        gg_theme +
        ylab("rlog (DESeq2)") +xlab("")


up = diff_NCoR1.SMRT.annotation.df %>% 
        filter(SYMBOL %in% Emp_DE.list[[1]] |SYMBOL %in% Emp_DE.list[[2]]) %>% 
        filter(SYMBOL %in% NCoR1_DE.list[[1]] | SYMBOL %in% NCoR1_DE.list[[2]]  |
                       SYMBOL %in% SMRT_DE.list[[1]] | SYMBOL %in% SMRT_DE.list[[2]]) %>% 
        filter(Clusters == "SMRT_Uns") %>%
        dplyr::select("SYMBOL") %>% pull(.,"SYMBOL") %>% unique()
#names(NCoR1_Uns.down) = rep("down",length(NCoR1_Uns.down))

NCoR1_SMRT_bound_DE.boxplot[[5]]= NCoR1_SMRT_rld_rep_merged %>% 
        filter(rownames(.) %in% up) %>%
        mutate(group= case_when(rownames(.) %in% Emp_DE.list[[1]] ~ "up",
                                rownames(.) %in% Emp_DE.list[[2]] ~ "down")) %>%
        dplyr::select(3,6,8,10,11) %>%
        melt() %>%
        mutate(Group= case_when(variable %in% c("Ctrl_CpG_6h","SMRT_CpG_6h") ~ "SMRT KD",
                                variable %in% c("Emp_CpG_6h","NCoR1_CpG_6h") ~ "NCoR1 KD")) %>%
        mutate(variable=factor(variable, levels = c("Emp_CpG_6h","NCoR1_CpG_6h","Ctrl_CpG_6h","SMRT_CpG_6h"),
                               labels = c("Ctrl 6h","NCoR1 KD 6h","Ctrl 6h","SMRT KD 6h"))) %>%
        ggplot(.,aes(x=Group,y=value,color=variable))+
        geom_boxplot(width=0.4)+
        facet_wrap(~group)+
        scale_color_manual(name = "Group",values = c("grey", "#DEB132","#058983"))+
        gg_theme +
        ylab("rlog (DESeq2)") +xlab("")




up = diff_NCoR1.SMRT.annotation.df %>% 
        filter(SYMBOL %in% Emp_DE.list[[1]] |SYMBOL %in% Emp_DE.list[[2]]) %>%
        filter(SYMBOL %in% NCoR1_DE.list[[3]] | SYMBOL %in% NCoR1_DE.list[[4]]  |
                       SYMBOL %in% SMRT_DE.list[[5]] | SYMBOL %in% SMRT_DE.list[[6]]) %>% 
        filter(Clusters == "SMRT_CpG") %>%
        dplyr::select("SYMBOL") %>% pull(.,"SYMBOL")%>% unique()

#names(NCoR1_Uns.down) = rep("down",length(NCoR1_Uns.down))

NCoR1_SMRT_bound_DE.boxplot[[6]]= NCoR1_SMRT_rld_rep_merged %>% 
        filter(rownames(.) %in% up) %>%
        mutate(group= case_when(rownames(.) %in% Emp_DE.list[[1]] ~ "up",
                                rownames(.) %in% Emp_DE.list[[2]] ~ "down")) %>%
        dplyr::select(3,6,8,10,11) %>%
        melt() %>%
        mutate(Group= case_when(variable %in% c("Ctrl_CpG_6h","SMRT_CpG_6h") ~ "SMRT KD",
                                variable %in% c("Emp_CpG_6h","NCoR1_CpG_6h") ~ "NCoR1 KD")) %>%
        mutate(variable=factor(variable, levels = c("Emp_CpG_6h","NCoR1_CpG_6h","Ctrl_CpG_6h","SMRT_CpG_6h"),
                               labels = c("Ctrl 6h","NCoR1 KD 6h","Ctrl 6h","SMRT KD 6h"))) %>%
        ggplot(.,aes(x=Group,y=value,color=variable))+
        geom_boxplot(width=0.4)+
        facet_wrap(~group)+
        scale_color_manual(name = "Group",values = c("grey", "#DEB132","#058983"))+
        gg_theme +
        ylab("rlog (DESeq2)") +xlab("")


up = diff_NCoR1.SMRT.annotation.df %>% 
        filter(SYMBOL %in% Emp_DE.list[[1]] |SYMBOL %in% Emp_DE.list[[2]]) %>%
        filter(SYMBOL %in% NCoR1_DE.list[[3]] | SYMBOL %in% NCoR1_DE.list[[4]]  |
                       SYMBOL %in% SMRT_DE.list[[5]] | SYMBOL %in% SMRT_DE.list[[6]]) %>% 
        filter(Clusters == "SMRT_Uns_CpG") %>%
        dplyr::select("SYMBOL") %>% pull(.,"SYMBOL")%>% unique()

#names(NCoR1_Uns.down) = rep("down",length(NCoR1_Uns.down))

NCoR1_SMRT_bound_DE.boxplot[[7]]= NCoR1_SMRT_rld_rep_merged %>% 
        filter(rownames(.) %in% up) %>%
        mutate(group= case_when(rownames(.) %in% Emp_DE.list[[1]] ~ "up",
                                rownames(.) %in% Emp_DE.list[[2]] ~ "down")) %>%
        dplyr::select(3,6,8,10,11) %>%
        melt() %>%
        mutate(Group= case_when(variable %in% c("Ctrl_CpG_6h","SMRT_CpG_6h") ~ "SMRT KD",
                                variable %in% c("Emp_CpG_6h","NCoR1_CpG_6h") ~ "NCoR1 KD")) %>%
        mutate(variable=factor(variable, levels = c("Emp_CpG_6h","NCoR1_CpG_6h","Ctrl_CpG_6h","SMRT_CpG_6h"),
                               labels = c("Ctrl 6h","NCoR1 KD 6h","Ctrl 6h","SMRT KD 6h"))) %>%
        ggplot(.,aes(x=Group,y=value,color=variable))+
        geom_boxplot(width=0.4)+
        facet_wrap(~group)+
        scale_color_manual(name = "Group",values = c("grey", "#DEB132","#058983"))+
        gg_theme +
        ylab("rlog (DESeq2)") +xlab("")

pdf("Figures/NCoR1_SMRT_binding_Control_DE_genes_rlog.pdf",height = 12,width = 5)
do.call(grid.arrange,c(NCoR1_SMRT_bound_DE.boxplot,ncol=1))
dev.off()
png("Figures/NCoR1_SMRT_binding_Control_DE_genes_rlog.png",height = 30,width = 12,units = "cm",res=200)
do.call(grid.arrange,c(NCoR1_SMRT_bound_DE.boxplot,ncol=1))
dev.off()





##########################################################################################################
# Overlap of CpG and pIC specific genes with NCoR1 KD genes list
gs.RNASeq = 25000
NCoR1_SMRT_binding_DEG.overlap = newGOM(diff_NCoR1.SMRT.annotation.df.list,Emp_DE.list, gs.RNASeq)
NCoR1_SMRT_binding_DEG.overlap_pval = getMatrix(NCoR1_SMRT_binding_DEG.overlap,name="intersection") %>% 
        as.data.frame() %>% .[c(4,1,2,3,6,5,7),] %>%
        #-log10(.) %>%
        #sapply(.,signif,2) %>% 
        as.data.frame() %>%
        sapply(.,as.character)

p=drawHeatmap(NCoR1_SMRT_binding_DEG.overlap,ncolused=5, grid.col="Blues", note.col="black",adj.p = TRUE)
p=t(p$carpet) %>% as.data.frame() 
p=p[c(4,7,6,5,2,3,1),]
colnames(p) = paste0(colnames(p),"(",lengths(NCoR1_SMRT_binding_DEG.overlap@gsetB),")")
pdf("Figures/Emp_Up_down_gene_overlap_NCoR1_SMRT_binding.pdf",height = 5,width = 3.5)
draw(Heatmap(p,name = "Odds Ratio",
             col= colorRamp2(c(0,5,10),c("grey","yellow","red")),
             heatmap_legend_param=list(at=c(0,5,10),color_bar="continuous", 
                                       legend_direction="horizontal", legend_width=unit(5,"cm"),
                                       title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             rect_gp = gpar(col = "black"),
             column_names_side = "top",
             column_names_rot = 30,
             column_names_gp = gpar(fontsize = 15),
             row_names_gp = gpar(fontsize = 15),
             cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf(NCoR1_SMRT_binding_DEG.overlap_pval[i, j]), x, y, gp = gpar(fontsize = 10))}),
     heatmap_legend_side="bottom")
dev.off()



gs.RNASeq = 25000
# Overlap of CpG and pIC specific genes with NCoR1 KD genes list
NCoR1_SMRT_binding_DEG_KD.overlap = newGOM(diff_NCoR1.SMRT.annotation.df.list,c(NCoR1_DE.list,SMRT_DE.list[c(1,2,5,6)]), 
                                           gs.RNASeq)
NCoR1_SMRT_binding_DEG.KD.overlap_pval = getMatrix(NCoR1_SMRT_binding_DEG_KD.overlap,name="intersection") %>% 
        as.data.frame() %>% .[c(4,1,2,3,6,5,7),] %>% 
        #-log10(.) %>% 
        #sapply(.,signif,2) %>%  
        as.data.frame() %>% 
        sapply(.,as.character) 

p=drawHeatmap(NCoR1_SMRT_binding_DEG_KD.overlap,ncolused=5, grid.col="Blues", note.col="black",adj.p = TRUE)
p=t(p$carpet) %>% as.data.frame() 
p=p[c(4,7,6,5,2,3,1),]
colnames(p)=  paste0(colnames(p),"(",lengths(NCoR1_SMRT_binding_DEG_KD.overlap@gsetB),")")
rownames(p) = paste0(rownames(p),"(",lengths(NCoR1_SMRT_binding_DEG_KD.overlap@gsetA)[c(4,1,2,3,6,5,7)],")")
pdf("Figures/KD_Up_down_gene_overlap_NCoR1_SMRT_binding.pdf",height = 5,width = 6)
draw(Heatmap(p,name = "Odds Ratio",
             col= colorRamp2(c(0,2.5,5),c("grey","yellow","red")),
             heatmap_legend_param=list(at=c(0,2.5,5),color_bar="continuous", 
                                       legend_direction="horizontal", legend_width=unit(5,"cm"),
                                       title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             rect_gp = gpar(col = "black"),
             column_names_side = "top",
             column_names_rot = 30,
             column_names_gp = gpar(fontsize = 15),
             row_names_gp = gpar(fontsize = 15),
             cell_fun = function(j, i, x, y, width, height, fill) {
             grid.text(sprintf(NCoR1_SMRT_binding_DEG.KD.overlap_pval[i, j]), x, y, gp = gpar(fontsize = 10))}),
     heatmap_legend_side="bottom")
dev.off()


getMatrix(NCoR1_SMRT_binding_DEG_KD.overlap,name="intersection") %>% 
        as.data.frame() %>% rownames_to_column("Binding") %>% melt() 
        
NCoR1_specific =  diff_NCoR1.SMRT.annotation.df.list[c(1,4)] %>% unlist() %>% unname() %>% unique()
NCoR1_SMRT_common =  diff_NCoR1.SMRT.annotation.df.list[c(2,3)] %>% unlist() %>% unname() %>% unique()
SMRT_specific =  diff_NCoR1.SMRT.annotation.df.list[c(5,6,7)] %>% unlist() %>% unname() %>% unique()

NCoR1_SMRT_binding_DEG_KD.overlap.list = getNestedList(NCoR1_SMRT_binding_DEG_KD.overlap, name="intersection")
NCoR1_SMRT_binding_DEG.overlap.list = getNestedList(NCoR1_SMRT_binding_DEG.overlap, name="intersection")
extract_gene = function(x){
        a = x[c(1,4)] %>% unlist() %>% unname() %>% unique() %>% length()
        b = x[c(2,3)] %>% unlist() %>% unname() %>% unique() %>% length()
        c = x[c(5,6,7)] %>% unlist() %>% unname() %>% unique() %>% length()
        return(c(a,b,c))
}


percent.plt = do.call(rbind,lapply(NCoR1_SMRT_binding_DEG_KD.overlap.list, extract_gene)) %>% 
                        set_colnames(c("NCoR1 ", "NCoR1 & SMRT", "SMRT")) %>% as.data.frame() %>% 
                        rownames_to_column("Genes") %>% 
                        melt(.) %>% 
                        arrange(Genes) %>% 
                        mutate(Comparison =(gsub("_down|_up","",Genes)),
                               Status = (gsub(".*_vs_.*_","",Genes))) %>% 
                        dplyr::mutate(Comparison = factor(Comparison,levels = c("NU_vs_EU","SU_vs_EU","NC_vs_EC","SC_vs_EC.6"),
                                                                 labels = c("NCoR1 KD 0hr","SMRT KD 0hr","NCoR1 KD 6hr","SMRT KD 6hr"))) %>%
                        dplyr::select(c(2,4,5,3)) %>% 
                        dplyr::group_by(variable,Comparison) %>% 
                        dplyr::mutate(Percent = value/sum(value)*100) %>% as.data.frame() %>%
                        ggplot(.,aes(x=Status,y=Percent,fill=variable))+
                                geom_bar(stat="identity",position = "dodge")+
                                facet_wrap(~Comparison,scales = "free_x")+
                        gg_theme+
                        ylab("Percent of Genes") +
                        scale_fill_manual(name = "Bound",values = c("#DEB132","grey","#058983"))

pdf("Figures/NCoR1_SMRT_bound_genes_percent_DE.pdf",height = 4,width =4.5)
percent.plt
dev.off()
png("Figures/NCoR1_SMRT__bound_genes_percent_DE.png",height = 12,width = 13,units = "cm",res=400)
percent.plt
dev.off()
