gs.RNASeq = 25000
# Overlap of CpG and pIC specific genes with NCoR1 KD genes list
NCoR1_SMRT_binding_DEG.overlap = newGOM(diff_NCoR1.SMRT.annotation.df.list,Emp_DE.list, gs.RNASeq)
NCoR1_SMRT_binding_DEG.overlap_pval = getMatrix(NCoR1_SMRT_binding_DEG.overlap,name="intersection") %>% 
                                  as.data.frame() %>% .[c(4,1,2,3,6,5,7),] %>%
                                  #-log10(.) %>%
                                  sapply(.,signif,2) %>% 
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
                                            sapply(.,signif,2) %>%  
                                            as.data.frame() %>% 
                                            sapply(.,as.character) 

p=drawHeatmap(NCoR1_SMRT_binding_DEG_KD.overlap,ncolused=5, grid.col="Blues", note.col="black",adj.p = TRUE)
p=t(p$carpet) %>% as.data.frame() 
p=p[c(4,7,6,5,2,3,1),]
colnames(p)=  paste0(colnames(p),"(",lengths(NCoR1_SMRT_binding_DEG_KD.overlap@gsetB),")")
rownames(p) = paste0(rownames(p),"(",lengths(NCoR1_SMRT_binding_DEG_KD.overlap@gsetA),")")
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
