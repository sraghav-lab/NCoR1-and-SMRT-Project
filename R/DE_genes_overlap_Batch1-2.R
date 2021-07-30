#Overlap of DE from batch 1 and batch 2

Rep2_Rep5.01 =newGOM(SMRT_DE.list[c(1,2,5,6)],SMRT_DE.list.5Rep[c(1,2,3,4)])
Rep2_Rep3.01 =newGOM(SMRT_DE.list[c(1,2,5,6)],SMRT_DE.list.3Rep[c(1,2,3,4)])
Rep2_Rep5.05 =newGOM(SMRT_DE.list[c(1,2,5,6)],SMRT_DE.list.5Rep[c(5,6,7,8)])
Rep2_Rep3.05 =newGOM(SMRT_DE.list[c(1,2,5,6)],SMRT_DE.list.3Rep[c(5,6,7,8)])

Rep2_Rep3_Rep5 = list()

Rep2_Rep3_Rep5[[1]] = getMatrix(Rep2_Rep5.01,name="intersection") %>% 
                    as.data.frame() %>%
                    rownames_to_column("two_Rep") %>% 
                    melt(variable.name = "five_Rep",) %>%
                    filter(two_Rep == five_Rep) %>%
                    mutate(two_Rep_total = lengths(SMRT_DE.list[c(1,2,5,6)]),
                           five_Rep_total = lengths(SMRT_DE.list.5Rep[c(1,2,3,4)]),
                           Difference_two = two_Rep_total - value,
                           Difference_five = five_Rep_total - value,
                           value.1 = value) %>%
                    dplyr::select(c(1,3,6,7,8)) %>%
                    set_colnames(c("Comparison","Common","Difference_two","Difference_five","value.1")) %>%
                    melt() %>%
                    mutate(variable = gsub("value.1", "Common",variable),
                           Regulated = gsub(".*_","",Comparison),
                           Stimulation = c(rep(c(rep("0hr",2),rep("6hr CpG" ,2)),4)),
                           Datasets = c(rep("2 Rep\n (0.05)",8),rep("5 Rep\n (0.01)",8)),
                           variable = factor(variable, levels = c("Difference_three","Common","Difference_two"))) %>%     
                    ggplot(.,aes(x=Datasets,y= value))+
                      geom_bar(stat= "identity", 
                               aes(fill=factor(variable,levels = c("Difference_two","Difference_five","Common")))) +
                      facet_grid(Regulated~Stimulation) +
                      geom_text(aes(label=value,y=value), color = "black",position=position_stack(vjust=0.5),size=5) +
                      ylab("No. of genes") +
                      gg_theme +
                      guides(fill = guide_legend(title = ""))


Rep2_Rep3_Rep5[[2]] = getMatrix(Rep2_Rep3.01,name="intersection") %>% 
                    as.data.frame() %>%
                    rownames_to_column("two_Rep") %>% 
                    melt(variable.name = "five_Rep",) %>%
                    filter(two_Rep == five_Rep) %>%
                    mutate(two_Rep_total = lengths(SMRT_DE.list[c(1,2,5,6)]),
                           five_Rep_total = lengths(SMRT_DE.list.3Rep[c(1,2,3,4)]),
                           Difference_two = two_Rep_total - value,
                           Difference_five = five_Rep_total - value,
                           value.1 = value) %>%
                    dplyr::select(c(1,3,6,7,8)) %>%
                    set_colnames(c("Comparison","Common","Difference_two","Difference_three","value.1")) %>%
                    melt() %>%
                    mutate(variable = gsub("value.1", "Common",variable),
                           Regulated = gsub(".*_","",Comparison),
                           Stimulation = c(rep(c(rep("0hr",2),rep("6hr CpG" ,2)),4)),
                           Datasets = c(rep("2 Rep\n (0.05)",8),rep("3 Rep\n (0.01)",8)), 
                           variable = factor(variable, levels = c("Difference_three","Common","Difference_two"))) %>%
                    ggplot(.,aes(x=Datasets,y= value))+
                           geom_bar(stat= "identity", 
                                    aes(fill=factor(variable,levels = c("Difference_two","Difference_three","Common")))) +
                           facet_grid(Regulated~Stimulation) +
                           geom_text(aes(label=value,y=value), color = "black",position=position_stack(vjust=0.5),size=5) +
                           ylab("No. of genes") +
                           gg_theme +
                           guides(fill = guide_legend(title = ""))

Rep2_Rep3_Rep5[[3]] = getMatrix(Rep2_Rep5.05,name="intersection") %>% 
                            as.data.frame() %>%
                            rownames_to_column("two_Rep") %>% 
                            melt(variable.name = "five_Rep",) %>%
                            filter(two_Rep == five_Rep) %>%
                            mutate(two_Rep_total = lengths(SMRT_DE.list[c(1,2,5,6)]),
                                   five_Rep_total = lengths(SMRT_DE.list.5Rep[c(5,6,7,8)]),
                                   Difference_two = two_Rep_total - value,
                                   Difference_five = five_Rep_total - value,
                                   value.1 = value) %>%
                            dplyr::select(c(1,3,6,7,8)) %>%
                            set_colnames(c("Comparison","Common","Difference_two","Difference_five","value.1")) %>%
                            melt() %>%
                            mutate(variable = gsub("value.1", "Common",variable),
                                   Regulated = gsub(".*_","",Comparison),
                                   Stimulation = c(rep(c(rep("0hr",2),rep("6hr CpG" ,2)),4)),
                                   Datasets = c(rep("2 Rep\n (0.05)",8),rep("5 Rep\n (0.05)",8)),
                                   variable = factor(variable, levels = c("Difference_three","Common","Difference_two"))) %>%     
                            ggplot(.,aes(x=Datasets,y= value))+
                            geom_bar(stat= "identity", 
                                     aes(fill=factor(variable,levels = c("Difference_two","Difference_five","Common")))) +
                            facet_grid(Regulated~Stimulation) +
                            geom_text(aes(label=value,y=value), color = "black",position=position_stack(vjust=0.5),size=5) +
                            ylab("No. of genes") +
                            gg_theme +
                            guides(fill = guide_legend(title = ""))


Rep2_Rep3_Rep5[[4]] = getMatrix(Rep2_Rep3.05,name="intersection") %>% 
                            as.data.frame() %>%
                            rownames_to_column("two_Rep") %>% 
                            melt(variable.name = "five_Rep",) %>%
                            filter(two_Rep == five_Rep) %>%
                            mutate(two_Rep_total = lengths(SMRT_DE.list[c(1,2,5,6)]),
                                   five_Rep_total = lengths(SMRT_DE.list.3Rep[c(5,6,7,8)]),
                                   Difference_two = two_Rep_total - value,
                                   Difference_five = five_Rep_total - value,
                                   value.1 = value) %>%
                            dplyr::select(c(1,3,6,7,8)) %>%
                            set_colnames(c("Comparison","Common","Difference_two","Difference_three","value.1")) %>%
                            melt() %>%
                            mutate(variable = gsub("value.1", "Common",variable),
                                   Regulated = gsub(".*_","",Comparison),
                                   Stimulation = c(rep(c(rep("0hr",2),rep("6hr CpG" ,2)),4)),
                                   Datasets = c(rep("2 Rep\n (0.05)",8),rep("3 Rep\n (0.05)",8)), 
                                   variable = factor(variable, levels = c("Difference_three","Common","Difference_two"))) %>%
                            ggplot(.,aes(x=Datasets,y= value))+
                            geom_bar(stat= "identity", 
                                     aes(fill=factor(variable,levels = c("Difference_two","Difference_three","Common")))) +
                            facet_grid(Regulated~Stimulation) +
                            geom_text(aes(label=value,y=value), color = "black",position=position_stack(vjust=0.5),size=5) +
                            ylab("No. of genes") +
                            gg_theme +
                            guides(fill = guide_legend(title = ""))

pdf("Figures/SMRT_batch1-2_DE_RNASeq_comparison.pdf", height = 12, width = 15)
do.call(gridExtra::grid.arrange,c(Rep2_Rep3_Rep5,ncol=2,nrow=2))

# The number as inside the bar plot is not reversed due to code problem 
SMRT_rep_comp
dev.off()

