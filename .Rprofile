library(Rsubread)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(tibble)
library(ggforce)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(dendextend)
library(RColorBrewer)
library(circlize)
library(reshape2)
library(ggpubr)
library(rtracklayer)
library(org.Mm.eg.db)
library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(tidyr)
library(gridExtra)
library(GeneOverlap)
library(tidytext)
library(Rmisc)
gg_theme <- 
  gg_theme <- theme_bw(15) + 
  theme(axis.text.x=element_text(size=15,color="black"),
        axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15,color="black"),
        axis.title.y=element_text(size=15),
        legend.justification=c(0,1),
        #legend.box.background = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=10))
gg_theme_xrot <- theme_bw(15) + 
  theme(axis.text.x=element_text(size=15,color="black",hjust = 1,vjust=1,angle = 45),
        axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15,color="black"),
        axis.title.y=element_text(size=15),
        legend.justification=c(0,1),
        #legend.box.background = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=10))

