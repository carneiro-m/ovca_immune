## STARTRAc analysis 

#load libraries
library(dplyr)
library(Seurat)
library(tidyr)
library(tibble)
library(Startrac)
library(ggplot2)

## Tcells obj -----

tcells <- readRDS("data/T_cell_cluster.seu.obj_meta.rds")#update in Dec 2021
cd8=unique(Idents(tcells))[grep('CD8', unique(Idents(tcells)))] %>% as.character()

star=(tcells@meta.data) %>% filter(patient %in% c("p09293", "p09454", "p09808" ,"p10329")) %>%
      select( clotype, patient, cell_mk, tissue) %>% 
      rownames_to_column(var='Cell_Name') %>% 
      rename( majorCluster=cell_mk, loc=tissue, clone.id=clotype) %>% 
      drop_na(clone.id) %>%
      filter(!majorCluster %in% c('myl', 'uk1', 'uk2', 'ukB', 'HS')) %>%
      mutate(across(where(is.factor), as.character))

star_cd8=(subset(tcells@meta.data, cell_mk %in% cd8.prolif)) %>% filter(patient %in% c("p09293", "p09454", "p09808" ,"p10329")) %>%
  select( clotype, patient, cell_mk, tissue) %>% 
  rownames_to_column(var='Cell_Name') %>% 
  rename( majorCluster=cell_mk, loc=tissue, clone.id=clotype) %>% 
  drop_na(clone.id) %>%
  mutate(across(where(is.factor), as.character))
          
out=Startrac.run(star, proj = 'ovca')   
out_cd8=Startrac.run(star_cd8, proj = 'ovca')

pdf("/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/startrac_all.pdf", width=6)
Startrac::plot(out,index.type="cluster.all",byPatient=T)
dev.off()
pdf("/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/startrac_blood-tumor.pdf", height = 5)
Startrac::plot(out,index.type="pairwise.migr",byPatient=T)+theme(
  axis.title.x = element_blank()
)
dev.off()
Startrac::plot(out,index.type="pairwise.tran",byPatient=F)

pdf("/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/startrac_CD8-all.pdf", width=6)
Startrac::plot(out_cd8,index.type="cluster.all",byPatient=T)
dev.off()

Startrac::plot(out_cd8,index.type="pairwise.tran",byPatient=F)

pdf("/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/startrac_blood-tumor_CD8.pdf", height = 5)
x= out_cd8@pIndex.sig.migr 
x$majorCluster=factor(x$majorCluster, levels=levels(Idents(cd8)))   
   ggplot(x,  aes(majorCluster, value))+
   geom_boxplot(outlier.shape = NA)+
   geom_jitter(aes( shape=aid), position=position_jitter(0.2))+
   theme_bw()+
   scale_colour_manual(values=met.brewer('Egypt'))+
   theme(axis.text.x=element_text(color=cd8.col, angle = 60, hjust=1), axis.title.x = element_blank())
  dev.off()

tissue_star=Startrac::calTissueDist(star_cd8) %>% as.data.frame() %>%
            rename(cluster=Var1, tissue=Var2) %>%
            pivot_wider(names_from = tissue, values_from=Freq)
 
ggplot(tissue_star, aes(Var2, Freq))+
  geom_line()+
  facet_grid(. ~Var1)

 ggplot(tissue_star) +
   geom_segment( aes(x=cluster, xend=cluster, y=blood, yend=tumor), color="grey") +
   geom_point( aes(x=cluster, y=blood, color='#00838F'), size=3 ) +
   geom_point( aes(x=cluster, y=tumor, color= '#C0392B'), size=3 ) +
   #coord_flip() + 
   theme_classic()+
   theme(panel.background = element_rect(color='darkgrey'), axis.line = element_blank(), 
         axis.text.x = element_text(angle=90))+
   xlab("")+ylab("") +scale_color_discrete(name = "Tissue", labels = c('blood', 'tumor'))
 
 
 ## Tumor obj ----
 colnames(tumor@meta.data)
 table(tumor$patient, useNA = 'always')
 star.tumor=(tumor@meta.data) %>% filter(patient %in% c("p09293", "p09454", "p09808" ,"p10329")) %>%
   select( clotype, patient, cell_mk, tissue) %>% 
   rownames_to_column(var='Cell_Name') %>% 
   rename( majorCluster=cell_mk, loc=tissue, clone.id=clotype) %>% 
   drop_na(clone.id) %>%
   filter(!majorCluster %in% c('myl', 'uk1', 'uk2', 'ukB', 'HS')) %>%
   mutate(across(where(is.factor), as.character))