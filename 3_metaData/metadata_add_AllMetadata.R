################################################################################
############################ T cells metadata ##################################
# organize Tcell metadata do add extra info: 
# 1. Add leiden cluster - from PAGA analysis on CD8 T cells in python - scanpy
# 2. TCR specificity - result from analysis in file "TCR_specificity.R"
# 3. Score cells with signatures from published articles: Gueguen, SciImmu2021, Szabo et al, 2019 - Nat Comm, and Savas, Nat Med 2018
# 4. Rename clusters using CD4 or CD8 + the name of top expressed gene

#load libraries -----
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)

#load data -----
tcells <- readRDS("~/analysis/data/T_cell_cluster.seu.obj_TCR.rds")
tcells$Cell.ID <- factor(tcells$Cell.ID, levels=c("CD4_1","CD4_2","CD4_3","CD4_4", "CD4_5","Treg", "CD8_1","CD8_2","CD8_3", "CD8_4","CD8_5","CD8_6","CD8_7","CD8_8", "HS","myl", "prolif", "ukB", "uk1",  "uk2"))
Idents(tcells) <- 'Cell.ID'
DefaultAssay(tcells)='SCT'
tcells$barcode=rownames(tcells@meta.data)


################################################################################
##### 1. Add leiden cluster
#load data
leiden <- read.csv("data/tcell_cd8s_bl_tu_leiden.csv") 
leiden$barcode=leiden$Unnamed..0
leiden=leiden %>% mutate(leiden=paste('k', leiden, sep='_'))
str(leiden)
n_k=length(unique(leiden$leiden))

tcells@meta.data = tcells@meta.data %>% left_join(leiden[c('leiden', 'barcode')], by='barcode') %>%
  mutate(new_id=coalesce(leiden,Cell.ID )) %>% column_to_rownames('barcode')
levels_leiden=paste(rep('k',n_k), seq(from=0, to=n_k-1), sep='_')
tcells$leiden=factor(tcells$leiden, levels=levels_leiden)
tcells$new_id=factor(tcells$new_id, levels = c("CD4_1","CD4_2","CD4_3","CD4_4", 
                                               "CD4_5","Treg",levels_leiden,"HS","myl", 
                                               "prolif", "ukB", "uk1",  "uk2" ))


################################################################################
##### 2. Add TCR specificity 
ov.ident=read.csv('TCR_specificity.csv')

tcells@meta.data$barcode=rownames(tcells@meta.data)
tcells@meta.data=tcells@meta.data %>% left_join(ov.ident %>% select(Pathology_groups, Pathology, barcode), by='barcode')
tcells$Pathology_groups=factor(tcells$Pathology_groups, levels = 
                                 c('tumor', 'virus', 'autoimmune',
                                   'others', 'unknown-TCR'))
tcells@meta.data =tcells@meta.data %>% column_to_rownames(var='barcode')


################################################################################
##### 3. Score cells with signatures from published articles

#subset cd8 for score 
cd8.id=levels(tcells$new_id)[7:30]
Idents(tcells)='new_id'
cd8=subset(tcells, idents=cd8.id)
til=subset(cd8, tissue='tumor')


#compare with Gueguen et al, 2021 signatures - Sci Imm  --------------------------------
gueguen <- read_csv("abd5778_Table_S2.csv")
unique(gueguen$cluster)
gueguen_cd8= gueguen %>% filter(cluster %in% c("CD8-SLC4A10","CD8-GZMH", "CD8-KLF2",  "CD8-XCL1",
                                               "CD8-GZMK" ,    "CD8-LAYN" , "CD8-FCGR3A")) %>%
                         filter(p_val_adj <0.05) %>% filter(avg_logFC>=0.4)
gueguen_list=list()
for (k in unique(gueguen_cd8$cluster)){
  gueguen_list[[k]]=(gueguen_cd8 %>% filter(cluster==k))$gene %>% list()
  
}
for (k in unique(gueguen_cd8$cluster)){
  cd8=AddModuleScore(cd8, features = gueguen_list[[k]], name=k)
}

### Compare with other signatures used in Gueguen et al, 2021
#fig 2 from Gueguen et al:
fig_2 <- read_csv("abd5778_Table_S3.csv")
fig_2.list=list()
for(f in colnames(fig_2)){
  fig_2.list[[f]]=fig_2[,f]%>% drop_na() %>% as.list()
}

colnames(fig_2)
for (f in colnames(fig_2)){
  cd8=AddModuleScore(cd8, features = fig_2.list[[f]], name=f) 
}


colnames(cd8@meta.data) #from column 34 to 41
FeaturePlot(cd8, features = colnames(cd8@meta.data)[38:45], order=T, 
            cols=rev(brewer.pal(3, 'RdBu')))&NoAxes()&NoLegend()

DotPlot(cd8, features = colnames(cd8@meta.data)[38:45])+theme(axis.text.x = element_text(angle = 60, hjust=1))


#compare with tissue resident signature from Szabo et al, 2019 - Nat Comm -------------------------------
szabo <- read.csv("Szabo_Tr_sig.csv")
head(szabo)
tail(szabo)
szabo=szabo %>% drop_na() %>% select(gene) %>% as.list()

cd8=AddModuleScore(cd8, features = szabo, name='Tr')
FeaturePlot(cd8, features = 'Tr1', cols=rev(brewer.pal(3, 'RdBu')))
VlnPlot(cd8, features = 'Tr1', split.by = 'tissue', split.plot = T, cols=brewer.pal(3, 'Set1')[1:2])
VlnPlot(cd8, features = 'Tr1', group.by = 'new_id', split.by = 'tissue', 
        cols=brewer.pal(3, 'Set1')[1:2], pt.size = 0.0001, split.plot = T)+
  geom_hline(yintercept = 0.1)


### Compare with Tr signature from Savas, Nat Med 2018: ------
savas <- read.csv("Savas_Tr.csv")
savas=savas %>% select(gene) %>% as.list()

cd8=AddModuleScore(cd8, features = savas, name='Tr_Savas')
FeaturePlot(cd8, features = 'Tr_Savas1', cols=rev(brewer.pal(3, 'RdBu')))
VlnPlot(cd8, features = 'Tr_Savas1', split.by = 'tissue', split.plot = T, cols=brewer.pal(3, 'Set1')[1:2])
VlnPlot(cd8, features = 'Tr_Savas1', group.by = 'new_id', 
        pt.size = 0.0001,)+geom_hline(yintercept = 0.1)


################################################################################
##### 4. rename clusters by top gene ----
unique(Idents(tcells))
tcells=RenameIdents(tcells,  'CD4_1'='CD4_CCR7','CD4_2'='CD4_KLRB1','CD4_3'='CD4_NEAT1',
                    'CD4_4'='CD4_FOS', 'CD4_5'='CD4_ISG15','Treg'='CD4_Treg',
                    'CD8_1'='CD8_GZMB', 'CD8_2'='CD8_GNLY','CD8_3'='CD8_GZMK',
                    'CD8_4'='CD8_CCL4','CD8_5'='CD8_GZMH','CD8_6'='CD8_ZNF683', 
                    'CD8_7'='CD8_XCL1' ,'CD8_8'='CD8_CXCL13')
tcells$cell_mk=Idents(tcells)

## SAVE IT! -----
saveRDS(tcells, "~/analysis/data/T_cell_cluster.seu.obj_meta.rds")
saveRDS(cd8, "~/analysis/data/cd8.seu.obj_all_meta.rds")

write.csv(tcells@meta.data, 'meta_tcells_All.csv')
write.csv(cd8@meta.data, 'meta_cd8_subsetFromTcells_All.csv')
