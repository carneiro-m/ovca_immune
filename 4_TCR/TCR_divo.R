
library(divo)
library(corrplot)
library(patchwork)
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)


#load data
tcells <- readRDS("data/T_cell_cluster.seu.obj_meta.rds")#update in Dec 2021
cd8.k=unique(Idents(tcells))[grep('CD8', unique(Idents(tcells)))] %>% as.character()
cd8= subset(tcells, idents=c(cd8.k, 'prolif'))
#RUN DIVO  -----
cd8$cell_mk=factor(cd8$cell_mk, levels=c("CD8_GNLY", "CD8_GZMB","CD8_GZMK","CD8_GZMH" ,
                                    "CD8_CCL4","CD8_ZNF683","CD8_XCL1" ,  "CD8_CXCL13"))

to_divo= cd8@meta.data %>% rownames_to_column(var="barcode") %>% 
  select(c( barcode,orig.ident,cell_mk,  clotype))  %>% 
  drop_na() %>% #filter(Cell.ID %in% unique(data$Cell.ID)[grep("CD", unique(data$Cell.ID))]) %>% #select inecessary informations, obtain barcode ID and remove cells that we do not have TCR-seq information; Filter only clusters in CD4 and CD8
  separate(orig.ident, sep="_", into=c("patient", "tissue")) %>%
  mutate_if(is.character, str_replace_all, pattern = "p09193", replacement = "p09293") %>%#since the label of patient 9293 is wrong in one of the samples, I will fix it 
  filter(patient != 'p09704') %>% #remove patient who data is only for tumor tissue 
  select(c( barcode,  clotype, cell_mk)) %>%
  add_count(clotype, cell_mk,  name="n") %>% #count the frequency of clone in each tissue of a sample/patient
  distinct(clotype, cell_mk, .keep_all=T) %>%
  select(n, cell_mk, clotype) %>% arrange(cell_mk) %>% pivot_wider(names_from='cell_mk', values_from='n')  %>% 
  replace(is.na(.), 0) %>%  column_to_rownames('clotype')  %>% as.matrix()

divo=mh(to_divo, CI = 0.95, graph =F, csv_output = FALSE,
   PlugIn = FALSE, size = 1, saveBootstrap = FALSE)


pdf("/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/Morisita.pdf")
corrplot(divo$Mean, is.corr = F,type='lower',col.lim = c(0,1),  
         col=COL1('Blues', 5), 
         addCoef.col = 'black', tl.col = 'black', tl.srt = 45,cl.ratio = 0.2, 
         method='color')

corrplot(divo$Mean, is.corr = F, method='color',  type='lower', tl.col = 'black', 
         tl.srt = 45, col=COL1('Blues', 5))

dev.off()


#plot with complexheatmap -----
library(ComplexHeatmap) 
Heatmap(divo$Mean, rect_gp = gpar(type = "none"), name='numbers',
        cluster_rows = FALSE, cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(i >= j) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          }
        }, col=hcl.colors(5, 'BluYl'))
Heatmap(divo$Mean, 
        cluster_rows = F, cluster_columns = F, col=hcl.colors(5, 'BluYl'))

od = hclust(dist(divo$Mean))$order
cor_mat = divo$Mean[od, od]
nm = rownames(cor_mat)
col_fun = circlize::colorRamp2(c(-1, 0, 1), c("green", "white", "red"))


Heatmap(cor_mat, name = "correlation", col = col_fun, rect_gp = gpar(type = "none"), 
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, 
                    gp = gpar(col = "grey", fill = NA))
          if(i == j) {
            grid.text(nm[i], x = x, y = y)
          } else if(i >= j) {
            grid.rect(x, y, width, height, gp = gpar(fill = fill, col = fill))
          }
        }, cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE, show_column_names = FALSE)



#RUN DIVO COMPARING BLOOD AND TUMOR IN EACH CLUSTER -----
to_divo_t= data_freq_cluster %>% filter(cell_mk %in% cd8.prolif) %>% #import data calcualted in the file "TCR_plots.R"
           select(clotype,tissue, cell_mk, n) %>% #select only imprtant data
           mutate(sample=paste(cell_mk, tissue, sep='-')) %>% select(-c(cell_mk, tissue)) %>% # create collumn with samples to look at: clusters/tissue
           pivot_wider(names_from='sample', values_from='n')  %>% 
           replace(is.na(.), 0) %>% column_to_rownames('clotype') %>% as.matrix()

  
divo_t=mh(to_divo_t, CI = 0.95, graph =F, csv_output = FALSE,
        PlugIn = FALSE, size = 1, saveBootstrap = FALSE)

corrplot(divo_t$Mean, is.corr = F,type='lower',col.lim = c(0,1),  
         col=COL1('Blues', 5), tl.cex=0.8,
          tl.col = 'black', tl.srt = 45,cl.ratio = 0.2) #addCoef.col = 'black'



#RUN DIVO COMPARING BLOOD AND TUMOR IN EACH leiden cluster (PAGA) DIDINT WORK -----
to_divo_l= tcells@meta.data  %>% rownames_to_column(var="barcode") %>% 
  select(c( barcode,orig.ident,tissue,new_id,  clotype))  %>% 
  drop_na() %>% #filter(Cell.ID %in% unique(data$Cell.ID)[grep("CD", unique(data$Cell.ID))]) %>% #select inecessary informations, obtain barcode ID and remove cells that we do not have TCR-seq information; Filter only clusters in CD4 and CD8
  separate(orig.ident, sep="_", into=c("patient", "tissue")) %>%
  mutate_if(is.character, str_replace_all, pattern = "p09193", replacement = "p09293") %>%#since the label of patient 9293 is wrong in one of the samples, I will fix it 
  filter(patient != 'p09704') %>%
  add_count(new_id, clotype, tissue, name="n") %>% #count the frequency of each clones per cluster
  distinct(clotype,tissue, new_id, .keep_all = T) %>% 
  filter(new_id %in% leiden) %>% #import data calcualted in the file "TCR_plots.R"
  select(clotype,tissue, new_id, n) %>% #select only imprtant data
  mutate(sample=paste(new_id, tissue, sep='-')) %>% select(-c(new_id, tissue)) %>% # create collumn with samples to look at: clusters/tissue
  pivot_wider(names_from='sample', values_from='n') %>% 
  replace(is.na(.), 0)   %>%
  column_to_rownames('clotype') %>% as.matrix()


divo_l=mh(to_divo_l)

corrplot(divo_l$Mean, is.corr = F,type='lower',col.lim = c(0,1),  
         col=COL1('Blues', 5), tl.cex=0.8,
          tl.col = 'black', tl.srt = 45,cl.ratio = 0.2)

#RUN DIVO FOR EACH PATIENT ----
data_type_cluster <- tcells@meta.data  %>% rownames_to_column(var="barcode") %>% 
  select(c( barcode,orig.ident,tissue,cell_mk,  clotype))  %>% 
  drop_na() %>% #filter(Cell.ID %in% unique(data$Cell.ID)[grep("CD", unique(data$Cell.ID))]) %>% #select inecessary informations, obtain barcode ID and remove cells that we do not have TCR-seq information; Filter only clusters in CD4 and CD8
  separate(orig.ident, sep="_", into=c("patient", "tissue")) %>%
  mutate_if(is.character, str_replace_all, pattern = "p09193", replacement = "p09293") %>%#since the label of patient 9293 is wrong in one of the samples, I will fix it 
  filter(patient != 'p09704') %>% #remove patient who data is only for tumor tissue 
  select(c( barcode, patient,tissue, clotype, cell_mk)) %>%
  add_count(patient, clotype, tissue,cell_mk,  name="n") %>% #count the frequency of clone in each tissue of a sample/patient
  select(-barcode) %>% distinct(clotype,tissue,cell_mk,.keep_all = T) %>%  # since we already counted, I will remove the barcode collumn and select unique row for a clonotype within tissue
  group_by(cell_mk) %>% spread(tissue,n ) %>% replace_na(list(blood=0, tumor=0))%>% 
  group_by(cell_mk, clotype,) %>% mutate(type=(ifelse(blood==1 & tumor==0,"b",ifelse(blood > 1 & tumor==0,"B", 
                                                                                     ifelse(blood >= 1 & tumor >=1,"D",ifelse(blood==0 & tumor==1,"t",
                                                                                                                              ifelse(blood==0 & tumor>1,"T","NS"))))))) %>% #classify by tissue expansion
  select(clotype, type, patient)

data_freq_cluster <- tcells@meta.data  %>% 
  select(c( orig.ident,tissue,cell_mk,  clotype))  %>% 
  drop_na() %>% #filter(Cell.ID %in% unique(data$Cell.ID)[grep("CD", unique(data$Cell.ID))]) %>% #select inecessary informations, obtain barcode ID and remove cells that we do not have TCR-seq information; Filter only clusters in CD4 and CD8
  separate(orig.ident, sep="_", into=c("patient", "tissue")) %>%
  mutate_if(is.character, str_replace_all, pattern = "p09193", replacement = "p09293") %>%#since the label of patient 9293 is wrong in one of the samples, I will fix it 
  filter(patient != 'p09704') %>%
  full_join(data_type_cluster, .by=clotype) %>% 
  add_count(cell_mk, clotype, tissue, name="n") %>% #count the frequency of each clones per cluster
  add_count(tissue,cell_mk, name="n_tissue") %>% #count total number of cells in tissue per cluster
  distinct(clotype,tissue, cell_mk, .keep_all = T) %>%  # since we already counted, I will remove the barcode collumn and select unique row for a clonotype within tissue
  mutate(freq=n/n_tissue) %>%  #calculate frequency withiin tissue
  mutate(level_1=ceiling(9/5 * log(n)))%>% mutate(level_2=case_when(level_1 >8 ~8)) %>%#level to choose color in matrix C
  mutate(level=coalesce(level_2, level_1, )) %>% select(-c(level_1, level_2)) %>% #clean up level: cannot be higher then 8 (numer of avialble colors)
  mutate(color=case_when(type =="b" ~ "#F7EF18", type=="t" ~ "#F59110")) %>% #add color for singletons
  mutate(level=recode(level,  `0` = 1))  #mutate(color=replace_na(color, "find")) %>% #space to replace the D color


to_patient= data_freq_cluster %>% filter(cell_mk %in% cd8.prolif) %>% #import data calcualted in the file "TCR_plots.R"
            group_by(patient) %>% group_split()

names(to_patient)= ((data_freq_cluster %>% filter(cell_mk %in% cd8.prolif) %>% #import data calcualted in the file "TCR_plots.R"
  group_by(patient) %>% group_keys())$patient)
  
to_divo_p=  map(to_patient, function (x) {
  x %>% select(clotype,tissue, cell_mk, n) %>% #select only imprtant data
    mutate(sample=paste(cell_mk, tissue, sep='-')) %>% select(-c(cell_mk, tissue)) %>% # create collumn with samples to look at: clusters/tissue
    pivot_wider(names_from='sample', values_from='n')  %>% 
    replace(is.na(.), 0) %>% column_to_rownames('clotype') %>% as.matrix()
  })
  

divo_p=map(to_divo_p, function (x) { mh(x, CI = 0.95, graph =F, csv_output = FALSE,
          PlugIn = FALSE, size = 1, saveBootstrap = FALSE)})

par(mfrow=c(2, 2))
plots_p=map(divo_p, function (x) {corrplot(x[['Mean']], is.corr = F,type='lower',col.lim = c(0,1),  
         col=COL1('Blues', 5), tl.cex=0.8,
          tl.col = 'black', tl.srt = 45,cl.ratio = 0.2)}) #addCoef.col = 'black',




plots_p=map(divo_p, function (x) {corrplot(x[['Mean']], is.corr = F,type='lower',col.lim = c(0,1),  
           col=COL1('Blues', 5), tl.cex=0.8,
           tl.col = 'black', tl.srt = 45,cl.ratio = 0.2)}) #addCoef.col = 'black',

pdf("/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/Morisita_pat.pdf")

par(mfrow=c(2,2))
corrplot(divo_p$p09293[['Mean']], is.corr = F,type='lower',col.lim = c(0,1),  
         col=COL1('Blues', 5), tl.cex=0.5,title = 'p9293',
          tl.col = 'black', tl.srt = 45,cl.ratio = 0.2) #addCoef.col = 'black'

corrplot(divo_p$p09454[['Mean']], is.corr = F,type='lower',col.lim = c(0,1),  
         col=COL1('Blues', 5), tl.cex=0.5,title = 'p9454',
          tl.col = 'black', tl.srt = 45,cl.ratio = 0.2)

corrplot(divo_p$p09808[['Mean']], is.corr = F,type='lower',col.lim = c(0,1),  
         col=COL1('Blues', 5), tl.cex=0.5,title='p9808',
          tl.col = 'black', tl.srt = 45,cl.ratio = 0.2)

corrplot(divo_p$p10329[['Mean']], is.corr = F,type='lower',col.lim = c(0,1),  
         col=COL1('Blues', 5), tl.cex=0.5,title = 'p10329',
          tl.col = 'black', tl.srt = 45,cl.ratio = 0.2)

dev.off()

#RUN DIVO ON SHARED CLONES
to_divo_shared= subset(tcells@meta.data, cell_mk %in% cd8.prolif) %>% rownames_to_column(var="barcode") %>% 
  select(c( barcode,orig.ident,cell_mk,  clotype))  %>% 
  drop_na() %>% #filter(Cell.ID %in% unique(data$Cell.ID)[grep("CD", unique(data$Cell.ID))]) %>% #select inecessary informations, obtain barcode ID and remove cells that we do not have TCR-seq information; Filter only clusters in CD4 and CD8
  separate(orig.ident, sep="_", into=c("patient", "tissue")) %>%
  mutate_if(is.character, str_replace_all, pattern = "p09193", replacement = "p09293") %>%#since the label of patient 9293 is wrong in one of the samples, I will fix it 
  filter(patient != 'p09704') %>% #remove patient who data is only for tumor tissue 
  select(c( barcode,  clotype, cell_mk)) %>%
  mutate_if(is.factor, as.character) %>%
  add_count(clotype, cell_mk,  name="n") %>% #count the frequency of clone in each tissue of a sample/patient
  distinct(clotype, cell_mk, .keep_all=T) %>%
  select(n, cell_mk, clotype) %>% pivot_wider(names_from='cell_mk', values_from='n')  %>% 
  replace(is.na(.), 0) %>% 
  filter(clotype %in% clo) %>% 
  column_to_rownames('clotype')  %>% as.matrix()

divo_shared=mh(to_divo_shared, CI = 0.95, graph =F, csv_output = FALSE,
        PlugIn = FALSE, size = 1, saveBootstrap = FALSE)

corrplot(divo_shared$Mean, is.corr = F,type='lower',col.lim = c(0,1),  
         col=COL1('Blues', 5), 
         addCoef.col = 'black', tl.col = 'black', tl.srt = 45,cl.ratio = 0.2)

ggplot(data_tcr_k %>% filter(cell_mk %in% cd8) %>% filter(type=='D')
       ,aes( axis2=freq_b,  axis1=freq_t))+
  geom_alluvium(aes(fill=cell_mk ))+ 
  geom_stratum()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void() + theme(axis.text.x = element_blank()) +
  ggtitle('' )

#RUN DIVO in leiden  -----
Idents(cd8)='leiden'
to_divo_leiden= subset(cd8, idents=c('k_23', 'k_18'), invert=T)@meta.data %>% rownames_to_column(var="barcode") %>% 
  select(c( barcode,orig.ident,leiden,  clotype))  %>% 
  drop_na() %>% #filter(Cell.ID %in% unique(data$Cell.ID)[grep("CD", unique(data$Cell.ID))]) %>% #select inecessary informations, obtain barcode ID and remove cells that we do not have TCR-seq information; Filter only clusters in CD4 and CD8
  separate(orig.ident, sep="_", into=c("patient", "tissue")) %>%
  mutate_if(is.character, str_replace_all, pattern = "p09193", replacement = "p09293") %>%#since the label of patient 9293 is wrong in one of the samples, I will fix it 
  filter(patient != 'p09704') %>% #remove patient who data is only for tumor tissue 
  select(c( barcode,  clotype, leiden)) %>%
  add_count(clotype, leiden,  name="n") %>% #count the frequency of clone in each tissue of a sample/patient
  distinct(clotype, leiden, .keep_all=T) %>%
  select(n, leiden, clotype) %>% arrange(leiden) %>% pivot_wider(names_from='leiden', values_from='n')  %>% 
  replace(is.na(.), 0) %>%  column_to_rownames('clotype')  %>% as.matrix()

divo_leiden=mh(to_divo_leiden, CI = 0.95, graph =F, csv_output = 'divo_leiden',
        PlugIn = FALSE, size = 1, saveBootstrap = FALSE)


pdf("/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/Morisita_leiden.pdf", width=16, height=10)
corrplot(divo_leiden$Mean, is.corr = F, method='color',  type='lower', tl.col = 'black', 
         tl.srt = 45, col=COL1('Blues', 5), addCoef.col = 'black',)
dev.off()


## check shared ZNF and GZMK -----
clones_znf=cd8@meta.data %>% rownames_to_column(var="barcode") %>% 
  select(c( barcode,leiden,  clotype))  %>% 
  drop_na() %>% #filter(Cell.ID %in% unique(data$Cell.ID)[grep("CD", unique(data$Cell.ID))]) %>% #select inecessary informations, obtain barcode ID and remove cells that we do not have TCR-seq information; Filter only clusters in CD4 and CD8
  filter(leiden %in% c('k_8', 'k_10', 'k_11', 'k_7')) %>% pull(clotype) %>% as.character()

clones_gzmk=cd8@meta.data %>% rownames_to_column(var="barcode") %>% 
  select(c( barcode,leiden,  clotype))  %>% 
  drop_na() %>% #filter(Cell.ID %in% unique(data$Cell.ID)[grep("CD", unique(data$Cell.ID))]) %>% #select inecessary informations, obtain barcode ID and remove cells that we do not have TCR-seq information; Filter only clusters in CD4 and CD8
  filter(leiden %in% c('k_22', 'k_5')) %>% pull(clotype) %>% as.character()

clones_all_fil=cd8@meta.data %>% rownames_to_column(var="barcode") %>% 
  select(c( barcode,leiden,  clotype))

list_4upset_gzmk = list()
for (k in levels(cd8$leiden)){
  list_4upset_gzmk[[k]]= clones_all_fil %>% filter(clotype %in% clones_gzmk) %>%
    filter(leiden==k) %>% pull(clotype) %>% as.character()
}

list_4upset_znf = list()
for (k in levels(cd8$leiden)){
  list_4upset_znf[[k]]= clones_all_fil %>% filter(clotype %in% clones_znf) %>%
    filter(leiden==k) %>% pull(clotype) %>% as.character()
}


library(UpSetR)
upset(fromList(list_4upset_gzmk), nsets = 24, order.by='freq',
      nintersects=27,mb.ratio = c(0.6 ,0.4))

upset(fromList(list_4upset_znf), nsets = 24,
      nintersects=30,mb.ratio = c(0.6 ,0.4))


