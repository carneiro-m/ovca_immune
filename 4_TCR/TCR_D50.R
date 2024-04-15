### Diversity index

#D50 is the percent of dominant T cell clones that account for the cumulative 50% of the total CDR3s counted in the sample. 
#The mathematical formulation of D50 is defined as follows:
#D50=(No.of uCDR3 that make up 50% of the total reads×100)/No. of uCDR3s

#Diversity 50 (D50) was calculated on Excel as the fraction of dominant clones that account for the cumulative 50% of the total paired CDR3s identified in each UMAP cluster.
#CDR3 similarity (TCR sharing/clonal overlap) was calculated using the Morisita-horn overlap index by using the divo package (Rempala and Seweryn, 2013)
#tcr data - arrange cells for D50
### IMPORTANT: calculate for each patient???



### Set colors
nb.cols <- 20
mycolors <- colorRampPalette(brewer.pal("Set3",n=9))(nb.cols)
names(mycolors) <- levels(tcells)

cd8.col=c(mycolors[7:14])
names(cd8.col)

#Set color condition
cols.tissue <- c('#C0392B', '#00838F')
names(cols.tissue) <- c("blood", "tumor")

#calculate ----
data_type_cluster <- tcells@meta.data  %>% rownames_to_column(var="barcode") %>% 
  select(c( barcode,orig.ident,tissue,cell_mk,  clotype))  %>% 
  drop_na() %>% #filter(Cell.ID %in% unique(data$Cell.ID)[grep("CD", unique(data$Cell.ID))]) %>% #select inecessary informations, obtain barcode ID and remove cells that we do not have TCR-seq information; Filter only clusters in CD4 and CD8
  separate(orig.ident, sep="_", into=c("patient", "tissue")) %>%
  mutate_if(is.character, str_replace_all, pattern = "p09193", replacement = "p09293") %>%#since the label of patient 9293 is wrong in one of the samples, I will fix it 
  filter(patient != 'p09704') %>% #remove patient who data is only for tumor tissue 
  select(c( barcode, patient,tissue, clotype, cell_mk)) %>%
  add_count(patient, clotype, tissue,  name="n") %>% #count the frequency of clone in each tissue of a sample/patient
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


data_freq_k_tumor=filter(data_freq_cluster, tissue=="tumor")
data_freq_k_blood=filter(data_freq_cluster, tissue=="blood")

#### D50
prep.tumor=data_freq_k_tumor %>% 
      group_by(cell_mk) %>% 
      arrange(desc(n), .by_group=T) %>% #order from greater to smaller number of cells of a clone within a cluster
      group_by(cell_mk) %>% 
      mutate(test=cumsum(n)) %>% #each row  of test will be the n+ the sum of previous clones, like that we can filter by amount of cells the clones represent
      filter(test<= n_tissue/2, .preserve=T) %>% #select only the dominant clones, which number of cells is equivalent to half of the cells in the cluster
      count(cell_mk) #now we have the mumber of clones that represent the half of the clones but we still need to calcuate their frequency in total number of cells to obtain the D50

d50.tumor=prep.tumor %>% 
    left_join(data_freq_k_tumor %>% select(cell_mk, n_tissue)) %>%
    mutate(d50_tumor=round(n/n_tissue, digits=2)) %>%
    distinct(cell_mk, .keep_all=T)

#blood
prep.blood=data_freq_k_blood %>% 
  group_by(cell_mk) %>% 
  arrange(desc(n), .by_group=T) %>% 
  group_by(cell_mk) %>% 
  mutate(test=cumsum(n)) %>% 
  filter(test<= n_tissue/2, .preserve=T) %>% 
  count(cell_mk)

d50.blood=prep.blood %>% 
  left_join(data_freq_k_blood %>% select(cell_mk, n_tissue)) %>%
  mutate(d50_blood=round(n/n_tissue, digits=2)) %>%
  distinct(cell_mk, .keep_all=T)

#final
d50.final=d50.blood %>% select( cell_mk,d50_blood) %>%
    left_join(d50.tumor %>% select(cell_mk, d50_tumor))

pdf("/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/d50.pdf", width=10)
d50.final %>% mutate(ratio=(d50_tumor/d50_blood)) %>% 
        ggplot(aes(cell_mk, ratio))+geom_col()+
        theme(axis.text.x = element_text(angle=60, hjust=1))

d50.final %>% filter(cell_mk %in% cd8.prolif) %>%
        pivot_longer(cols=2:3, names_to='tissue',values_to='d50') %>% 
        ggplot(aes(cell_mk, d50,fill=tissue))+geom_col(position='dodge')+
        theme(axis.text.x = element_text(angle=60, hjust=1, color=cd8.col))+
        scale_fill_manual(values =c("lightgrey", "#00838F"))
dev.off()



## D50 for each clusters without separating tissues -----
data_freq_cluster <- tcells@meta.data  %>% 
  select(c( orig.ident,cell_mk,  clotype))  %>% 
  drop_na() %>% #filter(Cell.ID %in% unique(data$Cell.ID)[grep("CD", unique(data$Cell.ID))]) %>% #select inecessary informations, obtain barcode ID and remove cells that we do not have TCR-seq information; Filter only clusters in CD4 and CD8
  separate(orig.ident, sep="_", into=c("patient", "tissue")) %>%
  mutate_if(is.character, str_replace_all, pattern = "p09193", replacement = "p09293") %>%#since the label of patient 9293 is wrong in one of the samples, I will fix it 
  filter(patient != 'p09704') %>%
  select(-tissue) %>%
  add_count( cell_mk, clotype, name="n") %>% #count the frequency of each clones per cluster
  add_count(cell_mk, name="n_k") %>% #count total number of cells in a cluster
  distinct(clotype, cell_mk, .keep_all = T) %>%  # since we already counted, I will remove the barcode collumn and select unique row for a clonotype within tissue
  mutate(freq=n/n_k)   #calculate frequency withiin tissue

#prep D50  
prep=data_freq_cluster %>% 
  group_by(cell_mk) %>% 
  arrange(desc(n), .by_group=T) %>% #order from greater to smaller number of cells of a clone within a cluster
  #group_by(cell_mk) %>% 
  mutate(test=cumsum(n)) %>% #each row  of test will be the n+ the sum of previous clones, like that we can filter by amount of cells the clones represent
  filter(test<= n_k/2, .preserve=T) %>% #select only the dominant clones, which number of cells is equivalent to half of the cells in the cluster
  count(cell_mk) #now we have the mumber of clones that represent the half of the clones but we still need to calcuate their frequency in total number of cells to obtain the D50

d50=prep %>% 
  left_join(data_freq_cluster %>% select(cell_mk, n_k)) %>%
  mutate(d50=round(n/n_k, digits=2)) %>%
  distinct(cell_mk, .keep_all=T)


(d50 %>% filter(cell_mk %in% cd8.prolif) %>%  
        ggplot()+geom_col(aes(cell_mk, d50), fill='lightgray'))/
  
  (d50.final %>% filter(cell_mk %in% cd8.prolif) %>%
     pivot_longer(cols=2:3, names_to='tissue',values_to='d50') %>% 
     ggplot(aes(cell_mk, d50,color=tissue))+geom_point()+
     scale_color_manual(values =c("lightgrey", "#00838F")))/

(d50 %>% left_join(d50.final) %>% filter(cell_mk %in% cd8.prolif) %>% ggplot() +
        geom_col(aes(cell_mk, d50), fill='lightgray')+
        geom_point(aes(cell_mk, d50_blood), color='red')+
        geom_point(aes(cell_mk, d50_tumor), color="#00838F"))

#they are correct! lets print it 
pdf("/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/d50_final.pdf", width=5, height = 5)
d50 %>% left_join(d50.final) %>% filter(cell_mk %in% cd8.prolif) %>% ggplot() +
  geom_col(aes(cell_mk, d50, width=0.75), fill='lightgray')+
  geom_point(aes(cell_mk, d50_blood), color='red', size=2)+
  geom_point(aes(cell_mk, d50_tumor), color="#00838F", size=2)+theme_bw()+theme(axis.title.x = element_blank(),
                                                                       axis.text.x = element_text(angle = 60, hjust = 1),
                                                                       axis.line = element_line(color='darkgray'),
                                                                       panel.border = element_blank(),
                                                                       panel.grid = element_blank()
                                                                        )
dev.off()

pdf("/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/d50_final_2.pdf", width=5, height = 6)
d50 %>% filter(cell_mk %in% names(cd8.col)) %>% ggplot() +
  geom_col(aes(cell_mk, d50, width=0.75), fill='lightgray')+
  theme_bw()+theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1),
          axis.line = element_line(color='darkgray'),
          panel.border = element_blank(),
          panel.grid = element_blank()
  )
dev.off()

## D50 for each clusters without separating tissues - LEIDEN -----
data_freq_cluster <- cd8@meta.data  %>% 
  select(c( orig.ident,leiden,  clotype))  %>% 
  drop_na() %>% #filter(Cell.ID %in% unique(data$Cell.ID)[grep("CD", unique(data$Cell.ID))]) %>% #select inecessary informations, obtain barcode ID and remove cells that we do not have TCR-seq information; Filter only clusters in CD4 and CD8
  separate(orig.ident, sep="_", into=c("patient", "tissue")) %>%
  mutate_if(is.character, str_replace_all, pattern = "p09193", replacement = "p09293") %>%#since the label of patient 9293 is wrong in one of the samples, I will fix it 
  filter(patient != 'p09704') %>%
  select(-tissue) %>%
  add_count( leiden, clotype, name="n") %>% #count the frequency of each clones per cluster
  add_count(leiden, name="n_k") %>% #count total number of cells in a cluster
  distinct(clotype, leiden, .keep_all = T) %>%  # since we already counted, I will remove the barcode collumn and select unique row for a clonotype within tissue
  mutate(freq=n/n_k)   #calculate frequency withiin tissue

#prep D50  
prep=data_freq_cluster %>% 
  group_by(leiden) %>% 
  arrange(desc(n), .by_group=T) %>% #order from greater to smaller number of cells of a clone within a cluster
  #group_by(leiden) %>% 
  mutate(test=cumsum(n)) %>% #each row  of test will be the n+ the sum of previous clones, like that we can filter by amount of cells the clones represent
  filter(test<= n_k/2, .preserve=T) %>% #select only the dominant clones, which number of cells is equivalent to half of the cells in the cluster
  count(leiden) #now we have the mumber of clones that represent the half of the clones but we still need to calcuate their frequency in total number of cells to obtain the D50

d50=prep %>% 
  left_join(data_freq_cluster %>% select(leiden, n_k)) %>%
  mutate(d50=round(n/n_k, digits=2)) %>%
  distinct(leiden, .keep_all=T)



#they are correct! lets print it 
pdf("/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/d50_final_leiden.pdf", width=8, height = 4)
d50 %>%  ggplot() +
  geom_col(aes(leiden, d50, width=0.75), fill='lightgray')+
  theme_bw()+theme(axis.title.x = element_blank(),
                   axis.text.x = element_text(angle = 60, hjust = 1),
                   axis.line = element_line(color='darkgray'),
                   panel.border = element_blank(),
                   panel.grid = element_blank()
  )
dev.off()
### D50 patient level  ----
patient.list=list()
prep.tumor.p=list()
d50.tumor.p=list()
prep.blood.p=list()
d50.blood.p=list()
d50.p=list()
patients=c('p09293', 'p09454', 'p09808', 'p10329')
names(patients)=patients


for (p in patients){
    prep.tumor.p[[p]]=data_freq_k_p_tumor %>%
    filter(patient==p) %>%
    group_by(cell_mk) %>% 
    arrange(desc(n), .by_group=T) %>% #order from greater to smaller number of cells of a clone within a cluster
    group_by(cell_mk) %>% 
    mutate(test=cumsum(n)) %>% #each row  of test will be the n+ the sum of previous clones, like that we can filter by amount of cells the clones represent
    filter(test<= n_tissue/2, .preserve=T) %>% #select only the dominant clones, which number of cells is equivalent to half of the cells in the cluster
    count(cell_mk) #now we have the mumber of clones that represent the half of the clones but we still need to calcuate their frequency in total number of cells to obtain the D50
  
  d50.tumor.p[[p]]=prep.tumor.p[[p]] %>% 
    left_join(data_freq_k_p_tumor %>% filter(patient==p) %>% select(cell_mk, n_tissue)) %>%
    mutate(d50_tumor=n/n_tissue) %>%
    distinct(cell_mk, .keep_all=T)
  
  #blood
  prep.blood.p[[p]]=data_freq_k_p_blood %>% 
    filter(patient==p) %>%
    group_by(cell_mk) %>% 
    arrange(desc(n), .by_group=T) %>% 
    group_by(cell_mk) %>% 
    mutate(test=cumsum(n)) %>% 
    filter(test<= n_tissue/2, .preserve=T) %>% 
    count(cell_mk)
  
  d50.blood.p[[p]]=prep.blood.p[[p]] %>% 
    left_join(data_freq_k_p_blood %>% filter(patient==p) %>% select(cell_mk, n_tissue)) %>%
    mutate(d50_blood=n/n_tissue) %>%
    distinct(cell_mk, .keep_all=T)
  
  #final
  d50.p[[p]]=d50.blood.p[[p]] %>% select( cell_mk,d50_blood) %>%
    left_join(d50.tumor.p[[p]] %>% select(cell_mk, d50_tumor)) %>%
    mutate(ratio=d50_tumor/d50_blood)
  patient.list[[p]]=d50.p[[p]]
  
}
 
 
 d50.pl.pat=list()
 for (p in patients){
   d50.pl.pat[[p]]=d50.p[[p]] %>% filter(!cell_mk %in% c('uk1', 'uk2', 'myl', 'ukB')) %>%
     pivot_longer(cols=2:3, names_to='tissue',values_to='d50') %>% 
     ggplot(aes(cell_mk, d50,fill=tissue))+geom_col(position='dodge')+
     theme(axis.text.x = element_text(angle=60, hjust=1))+ggtitle(p)
 }
 
 wrap_plots(d50.pl.pat)
 
 d50.pl.pat.2=list()
 for (p in patients){
   d50.pl.pat.2[[p]]= d50.p[[p]] %>% filter(!cell_mk %in% c('uk1', 'uk2', 'myl', 'ukB')) %>%
    mutate(ratio=(d50_tumor/d50_blood)) %>% 
    ggplot(aes(cell_mk, ratio))+geom_col()+
    theme(axis.text.x = element_text(angle=60, hjust=1))+ggtitle(p)
 }
 
 wrap_plots(d50.pl.pat.2)

cd8=levels(d50.p$p09293$cell_mk)[grep('CD8', levels(d50.p$p09293$cell_mk))]
  
pdf("/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/d50_ratio.pdf", height = 4)
bind_rows(d50.p, .id = 'patient') %>%  filter(cell_mk %in% cd8.prolif) %>%
  ggplot(aes(cell_mk, ratio))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes( shape=patient), position=position_jitter(0.2))+
  theme_bw()+
  scale_colour_manual(values=met.brewer('Egypt'))+
  theme(axis.text.x=element_text(color=cd8.col, angle = 60, hjust=1), axis.title.x = element_blank())+
  geom_hline(yintercept=1, linetype="dashed", color = "black")
dev.off()


pdf("/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/d50.pdf", width=5, height=6)
d50 %>% mutate(ratio=(d50_blood/d50_tumor)) %>% 
  ggplot(aes(cell_mk, ratio))+geom_col()+
  theme(axis.text.x = element_text(angle=60, hjust=1))


d50 %>% filter(cell_mk %in% cd8.prolif) %>%
  pivot_longer(cols=2:3, names_to='tissue',values_to='d50') %>% 
  ggplot(aes(cell_mk, d50,fill=tissue))+geom_col(position='dodge')+
  theme_light()+
  theme(axis.text.x = element_text(angle=60, hjust=1, color=cd8.col), 
        axis.title.x=element_blank(), legend.position = 'bottom')+
  scale_fill_manual(values = c('lightgrey', "#00838F" ))

pdf("/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/d50_ratio.pdf", width=4, height=4)
bind_rows(d50.p, .id = 'patient') %>%  filter(cell_mk %in% cd8.prolif) %>%
  ggplot(aes(cell_mk, ratio))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes( shape=patient), position=position_jitter(0.2))+
  theme_light()+
  scale_colour_manual(values=met.brewer('Egypt'))+
  theme(axis.text.x=element_text(color=cd8.col, angle = 60, hjust=1), 
        axis.title.x = element_blank(), legend.position = 'bottom')+
  geom_hline(yintercept=1, linetype="dashed", color = "black")

dev.off()


pdf("/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/d50_patients.pdf", width=8, height=4)
bind_rows(d50.p, .id = 'patient') %>%  filter(cell_mk %in% cd8.prolif) %>%
  pivot_longer(names_to='tissue', cols=starts_with('d50')) %>%
  ggplot(aes(cell_mk, value))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes( cell_mk, value,shape=patient))+#, position=position_jitter(0.2))+
  theme_light()+
  scale_colour_manual(values=met.brewer('Egypt'))+
  theme(axis.text.x=element_text( angle = 60, hjust=1), #color=cd8.col
        axis.title.x = element_blank(), legend.position = 'bottom')+
  facet_wrap(~ tissue, scales = "free_y", ncol = 5) 

  dev.off()



prep.geral=data_tcr_k %>% 
  group_by(cell_mk) %>% 
  arrange(desc(n), .by_group=T) %>% #order from greater to smaller number of cells of a clone within a cluster
  group_by(cell_mk) %>% 
  mutate(test=cumsum(n)) %>% #each row  of test will be the n+ the sum of previous clones, like that we can filter by amount of cells the clones represent
  filter(test<= n_tissue/2, .preserve=T) %>% #select only the dominant clones, which number of cells is equivalent to half of the cells in the cluster
  count(cell_mk) #now we have the mumber of clones that represent the half of the clones but we still need to calcuate their frequency in total number of cells to obtain the D50

d50.geral=prep.tumor %>% 
  left_join(data_freq_k_tumor %>% select(cell_mk, n_tissue)) %>%
  mutate(d50_tumor=round(n/n_tissue, digits=2)) %>%
  distinct(cell_mk, .keep_all=T)





## D50 for each clusters without separating tissues PATIET LEVEL----
data_freq_cluster_patient=list()
prep.p2=list()
d50.p2=list()
patients=c('p09293', 'p09454', 'p09808', 'p10329')
names(patients)=patients

for(p in names(patients)){
  data_freq_cluster_patient[[p]] <- tcells@meta.data  %>% 
    select(c( orig.ident,cell_mk,  clotype))  %>% 
    drop_na() %>% #filter(Cell.ID %in% unique(data$Cell.ID)[grep("CD", unique(data$Cell.ID))]) %>% #select inecessary informations, obtain barcode ID and remove cells that we do not have TCR-seq information; Filter only clusters in CD4 and CD8
    separate(orig.ident, sep="_", into=c("patient", "tissue")) %>%
    mutate_if(is.character, str_replace_all, pattern = "p09193", replacement = "p09293") %>%#since the label of patient 9293 is wrong in one of the samples, I will fix it 
    filter(patient == p) %>%
    select(-tissue) %>%
    add_count( cell_mk, clotype, name="n") %>% #count the frequency of each clones per cluster
    add_count(cell_mk, name="n_k") %>% #count total number of cells in a cluster
    distinct(clotype, cell_mk, .keep_all = T) %>%  # since we already counted, I will remove the barcode collumn and select unique row for a clonotype within tissue
    mutate(freq=n/n_k)   #calculate frequency withiin tissue
  
}

#prep D50  
for(p in names(patients)){
  prep.p2[[p]]=data_freq_cluster_patient[[p]] %>% 
    group_by(cell_mk) %>% 
    arrange(desc(n), .by_group=T) %>% #order from greater to smaller number of cells of a clone within a cluster
    #group_by(cell_mk) %>% 
    mutate(test=cumsum(n)) %>% #each row  of test will be the n+ the sum of previous clones, like that we can filter by amount of cells the clones represent
    filter(test<= n_k/2, .preserve=T) %>% #select only the dominant clones, which number of cells is equivalent to half of the cells in the cluster
    count(cell_mk) #now we have the mumber of clones that represent the half of the clones but we still need to calcuate their frequency in total number of cells to obtain the D50
}

for (p in names(patients)){
  d50.p2[[p]]=prep.p2[[p]] %>% 
    left_join(data_freq_cluster_patient[[p]] %>% select(cell_mk, n_k)) %>%
    mutate(d50=round(n/n_k, digits=2)) %>%
    distinct(cell_mk, .keep_all=T)
}

pdf("/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/d50_all_pat.pdf", width=4, height=4)
bind_rows(d50.p2, .id = 'patient') %>%  filter(cell_mk %in% cd8.k) %>%
  ggplot(aes(cell_mk, d50))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes( shape=patient), position=position_jitter(0.2))+
  theme_bw()+
  #scale_colour_manual(values=cd8.col)+
  theme(axis.text.x=element_text( angle = 60, hjust=1), axis.title.x = element_blank(), legend.position='top')

dev.off()


