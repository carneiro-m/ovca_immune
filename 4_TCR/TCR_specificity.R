### TCR specificity
# I have downloaded the genes TCR informations from 3 databases: 
#PIRD (accessed in May 10, 2021 - https://db.cngb.org/pird/tbadb/) 
#MPAS-TCR (accesed in May 10, 2021 - http://friedmanlab.weizmann.ac.il/McPAS-TCR/) - all downloaded
#VDJdb (accessed in May 10, 2021 - https://vdjdb.cdr3.net/search) - selected: human, alpha and beta chains

#updated the analysis with upddates in those database. - Acessed in Thu Jun 9 2:10PM


#load libraries ------
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(readr)
library(stringr)
library(ggplot2)  
library(RColorBrewer)
library(immunarch)

#load data -----
pird <- read_csv("data/TBAdb.csv")#pird <- read_excel("data/TBAdb.xlsx") 
mcpas=read.csv("data/McPAS-TCR.csv")
vdjdb=read_tsv('data/SearchTable-2022-06-09 18_04_22.322.tsv')
dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/TBAdb.xlsx", 
       "tbadb", .species = "Homo Sapiens", .chain = c("TRB", "TRA-TRB"), .pathology = "CMV")

#explore datasets -----
str(pird)
str(mcpas)
str(vdjdb)

#select info (id, cdr3 from alpha and beta, species, disease)
pird= pird %>% 
  rownames_to_column(var='id') %>%
  mutate(id=paste('pird', id, sep='_')) %>%
  rename(TRA=CDR3.alpha.aa, TRB=CDR3.beta.aa, Pathology=Disease.name) %>%
  pivot_longer(cols=c('TRA', 'TRB'), names_to='Gene', values_to='CDR3') %>%
  select(id, ICDname, Pathology, Category, Gene, Locus, CDR3) %>%
  filter(Pathology!='-') #some rows dont have the information of pathology

mcpas= mcpas %>% 
  rownames_to_column(var='id') %>%
  mutate(id=paste('mcpas', id, sep='_')) %>%
  filter(Species=='Human') %>%
  rename(TRA=CDR3.alpha.aa, TRB=CDR3.beta.aa) %>%
  pivot_longer(cols=c('TRA', 'TRB'), names_to='Gene', values_to='CDR3') %>%
  select(id, Pathology, Category, Gene, CDR3) 

vdjdb = vdjdb %>%
  rownames_to_column(var='id') %>%
  mutate(id=paste('vdjdb', id, sep='_')) %>%
  rename(Pathology=`Epitope species`) %>%
  select(id,Pathology, Gene, CDR3)%>%
  filter( Pathology!='HomoSapiens') #what is that?? I am excluding 

#join the data.frames -----
tcr.base=bind_rows(pird, mcpas, vdjdb) %>% distinct(CDR3, .keep_all=T) %>%
  filter(CDR3 != '-') %>%  #remove this row without CDR3 sequency
  add_row(Pathology='noTCR-seq') #add row with NA in CDR3 to be able to track the cells without TCR-seq afterwards

table(tcr.base$Locus)
tcr.base %>% filter(Pathology=='HomoSapiens')


#load our dataset -----
##load data
tcells <- readRDS("data/T_cell_cluster.seu.obj_meta.rds")#update in Dec 2021
tcells$barcode=rownames(tcells@meta.data)

#add leiden info to tcells obj -----
#leiden clusters
#leiden <- read.csv("~/analysis/results/tcell_cd8s_bl_tu_leiden.csv") 
#leiden$barcode=leiden$Unnamed..0
#leiden=leiden %>% mutate(leiden=paste('k', leiden, sep='_'))
#str(leiden)
#n_k=length(unique(leiden$leiden))

#tcells@meta.data = tcells@meta.data %>% left_join(leiden[c('leiden', 'barcode')], by='barcode') %>%
#  mutate(new_id=coalesce(leiden,Cell.ID )) %>% column_to_rownames('barcode')
#levels_leiden=paste(rep('k',n_k), seq(from=0, to=n_k-1), sep='_')
#tcells$leiden=factor(tcells$leiden, levels=levels_leiden)
#tcells$new_id=factor(tcells$new_id, levels = c("CD4_1","CD4_2","CD4_3","CD4_4", 
                                               #"CD4_5","Treg",levels_leiden,"HS","myl", 
                                               # "ukB", "uk1",  "uk2" ))
#table(tcells$new_id)


# looking at same TCRs in our dataset and bank of TCRs ------
ov.data= tcells@meta.data %>% select(clotype, cdr3s_aa, cell_mk, patient) %>%
  mutate(cluster='cell_mk') %>% select(-cell_mk) %>% #rename the cell_mk collumn to make it easier
  rownames_to_column(var='barcode') %>%
  drop_na(cdr3s_aa) %>% #remove cells without TCR info
  separate_rows(cdr3s_aa, sep=';') %>% #separate the chains
  separate(cdr3s_aa,c('Gene', 'CDR3'), sep=':') #put the chain info in another column to be able to match with CDR3 from tcr.base



#add the Pathology to clonotype in ovarian cancer dataset
ov.ident= ov.data %>% left_join(tcr.base, by='CDR3') 
ov.ident %>% select(Pathology, CDR3, clotype, barcode) %>% head()


#since some pathologies are gropued with different names (because it came from diff database), I will organize it
unique(ov.ident$Pathology)
hiv=unique(ov.ident$Pathology)[grep('HIV', unique(ov.ident$Pathology))]
cmv=unique(ov.ident$Pathology)[grep('CMV', unique(ov.ident$Pathology))]
ebv=unique(ov.ident$Pathology)[grep('EBV', unique(ov.ident$Pathology))]
hcv=c('HCV', 'Hepatitis C virus')
influenza=unique(ov.ident$Pathology)[grep('Influ', unique(ov.ident$Pathology))]
tuberculosis=unique(ov.ident$Pathology)[grep('ubercu', unique(ov.ident$Pathology))]
covid=c('SARS-CoV-2', 'COVID-19') 
yfv=c(unique(ov.ident$Pathology)[grep('Yello', unique(ov.ident$Pathology))], 'YFV')
tumor=c('Malignant neoplasm of brain', "Lymphoma","Neoantigen" , "TAA",
        "Melanoma","Cervical Cancer","Carcinoma", "Breast Cancer" , "Merkel cell carcinoma", "Colorectal cancer" )
virus=c('CMV', "Influenza" ,"EBV", "YFV","HIV","SARS-CoV-2","HCV" ,"HTLV-I-Associated Myelopathy/Tropical Spastic Paraparesis (HAM/TSP)",
        "MCPyV","Hepatitis E virus infection (cHEV)" ,"HTLV-1","Chronic viral hepatitis B"  ,
        "HPV", "Human herpes virus 1" ,"Herpes simplex virus 2 (HSV2)" )

auto.immune=c("Systemic lupus erythematosus(SLE)" , "Inflammatory bowel disease (IBD)","Psoriatic arthritis",
              "Alzheimer's disease","Celiac disease","Parkinson disease", "Diabetes Type 1",
              "Multiple sclerosis", "Rheumatoid Arthritis (RA)" , "Renal allograft" )

others=c("tuberculosis","Rasmussen encephalitis","DENV1","DENV3/4" ,
          "Sarcoidosis","Aplastic anemia", "Inclusion body myositis","Chronic inflammatory periodontal disease",
         "Sézary disease" , "Polymyositis","Polyradiculitis"  ,"Ataxia telangiectasia(ATM gene)" ,
          "Allergy","Encounter for immunization")


# change the names duplicated
ov.ident=ov.ident %>%
  mutate(Pathology= case_when(grepl(hiv[1], Pathology) ~ 'HIV',
                              grepl('Human immunodeficiency', Pathology) ~ 'HIV',
                              grepl(cmv[1], Pathology) ~ 'CMV',
                              grepl(cmv[2], Pathology) ~ 'CMV',
                              grepl(ebv[1], Pathology) ~ 'EBV',
                              grepl(ebv[2], Pathology) ~ 'EBV',
                              grepl(hcv[1], Pathology) ~ 'HCV',
                              grepl(hcv[2], Pathology) ~ 'HCV',
                              grepl(influenza[1], Pathology) ~ 'Influenza',
                              grepl(influenza[2], Pathology) ~ 'Influenza',
                              grepl(tuberculosis[1], Pathology) ~ 'tuberculosis',
                              grepl(tuberculosis[2], Pathology) ~ 'tuberculosis',
                              grepl(tuberculosis[3], Pathology) ~ 'tuberculosis',
                              grepl(tuberculosis[4], Pathology) ~ 'tuberculosis',
                              grepl(covid[1], Pathology) ~ 'SARS-CoV-2',
                              grepl(covid[2], Pathology) ~ 'SARS-CoV-2',
                              grepl(yfv[1], Pathology) ~ 'YFV',
                              grepl(yfv[2], Pathology) ~ 'YFV',
                              grepl(yfv[3], Pathology) ~ 'YFV',
                              grepl("Tumor", Pathology) ~ "TAA",
                              grepl("Inflammatory bowel", Pathology) ~ 'Inflammatory bowel disease (IBD)',
                              TRUE ~ Pathology
                              
  ))

to_change=data.frame(Pathology=c( tumor, virus, auto.immune, others), 
                     Pathology_groups=c(rep('tumor', 10), rep('virus', 15),
                              rep('autoimmune', 10), rep('others', 14)))

#add the groups of pathology (tumor, virus, )
ov.ident=ov.ident %>% left_join(to_change, by='Pathology') %>% 
         distinct(barcode, .keep_all=T)
ov.ident$Pathology_groups=  replace_na(ov.ident$Pathology_groups, 'unknown-TCR')
write.csv(ov.ident, 'TCR_specificity_upda.csv')

#check if all cells were preserved
length(unique(ov.ident$barcode)) 
summary(is.na(tcells@meta.data$cdr3s_aa)) # yes!


tcells@meta.data$barcode=rownames(tcells@meta.data)
tcells@meta.data=tcells@meta.data %>% left_join(ov.ident %>% select(Pathology_groups, Pathology, barcode), by='barcode')
tcells$Pathology_groups=factor(tcells$Pathology_groups, levels = 
                        c('tumor', 'virus', 'autoimmune',
                          'others', 'unknown-TCR'))



#Plot ----
nb.cols <- 20
mycolors <- colorRampPalette(brewer.pal("Set3",n=9))(nb.cols)
names(mycolors) <- levels(tcells)

cd8.id=levels(tcells$new_id)[7:30]
cols.path=c('#00838F', brewer.pal(3, 'RdPu'),'grey' )
names(cols.path)=levels(tcells$Pathology_groups)

cd8.k=levels(unique(Idents(tcells)))[grep('CD8',levels(Idents(tcells)) )]
names(cd8.k)=cd8.k

# Seurat clusters - Cell.ID
pdf('/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/specificity_cell_mk.pdf', width=10, height=7)
table(Pathology=tcells$Pathology_groups, cluster=tcells$cell_mk) %>% 
  data.frame()%>% #filter(cluster %in% cd8.k) %>%
  ggplot(aes(cluster, Freq, fill=Pathology))+geom_col(position='fill')+ #position='fill' for prportion instead of number of cells
  coord_flip()+
  scale_fill_manual(values = cols.path)+
  ggtitle('')+ylab('Freq')+xlab('')+theme_minimal() +
  labs(title="TCR specificity")+
  theme(title = element_text(face="bold"), axis.title.y = element_blank(), 
        axis.text.y = element_text(face="bold"),
        legend.position = "right", panel.grid.major.x =  element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"))
  
  
dev.off()


#calcualte the frequency to plot in boxplot 
pdf('/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/specificity_CD8_boxplot.pdf', width=11)
table(Pathology=tcells$Pathology_groups, cluster=tcells$cell_mk) %>% 
  data.frame()  %>% 
  filter(cluster %in% cd8.k) %>% 
  group_by(cluster) %>% 
  mutate(n_k=across(2, sum, na.rm = TRUE)) %>%
  mutate(prop=Freq/n_k$Freq) %>%
  ggplot(aes(Pathology, prop, color=cluster))+
    geom_boxplot(position=position_dodge(0.6), color="black")+
    geom_jitter(position=position_dodge(0.6), size=3.5)+
    theme_bw()+theme(axis.text.x = element_text(angle=60, hjust=1))+
    scale_color_manual(values=mycolors[cd8.k])

dev.off()
  


#tupes of tumor antigen
pdf('/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/specificity_tumor_cd8.pdf', width=12, height=7)
test=tcells@meta.data %>% filter(Pathology_groups=='tumor') 
table(Pathology=test$Pathology, cluster=test$cell_mk) %>% 
  data.frame()%>% filter(cluster %in% cd8.k) %>%
  ggplot(aes(cluster, Freq, fill=Pathology))+geom_col(position='fill')+ #position='fill' for prportion instead of number of cells
  coord_flip()+
  scale_fill_manual(name="sample",values = met.brewer('Tsimshian', 7))+
  ggtitle('')+ylab('Freq')+xlab('')+theme_minimal() +
  labs(title="TCR specificity")+
  theme(title = element_text(face="bold"), axis.title.y = element_blank(), 
        axis.text.y = element_text(face="bold"),
        legend.position = "right", panel.grid.major.x =  element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"))


dev.off()


pdf('/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/specificity_leiden.pdf', width=12, height=7)
table(Pathology=tcells$Pathology_groups, cluster=tcells$new_id) %>% 
  data.frame()%>% filter(cluster %in% cd8.id) %>%
  ggplot(aes(cluster, Freq, fill=Pathology))+geom_col(position='fill')+ #position='fill' for prportion instead of number of cells
  coord_flip()+
  scale_fill_manual(values = cols.path)+
  ggtitle('')+ylab('Freq')+xlab('')+theme_minimal() +
  labs(title="TCR specificity")+
  theme(title = element_text(face="bold"), axis.title.y = element_blank(), 
        axis.text.y = element_text(face="bold"),
        legend.position = "right", panel.grid.major.x =  element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"))


dev.off()


# leiden clusters
pdf('/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/TCR/specificity_leiden.pdf', width=12, height=7)
table(Pathology=tcells$Pathology_groups, cluster=tcells$new_id) %>% 
  data.frame()%>% filter(cluster %in% cd8.id) %>%
  ggplot(aes(cluster, Freq, fill=Pathology))+geom_col()+ #position='fill' for prportion instead of number of cells
  coord_flip()+
  scale_fill_manual(values = cols.path)+
  ggtitle('')+ylab('n_cells')+xlab('')+theme(
    panel.background = element_rect(fill = "white", colour = "grey50"),
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(0.9, 0.8)
  ) +  ylim(0, 1200) 
dev.off()

#calcualte the frequency to plot in boxplot 
library(ggrepel)
table(Pathology=tcells$Pathology_groups, cluster=tcells$leiden) %>% 
  data.frame()  %>% 
  filter(cluster %in% cd8.id) %>% 
  group_by(cluster) %>% 
  mutate(n_k=across(2, sum, na.rm = TRUE)) %>%
  mutate(prop=Freq/n_k$Freq) %>%
  ggplot(aes(Pathology, prop, color=cluster))+
  geom_boxplot(position=position_dodge(0.6), color="black")+
  geom_jitter(position=position_dodge(0.6))+
  theme_bw()+theme(axis.text.x = element_text(angle=60, hjust=1))+
  geom_label_repel(aes(label=cluster), position = position_jitter(seed = 1))
