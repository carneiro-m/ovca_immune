################################################################################
####### TCR - calculate frequency and classify by site of expansion ############

# In this study, the T cell clones were classified according to the site of expansion (blood, tumor, or both)
#and the frequency of a clone was calcualted within a tissue (blood or tumor).

# Because the initial VDJ analysis was performed in cellranger VDJ output file (file "TCR_metadata_VDJ_results.R"), 
#which didn't consider the cell filtering in the QC step, I have recalculated the frequency as well as
#the expansion pattern (dual expanded - blood and tumor, singleton, multiplet) 



#load libraries 
library(dplyr)
library(tibble)
library(Seurat)
library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(patchwork)
library(colorspace)
library(scales)
library (ComplexHeatmap)


##load data
tcells <- readRDS("data/T_cell_cluster.seu.obj_meta.rds")


### Obtain TCR and cluster data from seurat obj -----
data <- tcells@meta.data  %>% rownames_to_column(var="barcode")
#clen up data
data <- data %>% select(c( barcode,tissue,patient,cell_mk, cdr3s_aa, cdr3s_nt, clotype))  %>% 
  drop_na() %>% #filter(Cell.ID %in% unique(data$Cell.ID)[grep("CD", unique(data$Cell.ID))]) %>% #select inecessary informations, obtain barcode ID and remove cells that we do not have TCR-seq information; Filter only clusters in CD4 and CD8
  #separate(orig.ident, sep="_", into=c("patient", "tissue")) %>%
  mutate_if(is.character, str_replace_all, pattern = "p09193", replacement = "p09293") %>%#since the label of patient 9293 is wrong in one of the samples, I will fix it 
  filter(patient != 'p09704') #remove patient who data is only for tumor tissue 

### Recalculate the type classification (singleton and dual expanded) ------------------------------------------------------------------------------------------------------------------------
#I will calculate it again because some of cells that allowed some clone be classified as dual were filtered iut in Seurat piperline.

data_type <- data %>% select(c( barcode, patient,tissue, clotype)) %>%
  add_count(patient, clotype, tissue, name="n") %>% #count the frequency of clone in each tissue of a sample/patient
  select(-barcode) %>% distinct(clotype,tissue,.keep_all = T) %>%  # since we already counted, I will remove the barcode collumn and select unique row for a clonotype within tissue
  spread(tissue,n ) %>% replace(is.na(.), 0)%>% 
  mutate(freq=coalesce(blood, tumor)) %>% mutate(type=(ifelse(blood==1 & tumor==0,"b",ifelse(blood > 1 & tumor==0,"B", ifelse(blood >= 1 & tumor >=1,"D",
                                                                                                                              ifelse(blood==0 & tumor==1,"t",ifelse(blood==0 & tumor>1,"T","NS"))))))) %>% #classify by tissue expansion
  mutate(tumor=na_if(tumor, 0), blood=na_if(blood, 0)) %>% select(clotype, type, patient)


data=data %>% left_join(data_type, by=c('clotype', 'patient'))


