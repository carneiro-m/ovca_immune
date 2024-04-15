################################################################################
#################### Add VDJ results to Seurat obj  ############################

# After defining the final clustering of T cells from blood and tumor (integration with regression of heat-shock genes),
#we will here add the TCR-seq results to the gene expression data available in the Seurat obj

# The following steps will be performed:
# 1.rename the sample ID in the barcode to match TCR and GEX barcodes
# 2. Add TCR data to metadata in the Seurat obj
# 3. Save new seurat as: "/T_cell_cluster.seu.obj_TCR.rds"

#load libraries
library(dplyr)
library(tibble)
library(Seurat)

#load data
merged <- read.table("analysis/data/tcr.data.merge.tsv", sep="\t") #provided by Cheng
tcells <- readRDS("analysis/data/T_cell_cluster.seu.obj.rds")
tumor <- readRDS("analysis/data/tumor_integrated_f_r.rds")

################################################################################
###### 1.rename the sample ID in the barcode to match TCR and GEX barcodes

#remove "-1" in baarcode
merged$barcode <- gsub("-1", "", merged$barcode)
rownames(merged) <- merged$barcode

#change barcode names (patient label)
unique(merged$source)
merged$barcode <- gsub("PT9454_blood_tcr", "p09454_blood", merged$barcode)
merged$barcode <- gsub("PT9454_tumor_tcr", "p09454_tumor", merged$barcode)
merged$barcode <- gsub("PT10329_blood_tcr", "p10329_blood", merged$barcode)
merged$barcode <- gsub("PT9193_tumor_tcr", "p09293_tumor", merged$barcode)
merged$barcode <- gsub("PT9193_blood_tcr", "p09293_blood", merged$barcode)
merged$barcode <- gsub("PT9808_blood_tcr", "p09808_blood", merged$barcode)
merged$barcode <- gsub("PT9704_tumor_tcr", "p09704_tumor", merged$barcode)
merged$barcode <- gsub("PT10329_tumor_tcr", "p10329_tumor", merged$barcode)
merged$barcode <- gsub("PT9808_tumor_tcr", "p09808_tumor", merged$barcode)
merged <- merged %>% left_join(tcell_cluster.ID, by="barcode") %>% mutate(ID_type=paste(Tcell.ID, type, sep="_")) %>% column_to_rownames(var="barcode")

merged_cluster_type <- merged[c("type", "barcode")]
head(merged_cluster_type)

################################################################################
###### 2. Add TCR data to metadata in the Seurat obj
rownames(merged) <- merged$barcode
tcells <- AddMetaData(tcells, merged)
tail(tcells@meta.data)

################################################################################
###### 3. Save new seurat as: "/T_cell_cluster.seu.obj_TCR.rds"
saveRDS(tcells,"analysis/data/T_cell_cluster.seu.obj_TCR.rds")


### Some checks
as.data.frame(table(tissue=tcells$tissue, clon.type=tcells$ID_type)) %>% spread(clon.type, Freq)

nrow(filter(merged, type=="D"))
dual <- merged %>% filter(type=="D") %>% arrange(desc(frequency_ies))
dual %>% unique(dual$clotype)
table(merged$clonotype)
unique(dual$clotype)
merged$barcode <- gsub("PT9454")

