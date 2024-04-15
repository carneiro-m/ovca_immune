#load libraries
library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(scran)
library(scater)
library(pheatmap)
library(dplyr)
library(tibble)
library(ggplot2)

#load data
tcells <- readRDS("~/analysis/data/T_cell_cluster.seu.obj_TCR.rds")
tumor <- readRDS('data/tumor_integrated_fout_regrssoutHS_2_TCR.rds')
blood <- readRDS("data/blood_integrated_seurat.rds")

ref <- HumanPrimaryCellAtlasData()

#Organize objects
DefaultAssay(tumor) <- "SCT"
tumor$cell.id=factor(tumor$cell.id, levels=c( "CD4_1","CD4_2", "CD4_3", "CD4_Treg" , "CD4_Tfh", 
 "CD8_ctx.1", "CD8_ctx.2", "CD8_ctx.3", "CD8_ctx.4", "CD8_exh" ,"T_prolif" ,"T_XCL" , "T_IFIT","T_mito",
   "NK",  "B" , "PC_1" ,"PC_2" ,"myl_1", "myl_2" , "pDC" ,  "unk_1", "epith"))
Idents(tumor) <- 'cell.id'
tumor$patient=sub("_.*", "", tumor$orig.ident)

tcells$Cell.ID <- factor(tcells$Cell.ID, levels=c("CD4_1","CD4_2","CD4_3","CD4_4", "CD4_5","Treg", "CD8_1","CD8_2","CD8_3", "CD8_4","CD8_5","CD8_6","CD8_7","CD8_8", "HS","myl", "prolif", "ukB", "uk1",  "uk2"))
Idents(tcells) <- 'Cell.ID'

blood$pop=factor(blood$pop, levels=c("CD4 T_1", "CD4 T_2", "CD4 T_3", "CD4 T_4", "CD4 T_5", 'Treg', 'CD8 T_1', 'CD8 T_2', 'CD8 T_3', 'CD8 T_4', 'CD8 T_5', 'CD8 T_6', 'CD8 T_7', 'CD8 T_8','CD8 T_prol', 'NK', 'NK_CD16', 'B_1', 'B_2', 'PCs', 'mono_CD14_1', 'mono_CD14_2', 'mono_CD14_3', 'mono_CD14_4', 'mono_CD14_CD16', 'mono_CD16', 'pDCs', "Megakaryocytes"))
Idents(blood) <- 'pop'

#as sce obj
sce <- as.SingleCellExperiment(tcells, assay = "SCT")
sce.tumor <- as.SingleCellExperiment(tumor, assay="SCT")
sce.blood <- as.SingleCellExperiment(blood, assay = "SCT")

colData(sce)
pred <- SingleR(sce, ref=ref, labels=ref$label.fine)
table(pred$labels)
plotScoreHeatmap(pred)


pred.tumor <- SingleR(sce.tumor, ref=ref, labels=ref$label.fine)
table(pred.tumor$labels)


pred.blood <- SingleR(sce.blood, ref=ref, labels=ref$label.fine)


#Correlation with clusters previously identified
tab <- table(Assigned=pred$pruned.labels, Cluster=sce$ident)
pdf('/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/Assign/Blueprint_tcells.pdf', width=10, height=8)
pheatmap(log2(tab+10), color=colorRampPalette(c("white", "#EE4AC9"))(101))
dev.off()

tab.tumor <- table(Assigned=pred.tumor$pruned.labels, Cluster=sce.tumor$cell.id)
pdf('/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/Assign/Blueprint_tumor.pdf', width=15, height=12)
pheatmap(log2(tab.tumor+10), cluster_cols=FALSE,color=colorRampPalette(c("white", "#EE4AC9"))(101))
dev.off()

tab.blood <- table(Assigned=pred.blood$pruned.labels, Cluster=sce.blood$pop)
pdf('/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/Assign/Blueprint_blood.pdf', width=15, height=12)
pheatmap(log2(tab.blood+10), color=colorRampPalette(c("white", "#EE4AC9"))(101))
dev.off()


