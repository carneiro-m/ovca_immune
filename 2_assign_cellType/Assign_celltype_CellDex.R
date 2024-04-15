#load libraries
setwd("analysis")
library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(scran)
library(scater)
library(pheatmap)
library(dplyr)
library(tibble)

#load data
tcells <- readRDS("data/T_cell_cluster.seu.obj_TCR.rds")
tumor <- readRDS('data/tumor_integrated_fout_regrssoutHS_2_TCR.rds')
blood <- readRDS("data/blood_integrated_seurat.rds")

monaco <- readRDS("data/Monaco.rds")
unique(monaco$label.main)


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
pred <- SingleR(sce, ref=monaco, labels=monaco$label.fine, method="single")
plotScoreHeatmap(pred)

pred.tumor <- SingleR(sce.tumor, ref=monaco, labels=monaco$label.fine, method="single")
plotScoreHeatmap(pred.tumor)
saveRDS(pred.tumor, "data/labels_SingleR_tcells.rds")

pred.blood <- SingleR(sce.blood, ref=monaco, labels=monaco$label.fine, method="single")
plotScoreHeatmap(pred.blood)

#Correlation with clusters previously identified
tab <- table(Assigned=pred$pruned.labels, Cluster=sce$ident)
pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))+ggtitle("tcells")

tab.tumor <- table(Assigned=pred.tumor$pruned.labels, Cluster=sce.tumor$pop)
pheatmap(log2(tab.tumor+10), color=colorRampPalette(c("white", "blue"))(101))+ggtitle("tumor")

tab.blood <- table(Assigned=pred.blood$pruned.labels, Cluster=sce.blood$pop)
pheatmap(log2(tab.blood+10), color=colorRampPalette(c("white", "blue"))(101))

Idents(tumor) <- "pop"
tumor$barcode <- rownames(tumor@meta.data)
tfh.tumor <- rownames(subset(pred.tumor, pruned.labels=="Follicular helper T cells"))
cd4_4 <- WhichCells(subset(tumor, idents="CD4_4"))
tfh <-dplyr::intersect(tfh.tumor, cd4_4)
df=data.frame(barcode=tfh, new="Tfh")
test <- tumor@meta.data[, c("barcode", "pop")] %>% left_join(df, by="barcode") %>% mutate(name=coalesce(new, pop)) 
tumor$name <- test$name
table(tumor$pop, tumor$name)

(DimPlot(tumor, cells.highlight = tfh.tumor, sizes.highlight = 0.0001)+ggtitle("Tfh identified by SingleR in tumor"))|
(DimPlot(tcells, cells.highlight = tfh.tumor, sizes.highlight = 0.0001)+ggtitle("Tfh identified by SingleR in tcells"))



markers.th.cd4.4 <- FindMarkers(tumor.norm, ident.1 = "Tfh",ident.2="CD4_4",test.use="MAST")
markers.th.cd4.4 <- markers.th.cd4.4 %>% filter(p_val_adj<0.05) %>% rownames_to_column(var = "gene") %>% arrange(desc(avg_logFC))
EnhancedVolcano(markers.th.cd4.4, x=avg_logFC, y=p_val_adj, label=gene)








test.2 <- tumor@meta.data[, c("barcode", "groups")] %>% left_join(df, by="barcode") %>% mutate(name=coalesce(new, groups)) 
tumor$name2 <- test.2$name
table(tumor$pop, tumor$name2)
Idents(tumor) <- "name2"
markers.th.cd4 <- FindMarkers(tumor, ident.1 = "Tfh", ident.2 = "CD4",test.use="MAST", only.pos = T)
markers.th.cd4 <- markers.th.cd4 %>% filter(p_val_adj<0.05) %>% rownames_to_column(var = "gene")
markers.th.cd4$gene
write.csv(markers.th.cd4, "~/analysis/tcells/markers_tfh_vs_CD4.csv")

map=mapIds(org.Hs.eg.db, markers.th$gene, 'ENTREZID', 'SYMBOL')
columns(org.Hs.eg.db)
