### SingleR to assign celll types in T cells

#load libraries
setwd("analysis/tcells")
library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(scran)
library(scater)
library(pheatmap)
library(dplyr)
library(tibble)

#load data
tcells <- readRDS("~/analysis/TCR/Data/T_cell_cluster.seu.obj_TCR.rds")
tumor <- readRDS("~/analysis/2020Jan_all_tumor_samples/Petropoulos/tumor_integrated_f_r.rds")
blood <- readRDS("~/analysis/2020Jan_all_blod_samples/blood_integrated_seurat.rds")

monaco <- readRDS("data/Monaco.rds")
unique(monaco$label.main)


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
