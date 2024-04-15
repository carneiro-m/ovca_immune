################################################################################ 
################## Integration of T cells from tumor and blood #################

# We start the analysis of all cells from each tissue (blood and tumor) separately.
# We first tried to integrate all cells from blood and tumor, however, this analysis 
#displayed a different profile in gene signature of T cells. This could be more clear 
#in the exhausted T cluster. This cluster, very well characterized in the tumor, was 
#not evident in the integrated object.

# Considering that, we decided to integrated only T cell clusters from blood and tumor
# The following steps will be performed:
# 1. Subset T cells from each tissue
# 2. Get only count matrix from RNA assay
# 3. Split the seurat to have the count of each sample/patient separate
# 4. Normalize: scTransform
# 5. Integrate 
# 6. Dimension reduction and clustering
# 7. Perform DE analysis to identify markers for each cluster
# 8. Visualization to explore the results


#load libraries
setwd("~/analysis")
library(Seurat)
library(dplyr)
library(ggplot2)

################################################################################
##### 1. Subset T cells from each tissue
integrated <- readRDS("2020Jan_all_tumor_samples/Petropoulos/tumor_integrated_f_r.rds")
unique(Idents(integrated))
tcells.t <- subset(integrated, idents=c("CD4", "CD8", "CD4_CD8", "mito", "stressed", "prolif.lympho", "Treg")) 

integrated.blood <- readRDS("2020Jan_all_blod_samples/blood_integrated_seurat.rds")
unique(Idents(integrated.blood))
tcells.b <- subset(integrated.blood, idents = c("CD8 T_1","CD8 T_2","CD8 T_3", "CD8 T_4",
            "CD8 T_5", "CD8 T_6", "CD8 T_7", "CD8 T_8", "CD8 T_prol", "CD4 T_1", "CD4 T_2",
            "CD4 T_3", "CD4 T_4", "CD4 T_5", "Treg"))

################################################################################
##### 2. Get only RNA assay
DefaultAssay(tcells.t) <- "RNA"
tcells.t <- DietSeurat(tcells.t, assays = "RNA")

DefaultAssay(tcells.b) <- "RNA"
tcells.b <- DietSeurat(tcells.b, assays = "RNA")  

################################################################################
##### 3. Split seurat obj to have the count of each sample/patient separate
tcells.t.list <- SplitObject(tcells.t, split.by = "orig.ident")
tcells.b.list <- SplitObject(tcells.b, split.by = "orig.ident")
tcells.list <- c(tcells.t.list, tcells.b.list)
tcells.list  

################################################################################
##### 4. Normalization: scTransform
for (i in seq_along(tcells.list)){
  tcells.list[[i]] <- SCTransform(tcells.list[[i]], vars.to.regress = "percent.mt", verbose = F)
}

################################################################################
##### 5. Integrate
options(future.globals.maxSize = 5000 * 1024^2)
tcells.features <- SelectIntegrationFeatures(object.list = tcells.list, nfeatures = 3000)
tcells.list <- PrepSCTIntegration(object.list = tcells.list, anchor.features = tcells.features, verbose = FALSE)

tcells.anchors <- FindIntegrationAnchors(object.list = tcells.list, normalization.method = "SCT", anchor.features = tcells.features, verbose = FALSE)
tcells.integrated <- IntegrateData(anchorset = tcells.anchors, normalization.method = "SCT", verbose = FALSE)


################################################################################
##### 6. Dimension reduction and clustering
tcells.integrated <- RunPCA(tcells.integrated, verbose = FALSE)
tcells.integrated <- RunUMAP(tcells.integrated, dims = 1:30)
plots <- DimPlot(tcells.integrated, group.by = c("tissue", "orig.ident"))
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3))) & ggtitle("Integration of T cells")


tcells.integrated <- FindNeighbors(tcells.integrated, dims = 1:30)
tcells.integrated <- FindClusters(tcells.integrated)
DimPlot(tcells.integrated, label=T)+NoLegend()

#Increasing resolution for cluster 
DefaultAssay(tcells.integrated) <- "integrated"
tcells.integrated.hires <- FindClusters(tcells.integrated, resolution = 1.6)
DimPlot(tcells.integrated.hires, label=T)+NoLegend()

################################################################################
##### 7. Perform DEG analysis to identify markers for each cluster

DefaultAssay(tcells.integrated) <- "RNA"
tcells.integrated.norm <- NormalizeData(tcells.integrated)
markers.tcells.integrated <- FindAllMarkers(tcells.integrated.norm, only.pos = T, test.use = "MAST")
head(markers.tcells.integrated)
markers.tcells.integrated %>% count(cluster)
markers.tcells.integrated.fil <- filter(markers.tcells.integrated, p_val_adj<0.05)
markers.tcells.integrated.fil %>% count(cluster)
saveRDS(markers.tcells.integrated.fil, "~/analysis/tcells/data/markers_tcells_integrated_bloodtumor.rds")
write.csv(markers.tcells.integrated.fil, "~/analysis/markers_tcells_integrated_bloodtumor.csv")
markers.tcells.integrated.fil <- readRDS("~/analysis/markers_tcells_integrated_bloodtumor.rds")

cnserved.markers.tcells <-FindConservedMarkers(tcells.integrated.norm,ident.1 = "ribo", grouping.var = "tissue")
conserved.mito.markerr <-FindConservedMarkers(tcells.integrated.norm,ident.1 = "mito", grouping.var = "tissue")

DefaultAssay(tcells.integrated.hires) <- "RNA"
tcells.integrated.hires.norm <- NormalizeData(tcells.integrated.hires)
markers.tcells.integrated.hires <- FindAllMarkers(tcells.integrated.hires.norm, only.pos = T, test.use = "MAST")
head(markers.tcells.integrated.hires)
markers.tcells.integrated.hires %>% count(cluster)
markers.tcells.integrated.hires.fil <- filter(markers.tcells.integrated.hires, p_val_adj<0.05)
markers.tcells.integrated.hires.fil %>% count(cluster)
saveRDS(markers.tcells.integrated.hires.fil, "~/analysis/markers_tcells_integrated_hires_bloodtumor.rds")
markers.tcells.integrated.hires.fil <- readRDS("~/analysis/markers_tcells_integrated_bloodtumor.rds")

################################################################################
##### 8. Visualization to explore the results
DefaultAssay(tcells.integrated.hires) <- "SCT"
FeaturePlot(tcells.integrated.hires, features = c("CD3G", "CD4", "CD8A", "FOXP3"), order=T, ncol=4, label=T)/  
VlnPlot(tcells.integrated.hires, features = c("CD3G", "CD4", "CD8A", "FOXP3"), pt.size = 0, ncol = 4)+RotatedAxis()

FeaturePlot(tcells.integrated.hires, label=T, features = c("CD79A", "CD19", "CD8A", "MS4A1"), order=T, ncol=4)/  
  VlnPlot(tcells.integrated.hires, features = c("CD79A", "CD19", "CD8A", "MS4A1"), pt.size = 0, ncol = 4)+RotatedAxis()

FeaturePlot(tcells.integrated, features = c("CD3G", "FCGR3A", "NCAM1", "CD79A"), order=T, ncol=4)/  
  VlnPlot(tcells.integrated, features = c("CD3G", "FCGR3A", "NCAM1", "CD79A"), pt.size = 0, ncol = 4)

FeaturePlot(tcells.integrated, features = c("PDCD1", "HAVCR2", "ENTPD1", "CXCL13"), order=T, ncol=4, label=T)/  
  VlnPlot(tcells.integrated, features = c("PDCD1", "HAVCR2", "ENTPD1", "CXCL13"), pt.size = 0, ncol = 4)

FeaturePlot(tcells.integrated, features = c("PDCD1", "TOX", "TCF7", "SELL"), order=T, ncol=4)/  
  VlnPlot(tcells.integrated, features = c("PDCD1", "TOX", "TCF7", "SELL"), pt.size = 0, ncol = 4)


FeaturePlot(tcells.integrated, features = c("ENTPD1", "CXCL13"),
            order=T, ncol=4, label=T, split.by = "tissue")

##Rename clusters
new.cluster.ids <- c("CD4_1", "CD8_1", "CD4.CD8_1", "CD4_2", "Treg", "CD8_3", "CD8_4", "CD4_3", 
"CD4_4", "mito",  "CD8_5", "ribo",  "CD4_5" , "B","prolif","HS_Tcell", "CD8_6","CD4.CD8_2",
"CD4.CD8_3","PCs", "CD8_7")
names(new.cluster.ids) <- levels(tcells.integrated)
tcells.integrated <- RenameIdents(tcells.integrated, new.cluster.ids)
tcells.integrated$pop <- Idents(tcells.integrated)
DimPlot(tcells.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(tcells.integrated, "tcells_integrated_blood_tumor.rds")
saveRDS(tcells.integrated, "tcells_integrated_blood_tumor.rds.gz", compress="gzip")

## DEG for CD8 and CD4 separately
#DEG CD8
t.cd8 <- subset(tcells.integrated, idents=c("1", "2", "5", "6", "10","14","16", "17", "18", "20"))
DefaultAssay(tcd8) <- "RNA"
tcd8.norm <- NormalizeData(t.cd8)
markers.t.cd8 <- FindAllMarkers(tcd8.norm, only.pos = T, test.use = "MAST")
head(markers.t.cd8)
markers.t.cd8 %>% count(cluster)
markers.t.cd8.fil <- filter(markers.t.cd8, p_val_adj<0.05)
markers.t.cd8.fil %>% count(cluster)
saveRDS(markers.t.cd8.fil, "~/analysis/markers_tcd8_integrated_bloodtumor.rds")

#DEG CD8_5 tumor versus blood
new.cluster.ids <- c("CD4_1", "CD8_1", "CD4.CD8_1", "CD4_2", "Treg", "CD8_3", "CD8_4", "CD4_3", 
                     "CD4_4", "mito",  "CD8_5", "ribo",  "CD4_5" , "B","prolif","HS_Tcell", "CD8_6","CD4.CD8_2",
                     "CD4.CD8_3","PCs", "CD8_7")
names(new.cluster.ids) <- levels(tcells.integrated.norm)
tcells.integrated.norm <- RenameIdents(tcells.integrated.norm, new.cluster.ids)
tcells.integrated.norm$pop <- Idents(tcells.integrated.norm)
cd8.5 <- FindMarkers(tcells.integrated.norm, ident.1 = "tumor", ident.2 = "blood", group.by = "tissue",
                     subset.ident = "CD8_5", test.use = "MAST", only.pos = T)
cd8.5 <- filter(cd8.5, p_val_adj<0.05)

cd8.7 <- FindMarkers(tcells.integrated.norm, ident.1 = "tumor", ident.2 = "blood", group.by = "tissue",
                     subset.ident = "CD8_7", test.use = "MAST", only.pos = T)
cd8.7 <- filter(cd8.7, p_val_adj<0.05)

cd8.3 <- FindMarkers(tcells.integrated.norm, ident.1 = "tumor", ident.2 = "blood", group.by = "tissue",
                     subset.ident = "CD8_3", test.use = "MAST", only.pos = T)
cd8.3 <- filter(cd8.3, p_val_adj<0.05)

##Distribution of clusters by tissue
table.t=table(tcells.integrated$tissue, Idents(tcells.integrated))
prop.cluster.t <-  as.data.frame(table.t)
table.t

ggplot(prop.cluster.t, aes(Var2, Freq, fill=Var1))+geom_col(position="fill", colour="#F4F6F6")+theme(axis.text.x = element_text(hjust=1, angle=60), axis.title.x = element_blank())+labs(fill="tissue")+ylab("Frequency")+ggtitle("")



