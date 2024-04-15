################################################################################ 
################## Integration of T cells from tumor and blood #################
## Aditional: regressing out cell cycle and heat shock proteins

# We tried to integrated only T cell clusters from blood and tumor, however we observed
# a cluster with cells expressing cell-cycle related genes and another expressing heat-shock/stree-related genes
# therefore, the enxt step was to regress out these variables in the scTransform step

# The following steps will be performed:
# 1. Subset T cells from each tissue
# 2. Get only count matrix from RNA assay
# 3. Split the seurat to have the count of each sample/patient separate
# 4. Score heat shock and cell cycle-related genes
# 5. Normalize: scTransform regressing out heat-shock related genes and cell-cycle related-genes
# 6. Integrate 
# 7. Dimension reduction and clustering

setwd("~/analysis")

#load libraries
library(Seurat)
library(dplyr)
library(ggplot2)

################################################################################
##### 1. Subset T cells from each tissue
integrated <- readRDS("2020Jan_all_tumor_samples/Petropoulos/tumor_integrated_f_r.rds")
unique(Idents(integrated))
tcells.t <- subset(integrated, idents=c("CD4", "CD8", "CD4_CD8", "mito", "stressed", "prolif.lympho", "Treg")) 

integrated.blood <- readRDS("2020Jan_all_blod_samples/blood_integrated_seurat.rds")
Idents(integrated.blood) <- "pop"
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
##### 4.Score heat shock and cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

heat.go <- read.table("tcell_integration/data/GO_0034605_cellular_response_heat_shock.txt", header = F, sep="\t")
head(heat.go)
heat.gene <- list(heat.go$V1)
  
for (i in names(tcells.list)){
  tcells.list[[i]] <- AddModuleScore(tcells.list[[i]], features = heat.gene, name="heat.shock", nbin = 10)
  tcells.list[[i]] <- CellCycleScoring(tcells.list[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
}


################################################################################
##### 5. Normalization: scTransform
for (i in seq_along(tcells.list)){
  tcells.list[[i]] <- SCTransform(tcells.list[[i]], vars.to.regress = c("percent.mt","G2M.Score", "S.Score","Phase","heat.shock1"), verbose = F)
}


################################################################################
##### 6. Integrate
options(future.globals.maxSize = 5000 * 1024^2)
tcells.features <- SelectIntegrationFeatures(object.list = tcells.list, nfeatures = 3000)
tcells.list <- PrepSCTIntegration(object.list = tcells.list, anchor.features = tcells.features, verbose = FALSE)

tcells.anchors <- FindIntegrationAnchors(object.list = tcells.list, normalization.method = "SCT", anchor.features = tcells.features, verbose = FALSE)
tcells.integrated <- IntegrateData(anchorset = tcells.anchors, normalization.method = "SCT", verbose = FALSE)

################################################################################
##### 7. Dimension reduction and clustering
tcells.integrated <- RunPCA(tcells.integrated, verbose = FALSE)
tcells.integrated <- RunUMAP(tcells.integrated, dims = 1:30)
plots <- DimPlot(tcells.integrated, group.by = c("tissue", "orig.ident"))
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3))) & ggtitle("Integration of T cells")


tcells.integrated <- FindNeighbors(tcells.integrated, dims = 1:30)
tcells.integrated <- FindClusters(tcells.integrated, resolution = 1.6)
DimPlot(tcells.integrated, label=T)+NoLegend()

