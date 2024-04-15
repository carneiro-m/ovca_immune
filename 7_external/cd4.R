
#load librarues
library(Seurat)
library(inborutils)
library(dplyr)
library(cowplot)
library(harmony)

##vdownload data -----
doi="10.5281/zenodo.5461803"
download_zenodo(doi = doi, path='./data', parallel = T)
download.file(url, destfile = 'meta.tar.gz')
untar(tarfile='meta.tar.gz')

options(timeout=100000)
url.data='https://zenodo.org/record/5461803/files/data.expression.tar.gz?download=1'
download.file(url.data, destfile = 'data/data.tar.gz')

#explore CD4 dataset
untar('data/data.tar.gz')
cd4 <- readRDS("~/analysis/external/data/expression/CD4/integration/CD4.thisStudy_10X.seu.rds")
cd4@meta.data %>% head()
table(cd4$batchV) 
table(cd4$dataset)
table(cd4$patient) 
#cd4 = cd4 %>% SCTransform()

### Dimension resuction -----
cd4.norm <- cd4 %>%
       Seurat::NormalizeData(verbose = FALSE) %>%
       FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
       ScaleData(verbose = FALSE) %>% 
       RunPCA(pc.genes = pbmc@var.genes, verbose = FALSE, npcs=100 )

ElbowPlot(cd4.norm, ndims=100)  
DimHeatmap(cd4.norm, dims = 1:10, cells = 500, balanced = TRUE)

cd4.norm <- cd4.norm %>%
            RunUMAP(dims=1:50)

DimPlot(cd4.norm)
DimPlot(cd4.norm, group.by='meta.cluster.coarse')
DimPlot(cd4.norm, split.by = 'meta.cluster')&NoLegend()

#Run harmony -----

cd4.norm= cd4.norm %>% RunHarmony("batchV", plot_convergence = TRUE)
cd4.norm <- cd4.norm %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
DimPlot(cd4.norm, reduction = "umap")
DimPlot(cd4.norm, group.by = 'meta.cluster', label=T)
DimPlot(cd4.norm, group.by = 'dataset', label=T)


#change resolution cluster 
cd4.norm=cd4.norm %>% FindClusters()

#markers in harmony clusters
DefaultAssay(cd4.norm)
markers=FindAllMarkers(cd4.norm, only.pos = T, test.use = 'MAST')
markers_fil=markers %>% group_by(cluster) %>% filter(p_val_adj<0.05)
write.csv(markers_fil, 'markers_cd4_harmony.csv')



levels(Idents(cd4.norm))
levels(cd4.norm)

new.id <- c('C0_FOXP3','C1_CCR7', 'C2_ANXA1', 'C3_FOS', 'C4_CCR6','C5_CH25H', 'C6_TNFRSF4','C7_TXNIP','C8_NKG7','C9_ALOX5AP', 
            'C10_STMN1', 'C11_GZMK','C12_HSPA1A', 'C13_ISG15', 'C14_CXCL13',  'C15_IGKV3')

names(new.id) <- levels(cd4.norm)
cd4.norm <- RenameIdents(cd4.norm, new.id)
cd4.norm$cell_mk <- Idents(cd4.norm) 


#DA

as_tibble(table(sample=cd4.norm$cancerType, cluster=cd4.norm$cell_mk)) %>% 
  add_count(sample, wt=n) %>% mutate(prop=n/nn*100)  %>%  #rename(total_sample=nn) %>%
  ggplot(aes(cluster, prop, color=sample))+geom_boxplot(position=position_dodge(0.6), color="dark gray")+
  geom_jitter(position=position_dodge(0.6), size=2.5)+theme_bw()+theme(axis.text.x = element_text(angle=60, hjust=1), 
                                                                       title = element_text(face="bold"))+ggtitle("Proportion of cell type per patient- tumor") 

as_tibble(table(sample=cd4.norm$cancerType, cluster=cd4.norm$cell_mk)) %>% 
  add_count(sample, wt=n) %>% mutate(prop=n/nn*100)  %>%
  ggplot( aes(sample, prop,color=sample ))+geom_boxplot( outlier.colour = NA) +  
  geom_jitter() + facet_wrap(~ cluster, scales = "free_y", ncol = 5) +theme( axis.text.x=element_text(angle=60, hjust=1),
                                                                             title = element_text(face="bold"))+ggtitle("Proportion of cell type per condition")+ylab("%")+scale_color_brewer(palette="Set3")


### Integration Seurat NOT DONE -----
#exclude patients with very dew cells: PACA.P20190909 and PACA.P20190515
cd4=subset(cd4, patient %in% c('PACA.P20190909', 'PACA.P20190515'),invert=T )
cd4.list <- SplitObject(cd4, split.by = "patient")

#normalization
#for (i in names(cd4.list)){
  cd4.list[[i]] <-NormalizeData(cd4.list[[i]],vars.to.regress = c("percent.mito"), verbose = F)
  cd4.list[[i]] <- RunPCA(cd4.list[[i]])
#}


#Integration Seurat - NOT performed now ------

options(future.globals.maxSize = 9000 * 1024^2)
#excludesample with few cells
cd4.list= cd4.list[-16]
features <- SelectIntegrationFeatures(object.list = cd4.list, nfeatures = 3000)
cd4.list <- PrepSCTIntegration(object.list = cd4.list, anchor.features = features, verbose = FALSE)

anchors <- FindIntegrationAnchors(object.list = cd4.list, normalization.method = "SCT", 
                                  anchor.features = features, k.filter = 80)


cd4.int <- IntegrateData(anchorset = anchors, normalization.method = "SCT", 
                         k.weight = 20,verbose = FALSE)

cd4.int <- RunPCA(cd4.int, verbose = FALSE, npcs = 100)
ElbowPlot(cd4.int, ndims=50) # Verify the number of PCs to be  included in UMAP
cd4.int <- RunUMAP(cd4.int, dims = 1:)
plots <- DimPlot(cd4.int, group.by = c("TumorSite", "orig.ident"))
plots & theme(legend.position = "top") & NoAxes() & guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3)))

cd4.int <- FindNeighbors(cd4.int, dims = 1:)
cd4.int <- FindClusters(cd4.int)
DimPlot(cd4.int)
saveRDS(icd4.int, "data/integrated_cd4_seu.rds")


