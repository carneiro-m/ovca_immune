### NicheNet CD8 clusters - only dual expanded cells

#load libraries
library(Seurat)
library(dplyr)
library(nichenetr)
library(patchwork)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(viridisLite)
library(RColorBrewer)

###load files -----------------------------------------------------------------------------
#load tumor
tumor <- readRDS("~/analysis/tumor/data/tumor_final.rds") #updated 2022
#tumor$cell.id=factor(tumor$cell.id, levels=c( "CD4_1","CD4_2", "CD4_3", "CD4_Treg" , "CD4_Tfh", 
                                             # "CD8_ctx.1", "CD8_ctx.2", "CD8_ctx.3", "CD8_ctx.4", "CD8_exh" ,"T_prolif" ,"T_XCL" , "T_IFIT","T_mito",
                                            #  "NK",  "B" , "PC_1" ,"PC_2" ,"myl_1", "myl_2" , "pDC" ,  "unk_1", "epith"))
tumor
colnames(tumor@meta.data)
Idents(tumor) <- 'cell_mk'
DefaultAssay(tumor) <- "RNA"
tumor.norm <- tumor %>% NormalizeData() 
DefaultAssay(tumor) <- "SCT"

#load NicheNet files
ligand_target_matrix = readRDS("~/analysis/tumor/data/NicheNet_ligand_target_matrix.rds")
lr_network = readRDS("~/analysis/tumor/data/NicheNet_lr_network.rds")
weighted_networks = readRDS("~/analysis/tumor/data/NicheNet_weghted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

#load Tcells
tcells <- readRDS("data/T_cell_cluster.seu.obj_meta.rds")
tcells=subset(tcells, type=='D') #subset the dual expandeed clones
table(tcells$type)
tcells@meta.data %>% subset(type=='D') %>% nrow() #check


#CD8
cd8.k=levels(unique(Idents(tcells)))[grep('CD8',levels(Idents(tcells)) )]
names(cd8.k)=cd8.k
cd8=subset(tcells, idents=cd8.k)


# Nomalize tcells and CD8 -----
#tcells
DefaultAssay(tcells) <- "RNA"
tcells.norm <- tcells %>% NormalizeData() 
DefaultAssay(tcells) <- "SCT"

#cd8
DefaultAssay(cd8) <- "RNA"
cd8.norm <- cd8 %>% NormalizeData() 
DefaultAssay(tcells) <- "SCT"


### Prepare for NicheNet ------------------------------------------------------------------------------
## sender
sender_celltypes = unique(Idents(tumor.norm)) #all cells
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, tumor.norm, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
## Define gene set of interest
#markers=list()
#for(l in names(cd8.k)){
#  markers[[l]] = FindMarkers(subset(cd8.norm, idents=l), ident.1 = "tumor", min.pct = 0.1,
  #                           group.by = "tissue", test.use="MAST", assay="RNA") %>% 
 #                           rownames_to_column("gene")
  
#}

#saveRDS(markers, 'results/markers_dual_nicheNet.rds')
markers=readRDS('results/markers_dual_nicheNet.rds')

## run for each leiden cluster
dual_k=c('CD8_GZMB', 'CD8_GZMH')
names(dual_k)=dual_k
niche_r=list()
for (l in names(dual_k)){
  niche_r[[l]]=nichenet.my_function(l)
}

niche_top=list()
for (l in names(k)){
  niche_top[[l]]=nichenet.my_function_top25(l)
}



gzmk_val=nichenet.my_function_val('CD8_GZMK')
gzmh_val=nichenet.my_function_val('CD8_GZMH')
gzmb_val=nichenet.my_function_val('CD8_GZMB')

k=c('CD8_GZMK', 'CD8_GZMH', 'CD8_GZMB', 'CD8_CCL4', 'CD8_ZNF683', 'CD8_XCL1', 'CD8_CXCL13')
names(k)=k
niche_val=list()
for (l in names(k)){
  niche_val[[l]]=nichenet.my_function_val(receiver=l)
}

plot=list()
for (l in names(k)){
  plot[[l]]=niche_val[[l]]$ligand.person+niche_val[[l]]$dot.plot+niche_val[[l]]$ligand.network+
    niche_val[[l]]$ligand.receptor+niche_val[[l]]$ligand.receptor.strict+
    plot_layout(design=layout, guides ='collect')
}

plot$CD8_GZMK
plot$CD8_GZMH
plot$CD8_GZMB
plot$CD8_CCL4
plot$CD8_ZNF683

### save results in pdf ----

#check
niche_r$CD8_1$ligand.person+niche_r$CD8_1$dot.plot+plot_layout(widths = c(1, 5))
niche_r$CD8_5$ligand.receptor
niche_top$CD8_1$ligand.person+niche_top$CD8_1$dot.plot+plot_layout(widths = c(1, 5))


layout <- "
ABBBBBBBB
ABBBBBBBB
ABBBBBBBB
ABBBBBBBB
CCCCC####
CCCCC####
DDDDD####
DDDDD####
EEEEE####
EEEEE####
FFFFFFFFF
"


plot=list()
for (l in names(dual_k)){
plot[[l]]=niche_r[[l]]$ligand.person+niche_r[[l]]$dot.plot+niche_r[[l]]$ligand.network+
  niche_r[[l]]$ligand.receptor+niche_r[[l]]$ligand.receptor.strict+
  plot_layout(design=layout, guides ='collect')

}

plot_top=list()
for (l in names(k)){
  plot_top[[l]]=niche_top[[l]]$ligand.person+niche_top[[l]]$dot.plot+niche_top[[l]]$ligand.network+
    niche_top[[l]]$ligand.receptor+niche_top[[l]]$ligand.receptor.strict+
    plot_layout(design=layout, guides ='collect')
  
}

#CD8_1
pdf('/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/NicheNet/niche_dualCD8_GZMB.pdf', height=35, width=30)
plot$CD8_GZMB
dev.off()

pdf('/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/NicheNet/niche_dualtop25_CD8_GZMB.pdf', height=20, width=15)
plot_top$CD8_GZMB
dev.off()

#CD8_5
pdf('/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/NicheNet/niche_dualCD8_GZMH.pdf', height=40, width=30)
plot$CD8_GZMH
dev.off()

pdf('/Users/mayracarneiro/OneDrive/PhD/Data/scRNAseq/2021Jan_Figures/NicheNet/niche_dualtop25_CD8_GZMH.pdf', height=20, width=15)
plot_top$CD8_GZMH
dev.off()


### Function --------------------------------------------------------------------------
nichenet.my_function <- function(receiver=cluster_receiver){
  teste2 <- list()
  
## Receiver
receiver <- receiver
expressed_genes_receiver = get_expressed_genes(receiver, cd8.norm, pct = 0.10, assay_oi = "RNA")
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
DE_table_receiver =  markers[[l]]
geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

## Define a set of poteintial ligands
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

## Perform NicheNet
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities
best_upstream_ligands = ligand_activities %>% top_n(nrow(ligand_activities), pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

## Infer receptors and top-predicted target genes
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands",
                                                                    "Predicted target genes", color = "purple", x_axis_position = "top",
                                                                    legend_title = "Regulatory potential")  + theme(legend.text = element_text(size=7), legend.title = element_text(size=7),legend.position="bottom",
                                                                                                                    axis.text.x = element_text(face = "italic", size=3), axis.text.y = element_text(size=4)) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))

## Receptors of top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]
dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")+theme(legend.position = "bottom", legend.title = element_text(size=7), legend.text = element_text(size=7), axis.text.y = element_text(size=4))

## Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()
lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))
lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)
dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]
dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))
vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% 
  make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", 
                      x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")+
  theme(legend.position = "bottom",legend.text = element_text(size=7), legend.title = element_text(size=7))

## Summary
# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% subset(!rownames(vis_ligand_pearson) %in% c('HLA.DRA', 'HLA.DMA', 'HLA.A')) %>% 
  make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "#34495E",
                       legend_position = "bottom", x_axis_position = "top", legend_title = 
                        "Pearson correlation coefficient\ntarget gene prediction ability)") + 
                        theme(legend.text = element_text(size = 7, angle = 60, hjust = 1), 
                        legend.title = element_text(size = 7), axis.text.y = element_text(size=5))

# ligand expression Seurat dotplot
rotated_dotplot = DotPlot(subset(tumor, idents = sender_celltypes),assay="SCT", features = order_ligands, cols = "RdGy") + coord_flip() + 
  theme(legend.text = element_text(size = 7), legend.title = element_text(size = 7), legend.position = "bottom", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = 
          element_text( size = 5), axis.text.x = element_text(size = 9,  angle = 60,hjust = 0)) + ylab("Expression in Sender") + 
  xlab("") + scale_y_discrete(position = "right")



# ligand_target network 
return(test2 <- list(ligand.person=p_ligand_pearson, dot.plot=rotated_dotplot,ligand.network= p_ligand_target_network, ligand.receptor=p_ligand_receptor_network, 
                     ligand.receptor.strict=p_ligand_receptor_network_strict))
}
nichenet.my_function_top25 <- function(receiver=cluster_receiver){
  teste2 <- list()
  
  ## Receiver
  receiver <- receiver
  expressed_genes_receiver = get_expressed_genes(receiver, cd8.norm, pct = 0.10, assay_oi = "RNA")
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  DE_table_receiver =  markers[[l]]
  geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
  geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  
  ## Define a set of poteintial ligands
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
  
  ## Perform NicheNet
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
  ligand_activities
  best_upstream_ligands = ligand_activities %>% top_n(25, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
  
  ## Infer receptors and top-predicted target genes
  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands",
                                                                      "Predicted target genes", color = "purple", x_axis_position = "top",
                                                                      legend_title = "Regulatory potential")  + theme(legend.text = element_text(size=7), legend.title = element_text(size=7),legend.position="bottom",
                                                                                                                      axis.text.x = element_text(face = "italic", size=3), axis.text.y = element_text(size=4)) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
  
  ## Receptors of top-ranked ligands
  lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")+theme(legend.position = "bottom", legend.title = element_text(size=7), legend.text = element_text(size=7), axis.text.y = element_text(size=4))
  
  ## Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases
  lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
  receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()
  lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
  lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))
  lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
  lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)
  dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))
  vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
  p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% 
    make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", 
                        x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")+
    theme(legend.position = "bottom",legend.text = element_text(size=7), legend.title = element_text(size=7))
  
  ## Summary
  # ligand activity heatmap
  ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% subset(!rownames(vis_ligand_pearson) %in% c('HLA.DRA', 'HLA.DMA', 'HLA.A')) %>% 
    make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "#34495E",
                        legend_position = "bottom", x_axis_position = "top", legend_title = 
                          "Pearson correlation coefficient\ntarget gene prediction ability)") + 
    theme(legend.text = element_text(size = 7, angle = 60, hjust = 1), 
          legend.title = element_text(size = 7), axis.text.y = element_text(size=5))
  
  # ligand expression Seurat dotplot
  rotated_dotplot = DotPlot(subset(tumor, idents = sender_celltypes),assay="SCT", features = order_ligands, cols = "RdGy") + coord_flip() + 
    theme(legend.text = element_text(size = 7), legend.title = element_text(size = 7), legend.position = "bottom", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = 
            element_text( size = 5), axis.text.x = element_text(size = 9,  angle = 60,hjust = 0)) + ylab("Expression in Sender") + 
    xlab("") + scale_y_discrete(position = "right")
  
  
  
  # ligand_target network 
  return(test2 <- list(ligand.person=p_ligand_pearson, dot.plot=rotated_dotplot,ligand.network= p_ligand_target_network, ligand.receptor=p_ligand_receptor_network, 
                       ligand.receptor.strict=p_ligand_receptor_network_strict))
}
#validated interactions
nichenet.my_function_val <- function(receiver=cluster_receiver){
  teste2 <- list()
  
  ## Receiver
  receiver <- receiver
  expressed_genes_receiver = get_expressed_genes(receiver, cd8.norm, pct = 0.10, assay_oi = "RNA")
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  DE_table_receiver =  markers[[receiver]]
  geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
  geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  
  ## Define a set of poteintial ligands
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
  
  ## Perform NicheNet
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
  ligand_activities
  best_upstream_ligands = ligand_activities %>% top_n(nrow(ligand_activities), pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
  
  ## Infer receptors and top-predicted target genes
  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands",
                                                                      "Predicted target genes", color = "purple", x_axis_position = "top",
                                                                      legend_title = "Regulatory potential")  + theme(legend.text = element_text(size=7), legend.title = element_text(size=7),legend.position="bottom",
                                                                                                                      axis.text.x = element_text(face = "italic", size=3), axis.text.y = element_text(size=4)) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
  
  ## Receptors of top-ranked ligands
  lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")+theme(legend.position = "bottom", legend.title = element_text(size=7), legend.text = element_text(size=7), axis.text.y = element_text(size=4))
  
  ## Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases
  lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
  receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()
  lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
  lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))
  lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
  lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)
  dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))
  vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
  p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% 
    make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", 
                        x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")+
    theme(legend.position = "bottom",legend.text = element_text(size=7), legend.title = element_text(size=7))
  
  ## Summary
  # ligand activity heatmap
  ligand_pearson_matrix = ligand_activities %>% filter(test_ligand %in% colnames(vis_ligand_receptor_network_strict)) %>% 
    select(pearson, test_ligand) %>% column_to_rownames(var='test_ligand') %>% arrange(pearson) %>% as.matrix() 
  #rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  #colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  #vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = ligand_pearson_matrix  %>%
    make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "#34495E",
                        legend_position = "bottom", x_axis_position = "top", legend_title = 
                          "Pearson correlation coefficient\ntarget gene prediction ability)") + 
    theme(legend.text = element_text(size = 7, angle = 60, hjust = 1), 
          legend.title = element_text(size = 7), axis.text.y = element_text(size=5))
  
  # ligand expression Seurat dotplot
  targets=ligand_activities %>% filter(test_ligand %in% colnames(vis_ligand_receptor_network_strict)) %>% 
    arrange(pearson) %>% select(test_ligand)
  rotated_dotplot = DotPlot(subset(tumor, idents = sender_celltypes),assay="SCT", features = targets$test_ligand,
                              cols = "RdGy") + coord_flip() + 
    theme(legend.text = element_text(size = 7), legend.title = element_text(size = 7), legend.position = "bottom", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = 
            element_text( size = 5), axis.text.x = element_text(size = 9,  angle = 60,hjust = 0)) + ylab("Expression in Sender") + 
    xlab("") + scale_y_discrete(position = "right")
  
  
  
  # ligand_target network 
  return(test2 <- list(ligand.person=p_ligand_pearson, dot.plot=rotated_dotplot,ligand.network= p_ligand_target_network, ligand.receptor=p_ligand_receptor_network, 
                       ligand.receptor.strict=p_ligand_receptor_network_strict))
}
