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



### Prepare for NicheNet ------------------------------------------------------------------------------
## sender
sender_celltypes = levels(Idents(tumor.norm)) #all cells
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, tumor.norm, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

## Define gene set of interest
markers=read.csv("/Users/mayracarneiro/OneDrive/PhD/manuscript/2021Jan_Figures/DEG/markers_optmized_tumor.csv") %>% 
        filter(p_val_adj<0.05) %>% 
        group_by(cluster) %>% 
        arrange(desc(avg_logFC), .by_group = T) %>%
        mutate(avg_log2FC=log2(10^avg_logFC)) %>% #transform from log to log2
        mutate(cluster= case_when(cluster == "CD4_1" ~ "CD4_LTB", 
                          cluster == "CD4_2" ~ "CD4_TMEM2",
                          cluster == "CD4_3" ~ "CD4_CD69",
                          cluster == "CD4_Tfh" ~ "CD4_CXCL13",
                          cluster == "CD4_Treg" ~ "CD4_Treg",
                          cluster == "CD8_exh" ~ "CD8_HAVCR2",
                          cluster == "CD8_ctx.1" ~ "CD8_GNLY",
                          cluster == "CD8_ctx.2" ~ "CD8_GZMB",
                          cluster == "CD8_ctx.3" ~ "CD8_GZMK",
                          cluster == "CD8_ctx.4" ~ "CD8_CCL4",
                          cluster == "T_prolif" ~ "T_prolif",
                          cluster == "T_IFIT" ~ "T_IFIT",
                          cluster == "T_mito" ~ "T_mito",
                          cluster == "unk_1" ~ "unk_1",
                          cluster == "PC_1" ~ "PC_1",
                          cluster == "PC_2" ~ "PC_2",
                          cluster == "B" ~ "B",
                          cluster == "T_XCL" ~ "T_XCL1",
                          cluster == "NK" ~ "NK",
                          cluster == "epith" ~ "epith",
                          cluster == "myl_2" ~ "myl_2",
                          cluster == "pDC" ~ "pDC",
                          cluster == "myl_1" ~ "myl_1"))#rename cluster with the new labels of cell clusters


## run for each leiden cluster
niche_r=list()
for (l in levels(Idents(tumor))[1:17]){
  niche_r[[l]]=nichenet.my_function(l)
}

niche_r_b=list()
for (l in levels(Idents(tumor))[19:23]){
  niche_r_b[[l]]=nichenet.my_function(l)
}

niche_top=list()
for (l in levels(Idents(tumor))[1:17]){
  niche_top[[l]]=nichenet.my_function_top25(l)
}

niche_top_b=list()
for (l in levels(Idents(tumor))[19:23]){
  niche_top_b[[l]]=nichenet.my_function_top25(l)
}

niche_top= c(niche_top, niche_top_b)

niche_val=list()
for (l in names(cd8.k)){
  niche_val[[l]]=nichenet.my_function_val(receiver=l)
}



### save results in pdf ----

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
for (l in levels(Idents(tumor))[1:17]){
plot[[l]]=niche_r[[l]]$ligand.person+niche_r[[l]]$dot.plot+niche_r[[l]]$ligand.network+
  niche_r[[l]]$ligand.receptor+niche_r[[l]]$ligand.receptor.strict+
  plot_layout(design=layout, guides ='collect')

}

plot_b=list()
for (l in names(niche_r_b)){
  plot_b[[l]]=niche_r_b[[l]]$ligand.person+niche_r_b[[l]]$dot.plot+niche_r_b[[l]]$ligand.network+
    niche_r_b[[l]]$ligand.receptor+niche_r_b[[l]]$ligand.receptor.strict+
    plot_layout(design=layout, guides ='collect')
  
}

plot_top=list()
for (l in levels(Idents(tumor))[1:17]){
  plot_top[[l]]=niche_top[[l]]$ligand.person+niche_top[[l]]$dot.plot+niche_top[[l]]$ligand.network+
    niche_top[[l]]$ligand.receptor+niche_top[[l]]$ligand.receptor.strict+
    plot_layout(design=layout, guides ='collect')
  
}

plot_top_b=list()
for (l in levels(Idents(tumor))[19:23]){
  plot_top_b[[l]]=niche_top[[l]]$ligand.person+niche_top[[l]]$dot.plot+niche_top[[l]]$ligand.network+
    niche_top[[l]]$ligand.receptor+niche_top[[l]]$ligand.receptor.strict+
    plot_layout(design=layout, guides ='collect')
  
}

#print the plots
for (l in names(plot)){
  pdf(paste0('/Users/mayracarneiro/OneDrive/PhD/manuscript/2021Jan_Figures/NicheNet/niche_all_CD45_tumor_', l, '.pdf'), height=20, width=15)
   plot=niche_r[[l]]$ligand.person+niche_r[[l]]$dot.plot+niche_r[[l]]$ligand.network+
   niche_r[[l]]$ligand.receptor+niche_r[[l]]$ligand.receptor.strict+
   plot_layout(design=layout, guides ='collect')
   print(plot)
  dev.off()
}

for (l in names(plot_b)){
  pdf(paste0('/Users/mayracarneiro/OneDrive/PhD/manuscript/2021Jan_Figures/NicheNet/niche_all_CD45_tumor_', l, '.pdf'), height=20, width=15)
  plot=niche_r_b[[l]]$ligand.person+niche_r_b[[l]]$dot.plot+niche_r_b[[l]]$ligand.network+
    niche_r_b[[l]]$ligand.receptor+niche_r_b[[l]]$ligand.receptor.strict+
    plot_layout(design=layout, guides ='collect')
  print(plot)
  dev.off()
}

for (l in names(plot_top)){
  pdf(paste0('/Users/mayracarneiro/OneDrive/PhD/manuscript/2021Jan_Figures/NicheNet/niche_top25_CD45_tumor_', l, '.pdf'), height=20, width=15)
  plot=niche_top[[l]]$ligand.person+niche_top[[l]]$dot.plot+niche_top[[l]]$ligand.network+
    niche_top[[l]]$ligand.receptor+niche_top[[l]]$ligand.receptor.strict+
    plot_layout(design=layout, guides ='collect')
  print(plot)
  dev.off()
}

for (l in names(plot_top_b)){
  pdf(paste0('/Users/mayracarneiro/OneDrive/PhD/manuscript/2021Jan_Figures/NicheNet/niche_top25_CD45_tumor_', l, '.pdf'), height=20, width=15)
  plot=niche_top[[l]]$ligand.person+niche_top[[l]]$dot.plot+niche_top[[l]]$ligand.network+
    niche_top[[l]]$ligand.receptor+niche_top[[l]]$ligand.receptor.strict+
    plot_layout(design=layout, guides ='collect')
  print(plot)
  dev.off()
}



#################################################################################
########################### Representation ######################################
#plots for revision  - supplementary figure
#top10 ligands
dotplots_list=list()
for (k in names(niche_top)){
  dotplots_list[[k]]=niche_top[[k]]$'dot.plot' + ggtitle(k)+theme_bw()+ 
    theme(legend.position = 'right', axis.text.x.top = element_text(angle=90, hjust=0),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.title = element_blank())&NoLegend() 
}

for (k in names(niche_top)[-c(1, 2)]){
  dotplots_list[[k]]=dotplots_list[[k]] + theme(axis.text.x.top = element_blank(), 
                                                axis.title.x = element_blank())
}

for (k in names(niche_top)[-c(21,22)]){
  dotplots_list[[k]]=dotplots_list[[k]] + theme(axis.text = element_text(angle=90, hjust=1), 
                                                axis.title.x = element_blank())
}

pdf("/Users/mayracarneiro/OneDrive/PhD/manuscript/2021Jan_Figures/NicheNet/niche_all_CD45_tumor_top10.pdf", width = 9, height = 20)
wrap_plots(dotplots_list, ncol= 2)
dev.off()

pdf("/Users/mayracarneiro/OneDrive/PhD/manuscript/2021Jan_Figures/NicheNet/niche_all_CD45_tumor_top10_v2.pdf", width = 5, height = 40)
wrap_plots(dotplots_list, ncol= 1)
dev.off()

#################################################################################
### Function --------------------------------------------------------------------------
nichenet.my_function <- function(receiver=cluster_receiver){
  teste2 <- list()
  
## Receiver
receiver <- receiver
expressed_genes_receiver = get_expressed_genes(receiver, tumor.norm, pct = 0.10, assay_oi = "RNA")
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## Define gene set of interest  
DE_table_receiver =  markers %>% filter(cluster==receiver)  
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
  expressed_genes_receiver = get_expressed_genes(receiver, tumor.norm, pct = 0.10, assay_oi = "RNA")
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  DE_table_receiver =  markers %>% filter(cluster==receiver)
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
  best_upstream_ligands = ligand_activities %>% top_n(10, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
  
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
                       ligand.receptor.strict=p_ligand_receptor_network_strict
                       ))
}
