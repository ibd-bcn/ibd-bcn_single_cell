# Title: Cross-tissue single-cell landscape of
#        human monocytes and macrophages in health and disease

#
# Libraries and extras
#
library(Seurat)
cols_myeloids <- readRDS("~/000_GitHub/ibd-bcn_single_cell/source/colors/cols_myeloids.RDS")

#
# Our data
#

myeloids <- readRDS('~/000_GitHub/ibd-bcn_single_cell/Analysis of our data/02_Samples_Together/SUBSETS/ON_THEIR_OWN/myeloids_annotated.RDS')

data.query <- myeloids[,myeloids$annotation_intermediate %in% c("IDA macrophage",
                                                                    "M0"  ,  "M1" ,
                                                                    "Inflammatory monocytes", "M2",
                                                                    "Cycling myeloid")]
DimPlot(data.query)
data.query <- FindVariableFeatures(data.query, nfeatures = 3000)
Idents(data.query) <- data.query$annotation_intermediate

#
# MoMac dataset 
#

data.integrated <- readRDS('~/data_Albas/All_Together_31082021/08_VERSE/2021_MoMac_VERSE.rds')
data.integrated <- SetIdent(data.integrated, value = 'Clusters')
data.integrated <- FindVariableFeatures(data.integrated, nfeatures = 3000)

ElbowPlot(data.integrated, ndims = 100)
length(intersect(data.integrated@assays$RNA@var.features, data.query@assays$RNA@var.features))

data.anchors <- FindTransferAnchors(reference = data.integrated, query = data.query, 
                                    dims = 1:30, reference.reduction = "pca")

srt <- MapQuery( anchorset = data.anchors, query = data.query, reference = data.integrated,
                 reference.reduction = "pca",  reduction.model = "umap" )

a <- DimPlot(srt, reduction = 'ref.umap', group.by = 'annotation_refined') + 
  theme_void() + 
  scale_color_manual(values = cols_myeloids[as.character(unique(srt$annotation_refined))[order(as.character(unique(srt$annotation_refined)))]]) +
  theme(title = element_blank())

b <- DimPlot(data.integrated, label =T) + theme_void() + 
  theme(legend.position = 'none')

b+a

for(i in unique(srt$annotation_refined)){
  a <- DimPlot(srt, reduction = 'ref.umap', group.by = 'annotation_refined',
               cells.highlight =  list(cluster = colnames(srt[,srt$annotation_refined == i]))) + 
    theme_void() + labs(title = i) + NoLegend()
  k <- b+a
  ggsave(filename = gsub('/','__',gsub(' ', '_',
                                       paste('plot_MoMacProjection', i, '.jpeg', sep='_'))),
         plot = k, dpi = 300, device = 'jpeg')
}

