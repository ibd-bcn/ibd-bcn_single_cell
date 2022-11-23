# Title	          Single-cell analysis of Crohn’s disease lesions identifies a pathogenic
#                 cellular module associated with resistance to anti-TNF therapy
# Overall design	Single cell seq analysis of paired resection ileal samples from inflamed and uninflamed area
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134809
# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA556461&o=acc_s%3Aa

#
# Libraries: ------------------------------------------------------------------
#
library(Seurat)
library(plyr)
library(ggplot2)
library(DropletUtils)
library(beepr)
library(celda)
library(SingleCellExperiment)
library(scater)
library(scran)
library(scDblFinder)
require(formulaic)
source('source/functions_scrnaseq.R')

#
# data ------------------------------------------------------------------------
#
metadata <- read.table('~/data_Albas/GSE134809_RAW/SraRunTable.txt',
                       sep = ',', header=T)

metadata <- metadata[metadata$Tissue == 'ileal',]
setwd('~/data_Albas/GSE134809_RAW/')

list_data <- list()

for(i in unique(metadata$Sample.Name)){
  files <- grep(i,list.files('~/data_Albas/GSE134809_RAW/'), value = T)
  print(files)
  from1 = files
  to1 = c('barcodes.tsv', 'genes.tsv', 'matrix.mtx')
  file.rename(from1, to1)
  
  
  sce2 <- Read10X('~/data_Albas/GSE134809_RAW/')
  
  m4 <- data.frame('sample' = rep(i, ncol(sce2)),
                   'Health' = rep(metadata$status[metadata$Sample.Name == i][1], ncol(sce2))
  )
  
  rownames(m4) <- colnames(sce2)
  
  data <- CreateSeuratObject(
    sce2, min.features = 100,
    project = i,
    assay = "RNA",
    meta.data = m4
  )
  data[["percent.mt"]] <- PercentageFeatureSet(object = data, pattern = "^MT-")
  
  
  file.rename(from = to1, to = from1)
  list_data[[i]] <- data
}

seudata <- list_data[[1]]
for(i in 2:length(unique(metadata$Sample.Name))){
  seudata <- merge(seudata, list_data[[i]])
}

VlnPlot(seudata, features = c('percent.mt', 'nFeature_RNA'), group.by = 'Health', pt.size = 0)
VlnPlot(seudata[, seudata$percent.mt < 25 & seudata$nFeature_RNA > 100],
        features = c('percent.mt', 'nFeature_RNA'), group.by = 'Health', pt.size = 0)


seudata_f <- seudata[, seudata$percent.mt < 25 & seudata$nFeature_RNA > 100]

seudata_f <- seurat_to_pca(seudata_f)

p <- ElbowPlot(seudata_f, ndims = 100) +
  geom_vline(xintercept = 25, colour="#BB0000")+
  annotate(geom="text", x=20, y=4, label= paste("sdev > 1.6\nPCs =", PCS2),
           color="#BB0000")
ggsave(p, file = 'elbowplot_Cho.jpeg')

seudata_f <- FindNeighbors(seudata_f, dims = 1:25)
seudata_f <- RunUMAP(seudata_f, dims = 1:25)

DimPlot(seudata_f, group.by = 'sample') + labs(title = 'Paper J.Cho')
DimPlot(seudata_f, group.by = 'Health') + labs(title = 'Paper J.Cho')

FeaturePlot(seudata_f, features = c('EPCAM', 'CD3E', 'C1QB', 'DERL3'), order =T, cols = c('lightgray', 'red') )
FeaturePlot(seudata_f, features = c('C1QA', 'C1QB', 'TPSB2', 'AIF1'), order =T, cols = c('lightgray', 'red'))
seudata_f <- FindClusters(seudata_f, resolution = 0.1)
DimPlot(seudata_f, label = T)

myeloids_cho <- seudata_f[,seudata_f$seurat_clusters %in% c('3', '8')]

dir.create('myeloids_cho_with_masts')

FeaturePlot(myeloids_cho, features = c('CD4', 'MS4A1', 'DERL3', 'C1QB'), order =T, cols = c('lightgray', 'red'))

myeloids_cho <- seurat_to_pca(myeloids_cho)
setwd('~/data_Albas/GSE134809_RAW/')

PCS <- select_pcs(myeloids_cho, 2)
PCS2 <- select_pcs(myeloids_cho, 1.6) # 22
p <- ElbowPlot(myeloids_cho, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0(' myeloids_cho_with_masts'))
ggsave(p, file = 'elbowplot_myeloids_cho_with_masts.jpeg')

myeloids_cho <- FindNeighbors(myeloids_cho, dims = 1:PCS2)
myeloids_cho <- RunUMAP(myeloids_cho, dims = 1:PCS2)

DimPlot(myeloids_cho, group.by = 'sample') + labs(title = 'Paper J.Cho')
DimPlot(myeloids_cho, group.by = 'Health') + labs(title = 'Paper J.Cho')

FeaturePlot(myeloids_cho, features = c('CD3E', 'MS4A1', 'TPSB2', 'C1QB'), order =T, cols = c('lightgray', 'red'))
myeloids_cho <- FindClusters(myeloids_cho, resolution = 0.1)
DimPlot(myeloids_cho, label=T)

myeloids_cho <- resolutions(myeloids_cho, workingdir = 'myeloids_cho_with_masts', title = 'myeloids_cho')

# myeloids with masts filtering --------------------------------------------------------------------

Idents(myeloids_cho) <- myeloids_cho$RNA_snn_res.1.5
myeloids_cho <- subset(myeloids_cho, idents = c(0:4, 6:19, 22:24))

myeloids_cho <- seurat_to_pca(myeloids_cho)
setwd('~/data_Albas/GSE134809_RAW/')

PCS <- select_pcs(seudata_f, 2)
PCS2 <- select_pcs(seudata_f, 1.6)
p <- ElbowPlot(seudata_f, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0(' myeloids_cho_with_masts'))
ggsave(p, file = 'elbowplot_myeloids_cho_with_masts_filtering.jpeg')

myeloids_cho <- FindNeighbors(myeloids_cho, dims = 1:PCS2)
myeloids_cho <- RunUMAP(myeloids_cho, dims = 1:PCS2)

DimPlot(myeloids_cho, group.by = 'sample') + labs(title = 'Paper J.Cho')
DimPlot(myeloids_cho, group.by = 'Health') + labs(title = 'Paper J.Cho')

FeaturePlot(myeloids_cho, features = c('CD3E', 'MS4A1', 'TPSB2', 'C1QB'), order =T, cols = c('lightgray', 'red'))
dir.create('myeloids_cho_with_masts_filtering')
myeloids_cho <- resolutions(myeloids_cho, workingdir = 'myeloids_cho_with_masts_filtering', title = 'myeloids_cho_filtering')

saveRDS(myeloids_cho, file = 'myeloids_cho_with_masts_filtering.RDS')

# more filtering -----------------------------------------------------------------------

Idents(myeloids_cho) <- myeloids_cho$RNA_snn_res.1.5
myeloids_cho <- subset(myeloids_cho, idents = c(0:16, 18:20, 22:23))
DimPlot(myeloids_cho, label=T)

myeloids_cho <- seurat_to_pca(myeloids_cho)
setwd('~/data_Albas/GSE134809_RAW/')

PCS <- select_pcs(seudata_f, 2)
PCS2 <- select_pcs(seudata_f, 1.6)
p <- ElbowPlot(seudata_f, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0(' myeloids_cho_with_masts'))
ggsave(p, file = 'elbowplot_myeloids_cho_with_masts_more_filtering.jpeg')

myeloids_cho <- FindNeighbors(myeloids_cho, dims = 1:PCS2)
myeloids_cho <- RunUMAP(myeloids_cho, dims = 1:PCS2)

DimPlot(myeloids_cho, group.by = 'sample') + labs(title = 'Paper J.Cho')
DimPlot(myeloids_cho, group.by = 'Health') + labs(title = 'Paper J.Cho')

FeaturePlot(myeloids_cho, features = c('CD3E', 'MS4A1', 'TPSB2', 'C1QB'), order =T, cols = c('lightgray', 'red'))
dir.create('myeloids_cho_with_masts_more_filtering')
myeloids_cho <- resolutions(myeloids_cho, workingdir = 'myeloids_cho_with_masts_more_filtering', title = 'myeloids_cho_more_filtering')

saveRDS(myeloids_cho, file = 'myeloids_cho_with_masts_more_filtering.RDS')

## igsout --------------------------------------------------------------------------------
source('~/data_Albas/functions_scrseq.R')
myeloids_cho <- readRDS(file = 'myeloids_cho_with_masts_more_filtering.RDS')

gg <- rownames(myeloids_cho)[c(grep("^IGH",rownames(myeloids_cho)),
                               grep("^IGK", rownames(myeloids_cho)),
                               grep("^IGL", rownames(myeloids_cho)))]
genes <- setdiff(rownames(myeloids_cho),gg)
myeloids_cho <- subset(myeloids_cho,features = genes)
myeloids_cho

myeloids_cho <- seurat_to_pca(myeloids_cho)
setwd('~/data_Albas/GSE134809_RAW/')

PCS <- select_pcs(myeloids_cho, 2)
PCS2 <- select_pcs(myeloids_cho, 1.6)
p <- ElbowPlot(myeloids_cho, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0(' myeloids_cho_with_masts'))
ggsave(p, file = 'elbowplot_myeloids_cho_with_masts_more_filtering_igsout.jpeg')

myeloids_cho <- FindNeighbors(myeloids_cho, dims = 1:PCS2)
myeloids_cho <- RunUMAP(myeloids_cho, dims = 1:PCS2)

DimPlot(myeloids_cho, group.by = 'sample') + labs(title = 'Paper J.Cho')
DimPlot(myeloids_cho, group.by = 'Health') + labs(title = 'Paper J.Cho')

FeaturePlot(myeloids_cho, features = c('CD3E', 'MS4A1', 'TPSB2', 'C1QB'), order =T, cols = c('lightgray', 'red'))
FeaturePlot(myeloids_cho, features = c('percent.mt', 'nFeature_RNA', 'nCount_RNA'), order =T, cols = c('lightgray', 'red'))
dir.create('myeloids_cho_with_masts_more_filtering_igsout')
myeloids_cho <- resolutions(myeloids_cho, workingdir = 'myeloids_cho_with_masts_more_filtering_igsout', title = 'myeloids_cho_more_filtering_igsout')

saveRDS(myeloids_cho, file = 'myeloids_cho_with_masts_more_filtering_igsout.RDS')


# annotation -------------------------------------------------------------------------------------
setwd("~/000_GitHub/ibd-bcn_single_cell/GSE134809")
myeloids_cho <- readRDS(file = 'GSE134809_RAW/myeloids_cho_with_masts_more_filtering_igsout.RDS')
metadata <- data.frame(
        stringsAsFactors = FALSE,
                            clusters = c(0L,1L,2L,3L,4L,5L,6L,7L,8L,
                                         9L,10L,11L,12L,13L,14L,15L,16L,
                                         17L,18L,19L,20L,21L,22L),
                myeloids_cho_res_1.5 = c("Mast", "M2", "M0 Ribhi", "DCs CD1c", 
                                         "Mast", "M1", "Mast", "M2.2", "M0", 
                                         "M2.3", "DCs CCL19", "Mac Mthi", 
                                         "DCs Ribhi", "M1-M2", 
                                         "Inflammatory monocytes", "DCs Ribhi",
                                         "Macrophage CXCL10", "Macrophage MARCO",
                                         "Cycling monocytes", "DCs CCL19 RibHi",
                                         "M1.2", "Mast", "M1.3")
            )

myeloids_cho$updated_annotation <- mapvalues(myeloids_cho$RNA_snn_res.1.5,
                                     from = metadata$clusters,
                                     to = metadata$myeloids_cho_res_1.5)
DimPlot(myeloids_cho, group.by = 'updated_annotation', label=T)
saveRDS(myeloids_cho, file = 'GSE134809_RAW/myeloids_cho_with_masts_more_filtering_igsout.RDS')
Idents(myeloids_cho) <- myeloids_cho$updated_annotation
markers <- FindAllMarkers(myeloids_cho, only.pos = T,
                          min.pct = 0.25, thresh.use = 0.25)
markers <- markers[markers$p_val < 0.05,]
dat_list <- list()
for(i in 1:length(levels(markers$cluster))){
  genes <- as.character(markers$gene[markers$cluster == levels(markers$cluster)[i]])
  dat_list[[i]] <- genes
  names(dat_list)[i] <-  levels(markers$cluster)[i]
}

saveRDS(dat_list, file = 'GSE134809_RAW/markers_myeloids_cho.RDS')
saveRDS(myeloids_cho, 'GSE134809_RAW/myeloids_cho_with_masts_more_filtering_igsout.RDS')


# label transfer ---------------------------------------------------
library(matchSCore2)
library(nnet)
library(Matrix)
library(ggpubr)

setwd("~/000_GitHub/ibd-bcn_single_cell")
ref_a <- readRDS("Analysis of our data/02_Samples_Together/SUBSETS/ON_THEIR_OWN/myeloids_annotated.RDS")

cho_a <- readRDS('GSE134809/GSE134809_RAW/myeloids_cho_with_masts_more_filtering_igsout.RDS')

# gene_cl.ref: A named list of markers. 
ref_a@active.ident <- factor(ref_a$annotation_refined)
output.seu <- FindAllMarkers(ref_a, only.pos = T,  min.pct = 0.25, thresh.use = 0.25)
gene_cl.ref <- cut_markers(unique(output.seu$cluster),output.seu,ntop= 50)

cho_a@active.ident <- factor(cho_a$updated_annotation)
output.reg <- FindAllMarkers(cho_a, only.pos = T,  min.pct = 0.25, thresh.use = 0.25)
gene_cl <- cut_markers(unique(output.reg$cluster),output.reg,ntop=50)


### Training of the model  
ref_a <- ScaleData(ref_a, features = rownames(ref_a@assays$RNA@data))
mod <- train_model(scale.data = ref_a@assays$RNA@scale.data,
                   clus = ref_a@active.ident,
                   gene_cl.ref = gene_cl.ref,
                   prop = 0.75)

## Cell projection
cho_a <- ScaleData(cho_a, features = rownames(cho_a))
out <- identity_map(scale.data = cho_a$RNA@scale.data,
                    model = mod,
                    gene_cl.ref)

### cell identities
ids <- out$ids 
cho_a$label_transfer <- as.character(ids)
DimPlot(cho_a, group.by = 'label_transfer', label=T, repel=F)


ms <- matchSCore2(gene_cl.ref = gene_cl.ref[names(gene_cl.ref) %in% c('M2.2', 
                                                                      'M2',
                                                                      'M0',
                                                                      'IDA macrophage',
                                                                      'M1 ACOD1',
                                                                      'M1 CXCL5')],
                  gene_cl.obs = gene_cl[names(gene_cl) %in% c('M2', 
                                                              'M0 Ribhi', 
                                                              'M1',
                                                              'M2.2',
                                                              'M0',
                                                              'M2.3',
                                                              'M1-M2',
                                                              'Inflammatory monocytes',
                                                              'Macrophage CXCL10',
                                                              'Macrophage MARCO',
                                                              'M1.2',
                                                              'M1.3')], xlab = 'cho', ylab = 'Ours')

## The matchSCore heatmap is stored in the ggplot slot of ms. 
ms$ggplot + theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5),
        axis.text.y = element_text(size = 10), axis.title = element_text(size=8),
        legend.position = 'bottom')

