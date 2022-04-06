# Longitudinal single cell RNA-seq of the inflamed colon
# Mice were treated with 1.5%DSS for 6 days then switched to drinking water. 
# Colon from each time points (Day 0, 3, 6, 9, 12 and 15) or treatment 
# (PBS or Serpina3n) were collected. To limit mouse-to-mouse biased, 
# three biological replicates from each group were pooled and analyzed 
# using FACS-based, smart-seq2 method. The sequencing libraries were
# sequenced on the Illumina NextSeq500 platform.
library(ggplot2)
library(ggpubr)
library(viridis)
library(Seurat)
library(stringr)
library(data.table)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
source('~/data_Albas/functions_scrseq.R')


#
# Load data
#

setwd("~/000_GitHub/ibd-bcn_single_cell/GSE148794")
counts <- as.data.frame(fread('GSE148794_tc_ibd.count_table.tsv', sep = '\t'))
rownames(counts) <- counts$V1
counts <- counts[,-grep("V1", colnames(counts))]
counts <- as.matrix(counts)

#
# metadata 
#
metadata <- read.delim('GSE148794_tc_ibd.metadata.tsv', sep = '\t')
ref <- gsub('\\.','_',paste0('R19', str_split_fixed(pattern = '.19', str = metadata$Sample_name, n = 2)[,2]))
refok <- str_split_fixed(string = ref, pattern = '_S\\d', n = 2)[,1]
timepoint <- gsub(pattern = 'tc_ibd.',replacement = '',
                  x = str_split_fixed(pattern = '.19', 
                                      str = metadata$Sample_name,
                                      n = 2)[,1])
metadata$refok <- refok
metadata$timepoint <- timepoint
metadataok <- metadataoka.frame()
for(colname in colnames(counts)){
  ref <- substr(colname, 0,42)
  nn <-grep(ref, metametadataoka$refok)
  metadataok <- rbind(metadataok, metametadataoka[nn,])
  rownames(metadataok)[nrow(metadataok)] <- colname
}
saveRDS(metadataok, file = 'GSE148794_metadata_ok.RDS')
metadataok <- readRDS('GSE148794_metadata_ok.RDS')

#
# create seurat object
#

mice <- CreateSeuratObject(
  counts, min.features = 100,
  project = 'GSE148794',
  assay = "RNA",
  meta.data = metadataok
)

#
# percent MT
#
mice[["percent.mt"]] <- PercentageFeatureSet(object = mice, pattern = "^mt")

meta <- mice@meta.data
meta$density <- get_density(meta$percent.mt, meta$nFeature_RNA, n = 100)
jpeg(filename = 'percent_mt_density.jpeg', width = 1500, height = 1500, res = 150)
ggplot(meta) +
  geom_point(aes(percent.mt, nFeature_RNA, color = density)) +
  # geom_hline(yintercept = 100, color = 'gray', linetype="dashed") +
  # geom_vline(xintercept = 65, color = 'gray', linetype="dashed") +
  # geom_vline(xintercept = 25, color = 'gray', linetype="dashed") +
  scale_color_viridis() + theme_classic() + theme(plot.title = element_text(size = 25))+
  labs(title = 'GSE148794')
dev.off()

meta$density <- get_density(meta$nCount_RNA, meta$nFeature_RNA, n = 100)
jpeg(filename = 'count_feature_density.jpeg', width = 1500, height = 1500, res = 150)
ggplot(meta) +
  geom_point(aes(nCount_RNA, nFeature_RNA, color = density)) +
  # geom_hline(yintercept = 100, color = 'gray', linetype="dashed") +
  # geom_vline(xintercept = 65, color = 'gray', linetype="dashed") +
  # geom_vline(xintercept = 25, color = 'gray', linetype="dashed") +
  scale_color_viridis() + theme_classic() + theme(plot.title = element_text(size = 25))+
  labs(title = 'GSE148794')
dev.off()

# I do not remove cells by now.

#
# sample processing
#

mice <- seurat_to_pca(mice)
PCS <- select_pcs(mice, 2)
PCS2 <- select_pcs(mice, 1.6) # 41
p <- ElbowPlot(mice, ndims = 100) +
  geom_vline(xintercept = PCS,
             linetype = 'dotdash') +
  geom_vline(xintercept = PCS2, linetype = 'dashed',
             colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=8, 
           label= paste("sdev > 2;\nPCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=10, label= paste("sdev > 1.6;\nPCs =", PCS2),
           color="#BB0000")+
  labs(title = 'GSE148794')
ggsave(p, file = 'elbowplot.jpeg')

mice <- FindNeighbors(mice, dims = 1:PCS2)
mice <- RunUMAP(mice, dims = 1:PCS2)

DimPlot(mice, group.by = 'timepoint') + labs(title = 'GSE148794_together')

mice
# An object of class Seurat 
# 27407 features across 14634 samples within 1 assay 
# Active assay: RNA (27407 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap

dir.create('together_noharmony')
mice <- resolutions(mice, workingdir = 'together_noharmony/',
                    title = 'GSE148794_together_noharmony_')
saveRDS(mice, file = 'together_noharmony/mice.RDS')


library(harmony)
mice <- RunHarmony(mice, group.by = 'Sample_name', dims.use = 1:PCS2)
ElbowPlot(mice, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 50, linetype = 2) +
  labs(title = paste0('Harmony - 41pcs'))
mice<- FindNeighbors(mice, reduction = "harmony", dims = 1:50)
mice<-RunUMAP(mice, dims=1:50, reduction= "harmony")

DimPlot(mice, group.by = 'Sample_name') + labs(title = 'Mice - Harmony 41-50')

dir.create('mice_harmony_41_50')
mice <- resolutions(mice,
                   workingdir = 'mice_harmony_41_50',
                   title = 'mice_harmony_41_50')

saveRDS(mice, file = 'mice_harmony_41_50/mice_harmony_41_50.RDS')


#
# subsetting data --------------------------------------------------------------
#

mice <- readRDS('together_noharmony/mice.RDS')

anotation <- data.frame(
  stringsAsFactors = FALSE,
                 cluster_number = c(0L, 1L,
                                    2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,
                                    13L,14L,15L,16L,17L,18L,19L,20L,
                                    21L,22L,23L,24L,25L,26L,27L),
                   cluster_name = c("plasmas","myeloid","stroma","tcells","tcells",
                                    "plasmas","stroma","myeloid","myeloid","stroma",
                                    "epi","myeloid","myeloid","epi",
                                    "stroma","cycling","myeloid","epi","myeloid",
                                    "stroma","plasmas","stroma","epi",
                                    "stroma","stroma","stroma","epi","epi")
             )

mice$subset <- plyr::mapvalues(x = mice$RNA_snn_res.0.7, 
                               from = anotation$cluster_number, 
                               to= anotation$cluster_name)
DimPlot(mice, group.by = 'subset')
saveRDS(mice, file = 'together_noharmony/mice.RDS')

##
## epithelium ------------------------------------------------------------------
##
dir.create('~/000_GitHub/ibd-bcn_single_cell/GSE148794/epithelium')
setwd('~/000_GitHub/ibd-bcn_single_cell/GSE148794/epithelium')

epi <- mice[,mice$subset == 'epi'] # 1367 

epi <- seurat_to_pca(epi)
DimPlot(epi, group.by = 'timepoint', reduction = 'pca')
PCS <- select_pcs(epi, 2)
PCS2 <- select_pcs(epi, 1.6) # 41
p <- ElbowPlot(epi, ndims = 100) +
  geom_vline(xintercept = PCS,
             linetype = 'dotdash') +
  geom_vline(xintercept = 50,
             linetype = 'dashed', color = 'gray') +
  geom_vline(xintercept = PCS2, linetype = 'dashed',
             colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=8, 
           label= paste("sdev > 2;\nPCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=10, label= paste("sdev > 1.6;\nPCs =", PCS2),
           color="#BB0000")+
  labs(title = 'GSE148794 epi')

ggsave(p, file = 'elbowplot.jpeg')

epi <- FindNeighbors(epi, dims = 1:36)
epi <- RunUMAP(epi, dims = 1:36)

DimPlot(epi, group.by = 'timepoint') + labs(title = 'GSE148794_epi')

epi
# An object of class Seurat 
# 27407 features across 1367 samples within 1 assay 
# Active assay: RNA (27407 features, 2000 variable features)
# 3 dimensional reductions calculated: pca, umap, harmony

dir.create('epi_noharmony')
epi <- resolutions(epi, workingdir = 'epi_noharmony/',
                    title = 'GSE148794_epi_noharmony_')
saveRDS(epi, file = 'epi_noharmony/epi.RDS')

library(harmony)
epi <- RunHarmony(epi, group.by = 'Sample_name', dims.use = 1:36)
ElbowPlot(epi, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 40, linetype = 2) +
  labs(title = paste0('Harmony - 36pcs'))
epi<- FindNeighbors(epi, reduction = "harmony", dims = 1:40)
epi<-RunUMAP(epi, dims=1:40, reduction= "harmony")

DimPlot(epi, group.by = 'Sample_name', legend = F)
+ labs(title = 'epi - Harmony 36-40')

dir.create('epi_harmony_36_40')
epi <- resolutions(epi,
                    workingdir = 'epi_harmony_36_40',
                    title = 'epi_harmony_36_40')

saveRDS(epi, file = 'epi_harmony_36_41/epi_harmony_36_41.RDS')


##
## stroma ------------------------------------------------------------------
##
dir.create('~/000_GitHub/ibd-bcn_single_cell/GSE148794/stroma')
setwd('~/000_GitHub/ibd-bcn_single_cell/GSE148794/stroma')

stroma <- mice[,mice$subset == 'stroma'] #  3603

stroma <- seurat_to_pca(stroma)
DimPlot(stroma, group.by = 'timepoint', reduction = 'pca')
PCS <- select_pcs(stroma, 2)
PCS2 <- select_pcs(stroma, 1.6) # 41
p <- ElbowPlot(stroma, ndims = 100) +
  geom_vline(xintercept = PCS,
             linetype = 'dotdash') +
  geom_vline(xintercept = 50,
             linetype = 'dashed', color = 'gray') +
  geom_vline(xintercept = PCS2, linetype = 'dashed',
             colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=8, 
           label= paste("sdev > 2;\nPCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=10, label= paste("sdev > 1.6;\nPCs =", PCS2),
           color="#BB0000")+
  labs(title = 'GSE148794 stroma')

ggsave(p, file = 'elbowplot.jpeg')

stroma <- FindNeighbors(stroma, dims = 1:41)
stroma <- RunUMAP(stroma, dims = 1:41)

DimPlot(stroma, group.by = 'timepoint') + labs(title = 'GSE148794_stroma')

stroma
# An object of class Seurat 
# 27407 features across 3603 samples within 1 assay 
# Active assay: RNA (27407 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap

dir.create('stroma_noharmony')
stroma <- resolutions(stroma, workingdir = 'stroma_noharmony/',
                   title = 'GSE148794_stroma_noharmony_')
saveRDS(stroma, file = 'stroma_noharmony/stroma.RDS')

library(harmony)
stroma <- RunHarmony(stroma, group.by = 'Sample_name', dims.use = 1:41)
ElbowPlot(stroma, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 41, linetype = 2) +
  labs(title = paste0('Harmony - 41pcs'))
stroma<- FindNeighbors(stroma, reduction = "harmony", dims = 1:41)
stroma<-RunUMAP(stroma, dims=1:41, reduction= "harmony")

DimPlot(stroma, group.by = 'Sample_name', legend = F)
+ labs(title = 'stroma - Harmony 41-41')

dir.create('stroma_harmony_41_41')
stroma <- resolutions(stroma,
                   workingdir = 'stroma_harmony_41_41',
                   title = 'stroma_harmony_41_41')

saveRDS(stroma, file = 'stroma_harmony_41_41/stroma_harmony_41_41.RDS')

##
## myeloids ------------------------------------------------------------------
##
dir.create('~/000_GitHub/ibd-bcn_single_cell/GSE148794/myeloids')
setwd('~/000_GitHub/ibd-bcn_single_cell/GSE148794/myeloids')

myeloids <- mice[,mice$subset == 'myeloid'] # 4585 

myeloids <- seurat_to_pca(myeloids)
DimPlot(myeloids, group.by = 'timepoint', reduction = 'pca')
PCS <- select_pcs(myeloids, 2) # 18
PCS2 <- select_pcs(myeloids, 1.6) # 31
p <- ElbowPlot(myeloids, ndims = 100) +
  geom_vline(xintercept = PCS,
             linetype = 'dotdash') +
  geom_vline(xintercept = PCS2, linetype = 'dashed',
             colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=8, 
           label= paste("sdev > 2;\nPCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=10, label= paste("sdev > 1.6;\nPCs =", PCS2),
           color="#BB0000")+
  labs(title = 'GSE148794 myeloids')

ggsave(p, file = 'elbowplot.jpeg')

myeloids <- FindNeighbors(myeloids, dims = 1:31)
myeloids <- RunUMAP(myeloids, dims = 1:31)

DimPlot(myeloids, group.by = 'timepoint') + labs(title = 'GSE148794_myeloids')

myeloids
# An object of class Seurat 
# 27407 features across 4585 samples within 1 assay 
# Active assay: RNA (27407 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap

dir.create('myeloids_noharmony')
myeloids <- resolutions(myeloids, workingdir = 'myeloids_noharmony/',
                      title = 'GSE148794_myeloids_noharmony_')
saveRDS(myeloids, file = 'myeloids_noharmony/myeloids.RDS')

library(harmony)
myeloids <- RunHarmony(myeloids, group.by = 'Sample_name', dims.use = 1:31)
ElbowPlot(myeloids, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 30, linetype = 2) +
  labs(title = paste0('Harmony - 31pcs'))
myeloids<- FindNeighbors(myeloids, reduction = "harmony", dims = 1:30)
myeloids<-RunUMAP(myeloids, dims=1:30, reduction= "harmony")

DimPlot(myeloids, group.by = 'Sample_name') + NoLegend() + labs(title = 'myeloids - Harmony 31_30')

dir.create('myeloids_harmony_31_30')
myeloids <- resolutions(myeloids,
                      workingdir = 'myeloids_harmony_31_30',
                      title = 'myeloids_harmony_31_30')

saveRDS(myeloids, file = 'myeloids_harmony_31_30/myeloids_harmony_31_30.RDS')


