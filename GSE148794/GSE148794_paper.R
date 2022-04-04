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
