# Title	          Intra- and Inter-cellular Rewiring of the Human Colon during Ulcerative Colitis
#

source('~/000_GitHub/ibd-bcn_single_cell/source/functions_scrnaseq.R')

#
# Load data:
#

reg <- readRDS('~/single_cell/Regev/data/train.Imm.seur.rds')
myeloid_reg <- reg[,reg$Cluster %in% c('Macrophages',
                                       'DC1',
                                       'DC2',
                                       'Inflammatory Monocytes',
                                       'CD69+ Mast',
                                       'CD69- Mast',
                                       'Cycling Monocytes')]

# percent.mt filtering
myeloid_reg[["percent.mt"]] <- PercentageFeatureSet(myeloid_reg, pattern = "^MT-")
myeloid_reg <- myeloid_reg[, myeloid_reg$percent.mt < 25]

# PCA and UMAP calculation
myeloid_reg <- seurat_to_pca(myeloid_reg)

ElbowPlot(myeloid_reg, ndims = 100) +
  geom_vline(xintercept = 22, colour="#BB0000")+
  annotate(geom="text", x=22-5, y=4, label= paste("sdev > 1.6; PCs =", 22),
           color="#BB0000")+
  labs(title = paste0('Seurat- merge all myeloid regev'))

myeloid_reg <- FindNeighbors(myeloid_reg, dims = 1:22)
myeloid_reg <- RunUMAP(myeloid_reg, dims = 1:22)
Idents(myeloid_reg) <- myeloid_reg$Cluster

myeloid_reg <- resolutions(myeloid_reg, 
                           workingdir = getwd(), 
                           title = 'Myeloids_Regev')

myeloid_reg <- FindClusters(myeloid_reg, resolution = 0.7)

# filtered - cluster 14 from resolution 0.7 is stroma
myeloid_reg <- myeloid_reg[,myeloid_reg$RNA_snn_res.0.7 != '14']

# reanalyisis
myeloid_reg <- seurat_to_pca(myeloid_reg)

ElbowPlot(myeloid_reg, ndims = 100) +
  geom_vline(xintercept = 21, colour="#BB0000")+
  annotate(geom="text", x=21-5, y=4, label= paste("sdev > 1.6; PCs =", 21),
           color="#BB0000")

myeloid_reg <- FindNeighbors(myeloid_reg, dims = 1:21)
myeloid_reg <- RunUMAP(myeloid_reg, dims = 1:21)

# igs out -----------------
gg <- rownames(myeloid_reg)[c(grep("^IGH",rownames(myeloid_reg)),
                              grep("^IGK", rownames(myeloid_reg)),
                              grep("^IGL", rownames(myeloid_reg)))]
genes <- setdiff(rownames(myeloid_reg),gg)
myeloid_reg <- subset(myeloid_reg,features = genes)

myeloid_reg <- seurat_to_pca(myeloid_reg)

ElbowPlot(myeloid_reg, ndims = 100) +
  geom_vline(xintercept = 21, colour="#BB0000")+
  annotate(geom="text", x=21-5, y=4, label= paste("sdev > 1.6; PCs =", 21),
           color="#BB0000")

myeloid_reg <- FindNeighbors(myeloid_reg, dims = 1:21)
myeloid_reg <- RunUMAP(myeloid_reg, dims = 1:21)

DimPlot(myeloid_reg, group.by = 'Health') + labs(title = '')

FeaturePlot(myeloid_reg, features = c('CD3E', 'C1QA', 'TPSAB1', 'C1QB'), order =T, cols = c('lightgray', 'red'))
FeaturePlot(myeloid_reg, features = c('percent.mt', 'nFeature_RNA', 'nCount_RNA'), order =T, cols = c('lightgray', 'red'))
dir.create('myeloids_filtering_igsout')
myeloid_reg <- resolutions(myeloid_reg, 
                           workingdir = 'myeloids_filtering_igsout', 
                           title = 'myeloids_filtering_igsout')

saveRDS(myeloid_reg, file = 'myeloids_filtering_igsout.RDS')

# annotation ---------------------------------------------------------------
library(readr)
metadata <- data.frame(
  stringsAsFactors = FALSE,
  clusters = c(0L,1L,2L,3L,4L,5L,6L,7L,8L,
               9L,10L,11L,12L,13L,14L,15L,16L,
               17L),
  myeloids_Regev_res_1.3 = c("M2","Mast 1","DCs CD1c","M0",
                             "Mast 2","M0_2",
                             "Inflammatory monocytes-mac","Mast Ribhi","M2-like Ly96",
                             "M0 Ribhi","Cycling monocytes",
                             "Mac Mthi","Mac IFN pathway","DCs CD1c 2",
                             "Mast 3","DCs IDO1","Mast 4",
                             "M1 CXCL3")
)

myeloid_reg$annotation <- mapvalues(myeloid_reg$RNA_snn_res.1.3,
                                    from = metadata$clusters,
                                    to = metadata$myeloids_Regev_res_1.3)

#label transfer with our macrophages ----------------------------------------

library(matchSCore2)
library(nnet)
library(Matrix)

ref_a <- readRDS("~/000_GitHub/ibd-bcn_single_cell/Analysis of our data/02_Samples_Together/SUBSETS/ON_THEIR_OWN/myeloids_annotated.RDS")

regev_a <- readRDS('~/data_Albas/regev_myeloids/myeloids_regev_with_masts_more_filtering_igsout.RDS')

ref_a@active.ident <- factor(ref_a$annotation_refined)
output.seu <- FindAllMarkers(ref_a, only.pos = T)

gene_cl.ref <- cut_markers(unique(output.seu$cluster),output.seu,ntop= 50)

regev@active.ident <- factor(regev$annotation)
regev_a@active.ident <- factor(regev_a$annotation)
output.reg <- FindAllMarkers(regev_a, only.pos = T)


gene_cl <- cut_markers(unique(output.reg$cluster),output.reg,ntop=50)


### Training of the model  
ref <- ScaleData(ref, features = rownames(ref@assays$RNA@data))
mod <- train_model(scale.data = ref@assays$RNA@scale.data,
                   clus = ref@active.ident,
                   gene_cl.ref = gene_cl.ref,
                   prop = 0.75)

## Cell projection
regev_a <- ScaleData(regev_a, features = rownames(regev_a))
out <- identity_map(scale.data = regev_a$RNA@scale.data,
                    model = mod,
                    gene_cl.ref)

### cell identities
ids <- out$ids 
regev$label_transfer <- as.character(ids)
DimPlot(regev, group.by = 'label_transfer', label=T, repel=F)


ms <- matchSCore2(gene_cl.ref = gene_cl.ref[names(gene_cl.ref) %in% c('M2.2', 
                                                                      'M2',
                                                                      'M0',
                                                                      'IDA macrophage',
                                                                      'M1 ACOD1',
                                                                      'M1 CXCL5')],
                  gene_cl.obs = gene_cl[names(gene_cl) %in% c('M2', 
                                                              'M0',
                                                              'M0_2',
                                                              'Inflammatory monocytes-mac',
                                                              'M2-like Ly96',
                                                              'Mac IFN pathway',
                                                              'M1 CXCL3')],
                  xlab = 'Regev', ylab = 'Ours')

## The matchSCore heatmap is stored in the ggplot slot of ms. 
ms$ggplot + theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5),
        axis.text.y = element_text(size = 10), axis.title = element_text(size=8),
        legend.position = 'bottom')
regev$annotation <- droplevels.factor(regev$annotation)
table(regev$label_transfer, regev$annotation)

saveRDS(ms, file = 'ms_analysis.RDS')
