# UCB + MatureN cells

# UCB
# https://pubmed.ncbi.nlm.nih.gov/31049560/
#   Single-cell transcriptomic landscape of nucleated cells in umbilical cord blood - PubMed

# matureN cells
# paper: https://academic.oup.com/nsr/article/8/3/nwaa180/5896476

library(Seurat)
library(stringr)
source('~/data_Albas/functions_scrseq.R')

# UCB analysis -----------------------------------------------------------------
fil2 <- read.delim(file = '~/data_Albas/UCB_paper/fig1_scale_exp.tsv', sep = '\t')
fil <- read.delim(file = '~/data_Albas/UCB_paper/fig1_umi.tsv', sep = '\t')
coln <- as.data.frame(str_split_fixed(string = gsub("\\.", '_', colnames(fil)), 
                                      pattern = '_', n = 2))
colnames(coln) <- c('barcode', 'number')
rownames(coln) <- colnames(fil)

ucb <- CreateSeuratObject(
  fil, min.features = 100,
  project = 'UCB',
  assay = "RNA",
  meta.data = coln
)

ucb[["percent.mt"]] <- PercentageFeatureSet(object = ucb, pattern = "^MT")

ucb <- seurat_to_pca(ucb)

PCS <- select_pcs(ucb, 2)
PCS2 <- select_pcs(ucb, 1.6)
p <- ElbowPlot(ucb, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  geom_vline(xintercept = 25, colour="gray", linetype = 2)+
  annotate(geom="text", x=PCS-5, y=20, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=14, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('UCB_paper'))
ggsave(p, file = '~/data_Albas/UCB_paper/elbowplot.jpeg')

ucb <- FindNeighbors(ucb, dims = 1:25)
ucb <- RunUMAP(ucb, dims = 1:25)

DimPlot(ucb, group.by = 'number')
FeaturePlot(ucb, features = c('CMTM2', 'FCGR3B', 'PROK2', 'IFIT2'))

dir.create('~/data_Albas/UCB_paper/resolutions_umi')
ucb <- resolutions(ucb, 
                   workingdir = '~/data_Albas/UCB_paper/resolutions_umi', 
                   title = 'UCB_paper')
saveRDS(ucb, file = '~/data_Albas/UCB_paper/resolutions_umi/ucb_res.RDS') 

# matureN ----------------------------------------------------------------------
setwd('~/data_Albas/nwaa180/GSE149938_Human Blood Cells/')

mat <- read.table('GSE149938_umi_matrix.csv', sep=',')
mat_t <- t(mat)
colnames(mat_t)[1:3]

metadata <- str_split_fixed(colnames(mat_t), "_", 3)
rownames(metadata) <- colnames(mat_t)
colnames(metadata) <- c('annotation', 'patient', 'meh')
metadata <- as.data.frame(metadata)

seu_mydat_blood <- CreateSeuratObject(
  mat_t, min.features = 100,
  project = 'Blood',
  assay = "RNA",
  meta.data = metadata
)
seu_mydat_blood[["percent.mt"]] <- PercentageFeatureSet(object = seu_mydat_blood, pattern = "^MT")

seu <- seu_mydat_blood
seu <- seurat_to_pca(seu)

PCS <- select_pcs(seu, 2)
PCS2 <- select_pcs(seu, 1.6)
p <- ElbowPlot(seu, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('NWAA180'))
ggsave(p, file = '~/data_Albas/nwaa180/elbowplot.jpeg')

dir.create('~/data_Albas/nwaa180/resolutions')
seu <- FindNeighbors(seu, dims = 1:PCS)
seu <- RunUMAP(seu, dims = 1:PCS)
seu <- resolutions(seu, 
                   workingdir = '~/data_Albas/nwaa180/resolutions', 
                   title = 'nwaa180')
saveRDS(seu, file = '~/data_Albas/nwaa180/seurat_res.RDS') 

B <- read.table('GSM4793029_B.txt.gz')
NK <- read.table('GSM4793030_NK.txt.gz')
t <- read.table('GSM4793031_T.txt.gz')
Mo <- read.table('GSM4793032_Mo.txt.gz')
Neu<- read.table('GSM4793033_Neu.txt.gz')
Ery<- read.table('GSM4793034_Ery.txt.gz')
BM <- read.table('~/data_Albas/nwaa180/GSE137864_Bone Marrow/barcodes.tsv.gz')

intersect(colnames(seu), B$V1)

seu$annotation_group <- ifelse(colnames(seu) %in% B$V1, yes = 'B', no = 'BM')
seu$annotation_group <- ifelse(colnames(seu) %in% NK$V1, yes = 'NK', no = seu$annotation_group)
seu$annotation_group <- ifelse(colnames(seu) %in% t$V1, yes = 'T', no = seu$annotation_group)
seu$annotation_group <- ifelse(colnames(seu) %in% Mo$V1, yes = 'Mo', no = seu$annotation_group)
seu$annotation_group <- ifelse(colnames(seu) %in% Neu$V1, yes = 'Neu', no = seu$annotation_group)
seu$annotation_group <- ifelse(colnames(seu) %in% Ery$V1, yes = 'Ery', no = seu$annotation_group)
seu$annotation_group <- ifelse(colnames(seu) %in% BM$V1, yes = 'BM', no = seu$annotation_group)

seu$annotation_2 <- paste(seu$annotation_group, seu$annotation, sep = '_')
DimPlot(seu, group.by = 'annotation', label=T, repel = T)
DimPlot(seu, group.by = 'annotation_2', label=T, repel = T)

saveRDS(seu, file = '~/data_Albas/nwaa180/seurat_res_anot.RDS') 
Idents(seu) <- seu$annotation
markers_NeuMature <- FindMarkers(seu, ident.1 = 'matureN',  only.pos = TRUE,
                                 min.pct = 0.25, 
                                 thresh.use = 0.25)

# Heatmap ----------------------------------------------------------------------

genelist <- c("S100A8",   "S100A9" ,    "IFITM2" , "G0S2", 'SRGN',
              "SAT1", "NAMPT", "LITAF" ,  "SOD2",  "AQP9"  ,
              "PLEK" ,     "IVNS1ABP" , 'CXCL8', 'BCL2A1',
              "LTF" , "CYP4F3" ,"LCN2" ,"CAMP" ,"DEFA3" ,  "PADI4" , 
              "PGLYRP1",  "ANXA3" ,"MMP8"  , "CRISP3"  , 
              "PMS2"   ,  "BPI"    ,   "CEACAM8" ,"CD24", 
              "CHI3L1",  "HIST2H2BE" , "PHOSPHO1",
              'FTH1', 'IL1B', 'IL1RN',"PDE4B",
              "HCAR3", "PPIF" ,    "HCAR2",    "GK" ,
              "ACTB", "CD55" ,     "PTGS2" ,    "GCA"   , 
              "SOCS3"  ,   "NFKBIZ"  )

library(ComplexHeatmap)
library(circlize)
mat1 <- ucb
mat1 <- t(FetchData(mat1, vars = genelist))
mat1 <- matrix(data = rowMeans(mat1), ncol = 1, dimnames = list(rownames(mat1), 
                                                                'UCB Neutrophil'))

mat1a <- seu
mat1a <- t(FetchData(mat1a, vars = genelist))
mat1a <- matrix(data = rowMeans(mat1a), ncol = 1, 
                dimnames = list(rownames(mat1a),
                                'BM matureN'))

mat1 <- cbind(mat1, mat1a)

mat2 <- myeloids[,myeloids$annotation_refined %in% c('Neutrophil 1',
                                                         'Neutrophil 2',
                                                         'Neutrophil 3')]
mat2 <- FetchData(mat2, vars = c(genelist, 'annotation_refined'))
mat2 <- mat2[order(mat2$new_annotation_refined),]
mat_annot <- c('Neutrophil 1', 'Neutrophil 2', 'Neutrophil 3')
mat2 <- aggregate(mat2[, 1:(ncol(mat2)-1)], list(mat2$new_annotation_refined), mean)
mat2 <- t(mat2[,2:ncol(mat2)])
colnames(mat2) <- c('Neutrophil 1', 'Neutrophil 2', 'Neutrophil 3')


ha1 = HeatmapAnnotation(cluster = rep('UCB Neutrophil', ncol(mat1)), 
                        show_annotation_name = F, 
                        col = list(
                          cluster = c('UCB Neutrophil' =  '#71AD96',
                                      'BM matureN' = '#D00CC9')))

col_rnorm = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_runif = colorRamp2(c(0, 3), c("white", "orange"))
col_fun <- colorRamp2(c(-2, 0.5,2), c("#3288BD", "#FFFFBF", "#D53E4F"))

ht1 = Heatmap(mat1, name = "col_fun", top_annotation = ha1,
              column_order = colnames(mat1), row_order = rownames(mat1),
              row_names_gp = gpar(col = c("black"), fontsize = c(7)),
              column_names_gp = gpar(col = c("white"), fontsize = c(0)),
              col = col_fun, column_title = "UCB", row_km = 6)

ha2 = HeatmapAnnotation(cluster = mat_annot, show_annotation_name = F,
                        col = list(
                          cluster = c('Neutrophil 1' = "#A6CEE3",
                                      'Neutrophil 2' = "#1F78B4", 
                                      'Neutrophil 3' = "#B2DF8A")))

ht2 = Heatmap(mat2, name = "col_fun", col = col_fun, 
              clustering_distance_rows = 'canberra',
              top_annotation = ha2, column_order = colnames(mat2),  
              row_names_gp = gpar(col = c("black"), fontsize = c(7)),
              column_names_gp = gpar(col = c("white"), fontsize = c(0)),
              column_title = "ours")
ht_list = ht2 + ht1 

draw(ht_list, cluster_rows=T,row_km = 6)

