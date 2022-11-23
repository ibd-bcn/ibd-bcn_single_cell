
# Genesets from the following datasests are obtained to test
# our own macrophage data. We generate different csv files for each analysis
#
# M-MO & GM-MO with LPS: GSE156921 y GSE99056
# 
# Ref:J Innate Immun. 2022;14(3):243-256.
# 
# M-MO with 5HT: GSE121825, GSE161774 y GSE94608
# 
# Ref: J Immunol. 2020 May 15;204(10):2808-2817; Sci Rep 2017 Nov 7;7(1):14761.
# 
# siRNA for MAF or MAFB in M-MO: GSE155719
# 
# Ref: Front Immunol. 2020 Nov 18;11:603507.

library(Seurat)
library(readr)
source('~/000_GitHub/ibd-bcn_single_cell/source/colors.R')
myeloids <- readRDS('~/000_GitHub/ibd-bcn_single_cell/Analysis of our data/02_Samples_Together/SUBSETS/ON_THEIR_OWN/myeloids_annotated.RDS')

data.query <- myeloids[,myeloids$annotation_intermediate %in% c("IDA macrophage",
                                                                "M0"  ,  "M1" ,
                                                                "Inflammatory monocytes", "M2",
                                                                "DCs", "Cycling myeloid")]
DimPlot(data.query, group.by = 'annotation_refined', 
        cols = cols_myeloids[unique(data.query$annotation_refined)]) + labs(title = 'keepers') 


## modules ------------------------------------------------------------------------

files <- list.files(recursive = T)
files <- files[-grep('R$', files)] # avoid the code
names <- stringr::str_split_fixed(gsub(pattern = '.csv', replacement = '', x= files), pattern = "/", n = 2)[,2]
i = 0
list_modules <- vector(mode = 'list', length = length(names))
names(list_modules) <- names
for(file in files){
  i = i+1
  fi <- read_csv(file, col_names = FALSE)
  list_modules[[i]] <- fi$X1 
}

data.query <- AddModuleScore(data.query, features = list_modules, name = 'Module', search = T)

data.query@active.ident <- data.query$annotation_refined
saveRDS(data.query, file = 'data.query.RDS')


# Plot scores ------------------------------------------------------------------

a <- DimPlot(data.query, label = TRUE, repel = TRUE) + labs(title='Clusters') +
  NoLegend()

plots <- vector('list', length = length(names(list_modules)))
for(i in 1:length(names(list_modules))){
  b <- FeaturePlot(data.query,
                   features = paste("Module", i, sep=''),
                   label = FALSE, repel = TRUE, order =T) +
    scale_colour_gradientn(colours = rev(
      RColorBrewer::brewer.pal(n = 11, name = "RdBu"))) +
    labs(title = names(list_modules)[i]) +
    theme(title = element_text(size = 8))
  
  plots[[i]] <- b
}

names(plots) <- names(list_modules)
dir.create('PLOTS')
for(i in 1:length(plots)){
  ggsave(filename = paste(names(plots)[i], '.jpeg', sep=''),
         plot = plots[[i]], device = 'jpeg',path = 'PLOTS', dpi = 300
  )
}
