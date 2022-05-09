cols <- c("#A6CEE3","#1F78B4", 
          "#B2DF8A", "#33A02C",
          "#FB9A99", "#E31A1C",
          "#FDBF6F","#FF7F00",
          "#CAB2D6", "#6A3D9A", 
          "#FFFF99", "#B15928", 
          "#1B9E77", "#D95F02",
          "#7570B3", "#E7298A", 
          "#66A61E", "#E6AB02", 
          "#A6761D", "#666666", 
          "aquamarine", "aquamarine4", 
          "beige", "bisque4",
          "lightpink", "lightpink3",
          "slateblue1", "slateblue4", 
          "tan1", "tan4", "purple4", 
          "purple1", "lightsalmon4",
          "lightsalmon","maroon4",
          "maroon1", "lemonchiffon4", 
          "lemonchiffon")
           
colors_subset <- readRDS('source/colors_subset.RDS')
colors_health <- readRDS('source/colors_health.RDS')
cols_epi <- readRDS('source/cols_epi.RDS')
cols_myeloids <- readRDS('source/cols_myeloids.RDS')
cols_plasmas <- readRDS('source/cols_plasmas.RDS')
cols_stroma <- readRDS('source/cols_stroma.RDS')
cols_tcells <- readRDS('source/cols_tcells.RDS')
cols_myeloids_intermediate <- readRDS('source/cols_myeloids_intermediate.RDS')

colors_health_alpha <- c('HC' =  '#71AD9680',
                         'CDa' = '#52517480',
                         'UCa' = '#CBA55280')


