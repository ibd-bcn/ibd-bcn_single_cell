# Functions used in the analysis

#
# seurat object - steps to PCA
#

# standardized method for processing the samples to obtain the PCA reduction
seurat_to_pca <- function(object){
  print('NORMALIZATION')
  object <- NormalizeData(object, normalization.method = "LogNormalize")
  print('VARIABLE FEATURES')
  object <- FindVariableFeatures(object, nfeatures = 2000)
  print('SCALE DATA')
  object <- ScaleData(object)
  print('RUN PCA')
  object <- RunPCA(object, npcs = 100, ndims.print = 1, nfeatures.print = 5)
  return(object)
}


#
# Get the resolutions
#
# This function tests different resolutions, saves the markers as csv files
# and generates a cluster tree with the different resolutions tested.

resolutions <- function(object,
                        resolutions = c(0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5), 
                        workingdir,
                        title = object$orig.ident[1],
                        res = 300,
                        width = 25,
                        height = 20){
  require(clustree)
  require(Seurat)
  if(!('Seurat' %in% is(object))){stop('object is not Seurat')}
  oldworkingdir <- getwd()
  setwd(workingdir)
  for(i in 1:length(resolutions)){
    object <- FindClusters(object, resolution = resolutions[i])
    
    markers <- FindAllMarkers(object = object,
                              only.pos = TRUE,
                              min.pct = 0.25, 
                              thresh.use = 0.25)
    
    write.table(x = markers, file = paste0(title,'_markers_resolution_',resolutions[i],'.csv'),
                row.names = F, sep = ';', dec= ',', col.names = T)
    
    png(filename = paste0(title,'_resolution_',resolutions[i],'.png'),
        res = res, width = width, height = height, units = 'cm')
    p <- DimPlot(object, label = T) +
      theme(legend.position = 'right') + labs(title = title)
    print(p)
    dev.off()
    

    dev.off()
  }
  
  png(filename = paste0(title,'_clustree_resolutions.png'),
      res = res, width = width+5, height = height+5, units = 'cm')
  p <- clustree(object, prefix = 'RNA_snn_res.') + labs(subtitle = title )
  print(p)
  dev.off()
  
  setwd(oldworkingdir)
  
  rm(markers, oldworkingdir, p)
  return(object)
}

#
# Custom violin plot and ridge plot functions
#
# designed for seurat objects.

violin_plot <- function(object, features,
                        split.by = NULL, 
                        idents = NULL,
                        colors=NULL,
                        group.by=NULL) {
  require(ggplot2)
  require(Seurat)
  require(patchwork)
  require(formulaic)
  require(gplots)
  require(randomcoloR)
  require(scales)
  require(ggridges)
  require(cowplot)
  
  if(is.null(colors)){nucolor <- T}else{nucolor <- F}
  
  if (sum(features %in% rownames(object)) != length(features)){
    print(paste('Feature(s)', features[!(features %in% rownames(object))], 'not present in the object'))
    features <- features[features %in% rownames(object)]
  }
  
  if(!is.null(idents)){
    if (!is.null(x = group.by)) {
      object <- SetIdent(object, value = group.by)
      object <- subset(object, idents = idents)
    }else{
      object <- subset(object, idents = idents)
    }
    }
  
  #
  # select cells
  #
  if (is.null(x = idents)) {
    cells <- colnames(x = object)
  } else {
    cells <- names(x = Idents(object = object)[Idents(object = object) %in% idents])
  }
  
  idents <- if (is.null(x = group.by)) {
    Idents(object = object)[cells]
  } else {
    object[[group.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = idents)) {
    idents <- droplevels(factor(x = idents))
  }
  
  idents <- droplevels(factor(x = idents, levels = levels(idents)[order(levels(idents))]))
  # print(idents)
  if (is.null(x = colors)) {
    colors <- scales::hue_pal()(length(x = levels(x = idents)))
    # colors <- alpha(colors, alpha = 0.5)
  } else {
    if(length(colors) == 1) {
      colors <- rep(colors, length(levels(idents)))
    }
    if(length(colors) < length(levels(idents))){
      ll <- length(levels(idents)) - length(colors)
      colors <- c(colors, randomcoloR::randomColor(count = ll))
    }
  }
  
  y <- 'ident'
  xlab <- 'Expression Level'
  ylab <- 'Identity'
  
  data <- FetchData(object, vars = c(features, split.by, group.by),
                    cells = cells, slot = 'data')
  
  data[,group.by] <- factor( data[,group.by] , levels = levels(idents))
  if(is.null(group.by)){
    data <- cbind(data, object@active.ident)
    group.by <- 'object@active.ident'
  }
  if(is.null(split.by)){
    if(length(features) == 1){
      plot <- ggplot(data, aes_string(y=add.backtick(features), 
                                      x=group.by,
                                      fill = group.by)) + 
        geom_violin() +
        theme_classic() +
        theme(axis.text.x = element_text(angle=90, vjust = 0.5),
              axis.title.x = element_blank())
    }else{
      plot <- NULL
      list_plot <- vector(mode = "list", length = length(features))
      names(list_plot) <- features
      for(feature in features){
        if(feature != features[length(features)]){
          a <- ggplot(data, aes_string(y=add.backtick(feature), 
                                       x=group.by,
                                       fill = group.by))  + 
            geom_violin() +
            theme_classic() +
            theme(axis.text.x = element_blank(),
                  axis.title.x = element_blank())
          list_plot[[feature]] <- a}
        else{
          a <- ggplot(data, aes_string(y=add.backtick(feature), 
                                       x=group.by,
                                       fill = group.by))  + 
            geom_violin() +
            theme_classic() +
            theme(axis.text.x = element_text(angle=90, vjust = 0.5),
                  axis.title.x = element_blank())
          list_plot[[feature]] <- a
          
        }
      }
      plot <- wrap_plots(list_plot, nrow = length(features),
                         guides = 'collect')
    }
    
  }
  #
  # split by!
  #
  if(!is.null(split.by)){
    if(length(split.by) == 1){
      #
      # ONE SPLIT.BY
      #
      colnames(data)[colnames(data) %in% split.by] <- make.unique(rep('split.by', length(split.by)))
      if(length(features) == 1){
        plot <- ggplot(data, aes_string(y=add.backtick(features), 
                                        x=group.by, fill = group.by)) + 
          geom_violin() +
          theme_classic() +
          theme(axis.text.x = element_text(angle=90, vjust = 0.5),
                axis.title.x = element_blank())+
          facet_grid(cols = vars(split.by), scales = 'free')
      }
      if(length(features) > 1){
        plot <- NULL
        list_plot <- vector(mode = "list", length = length(features))
        names(list_plot) <- features
        for(feature in features){
          if(feature == features[1]){
            a <-ggplot(data, aes_string(y=add.backtick(feature), 
                                        x=group.by,
                                        fill = group.by)) + 
              geom_violin() +
              theme_classic() +
              facet_grid(cols = vars(split.by), scales = 'free') +
              theme(axis.text.x = element_blank(),
                    axis.title.x = element_blank())
            
            list_plot[[feature]] <- a
          }
          if(feature != features[length(features)] &
             feature != features[1]){
            a <-ggplot(data, aes_string(y=add.backtick(feature), 
                                        x=group.by,
                                        fill = group.by)) + 
              geom_violin() +
              theme_classic() +
              facet_grid(cols = vars(split.by), scales = 'free')+
              theme(axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    strip.background = element_blank(),
                    strip.text = element_blank()
              )
            
            list_plot[[feature]] <- a
            
          }
          if(feature == features[length(features)]){
            a <-ggplot(data, aes_string(y=add.backtick(feature), 
                                        x=group.by,
                                        fill = group.by)) + 
              geom_violin() +
              theme_classic() +
              facet_grid(cols = vars(split.by), scales = 'free')+
              theme(axis.text.x = element_text(angle=90, 
                                               vjust = 0.5),
                    axis.title.x = element_blank(),
                    strip.background = element_blank(),
                    strip.text = element_blank()
              )
            list_plot[[feature]] <- a
          }
        }
        plot <- wrap_plots(list_plot, nrow = length(features),
                           guides = 'collect')
      }
    }
    if(length(split.by) > 1){
      #
      # MULTIPLE SPLIT.BY
      #
      colnames(data)[colnames(data) %in% split.by] <- make.unique(rep('split.by', length(split.by)))
      if(length(features) == 1){
        plot <- ggplot(data, aes_string(y=add.backtick(features), 
                                        x=group.by, fill = group.by)) + 
          geom_violin() +
          theme_classic() +
          theme(axis.text.x = element_text(angle=90, vjust = 0.5),
                axis.title.x = element_blank())+
          facet_grid(cols = vars(split.by),
                     rows = vars(split.by.1),
                     scales = 'free')
      }
      if(length(features) > 1){
        plot <- NULL
        list_plot <- vector(mode = "list", length = length(features))
        names(list_plot) <- features
        for(feature in features){
          if(feature == features[1]){
            a <-ggplot(data, aes_string(y=add.backtick(feature), 
                                        x=group.by,
                                        fill = group.by)) + 
              geom_violin() +
              theme_classic() +
              facet_grid(cols = vars(split.by),
                         rows = vars(split.by.1),
                         scales = 'free')+
              theme(axis.text.x = element_blank(),
                    axis.title.x = element_blank())
            
            list_plot[[feature]] <- a
          }
          if(feature != features[length(features)] &
             feature != features[1]){
            a <-ggplot(data, aes_string(y=add.backtick(feature), 
                                        x=group.by,
                                        fill = group.by)) + 
              geom_violin() +
              theme_classic() +
              facet_grid(cols = vars(split.by),
                         rows = vars(split.by.1),
                         scales = 'free')+
              theme(axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    strip.background.x = element_blank(),
                    strip.text.x = element_blank()
              )
            
            list_plot[[feature]] <- a
            
          }
          if(feature == features[length(features)]){
            a <-ggplot(data, aes_string(y=add.backtick(feature), 
                                        x=group.by,
                                        fill = group.by)) + 
              geom_violin() +
              theme_classic() +
              facet_grid(cols = vars(split.by),
                         rows = vars(split.by.1),
                         scales = 'free')+
              theme(axis.text.x = element_text(angle=90, 
                                               vjust = 0.5),
                    axis.title.x = element_blank(),
                    strip.background.x = element_blank(),
                    strip.text.x = element_blank()
              )
            list_plot[[feature]] <- a
          }
        }
        plot <- wrap_plots(list_plot, nrow = length(features),
                           guides = 'collect')
      }
    }
    
  }
  
  if(!isTRUE(nucolor)){plot <- plot & scale_fill_manual(values = colors)}
  
  if(group.by == 'object@active.ident'){
    names(colors) <- NULL
    plot <- plot+
      scale_fill_manual(values = colors) + 
      labs(y = 'Idents')+
      NoLegend()
  }
  return(plot)
  
}


ridge_plot <- function(object, features,
                       split.by = NULL, 
                       idents = NULL,
                       colors=NULL,
                       group.by=NULL) {
  require(ggplot2)
  require(Seurat)
  require(patchwork)
  require(formulaic)
  require(gplots)
  require(scales)
  require(ggridges)
  require(cowplot)
  require(randomcoloR)
  
  if (sum(features %in% rownames(object)) != length(features)){
    print(paste('Feature(s)', features[!(features %in% rownames(object))], 'not present in the object'))
    
  }
  if (length(features) > 1){
    stop('Only one gene please!')
  }
  if(!is.null(idents)){
    if (!is.null(x = group.by)) {
      object <- SetIdent(object, value = group.by)
      object <- subset(object, idents = idents)
    }else{
      object <- subset(object, idents = idents)
    }
  }
  
  if (is.null(x = idents)) {
    cells <- colnames(x = object)
  } else {
    cells <- names(x = Idents(object = object)[Idents(object = object) %in% idents])
  }
  
  
  idents <- if (is.null(x = group.by)) {
    Idents(object = object)[cells]
  } else {
    object[[group.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = idents)) {
    idents <- factor(x = idents)
  }
  idents <- factor(x = idents, levels = levels(idents)[order(levels(idents))])
  if (is.null(x = split.by)) {
    split <- NULL
  } else {
    split <- object[[split.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = split)) {
    split <- factor(x = split)
  }
  if (is.null(x = colors)) {
    colors <- hue_pal()(length(x = levels(x = idents)))
    colors <- alpha(colors, alpha = 0.5)
  } 
  if (length(x = colors) < length(x = levels(x = idents))) {
    ll <- length(levels(idents)) - length(colors)
    colors <- c(colors, randomcoloR::randomColor(count = ll))
  }
  if (length(x = colors) < length(x = levels(x = split))) {
    ll <- length(levels(split)) - length(colors)
    colors <- c(colors, randomcoloR::randomColor(count = ll))
  }
  if(length(split)>0){
    colors <- rep_len(x = colors, length.out = length(x = levels(x = split)))
    names(x = colors) <- levels(x = split)
  }
  if(length(split)==0){
    names(x = colors) <- levels(idents)
  }
  
  y <- 'ident'
  xlab <- 'Expression Level'
  ylab <- 'Identity'
  
  if(sum(features %in% rownames(object)) == length(features)){
    data <- FetchData(object, vars = c(features, split.by, group.by),
                      cells = cells, slot = 'data')
    data[,group.by] <- factor( data[,group.by] , levels = levels(idents))
    if(is.null(group.by)){
      data <- cbind(data, object@active.ident)
      group.by <- 'object@active.ident'
    }
    if(is.null(split.by)){
      plot <- ggplot(data, aes_string(x=features, 
                                      y=group.by,
                                      fill = group.by)) +
        geom_density_ridges()+
        theme_cowplot() +
        labs(title = features)+ 
        scale_fill_manual(values = colors, labels = names(colors))
    }
    if(!is.null(split.by)){
      plot <- ggplot(data, aes_string(x=features, 
                                      y=group.by, fill = split.by)) +
        geom_density_ridges()+
        theme_cowplot() +
        labs(title = features)  + 
        scale_fill_manual(values = colors, labels = names(colors))
    }
  }
  if(group.by == 'object@active.ident'){
    names(colors) <- NULL
    plot <- plot+
      scale_fill_manual(values = colors) + 
      labs(y = 'Idents')
    if(is.null(split.by)){
      plot <- plot + NoLegend()
    }
  }
  
  return(plot)
  
}


#
# Barplot figure
#

barplot_figure <- function(object, features, group.by,
                           idents= NULL, cols = NULL){
  require(Seurat)
  require(gridExtra)
  require(ggpubr)
  require(formulaic)
  
  if(is.null(cols)){nucolor <- T}else{nucolor <- F}
  
  if (sum(features %in% rownames(object)) != length(features)){
    print(paste('Feature(s)', features[!(features %in% rownames(object))], 'not present in the object'))
    features <- features[features %in% rownames(object)]
  }
 
  #
  # select cells
  #
  if (is.null(x = idents)) {
    cells <- colnames(x = object)
  } else {
    cells <- names(x = Idents(object = object)[Idents(object = object) %in% idents])
  }
  
  idents <- if (is.null(x = group.by)) {
    Idents(object = object)[cells]
  } else {
    object[[group.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = idents)) {
    idents <- droplevels(factor(x = idents))
  }
  
  idents <- droplevels(factor(x = idents, levels = levels(idents)[order(levels(idents))]))
  # print(idents)
  
  data <- FetchData(object, vars = c(features, group))
  colnames(data)[ncol(data)] <- 'class'
  data <- data[order(data$class),]
  data$cells <- rownames(data)
  data$cells <- factor(data[,ncol(data)], levels =data[,ncol(data)])
  head(data)
  plot_list <- list()
  if(is.null(cols)){
    for(i in 1:length(features)){
      p = ggplot(data, aes_string(x='cells',
                                  y=add.backtick(features[i], 
                                                 include.backtick = 'all'), 
                                  fill = 'class', color = 'class')) + 
        geom_bar(stat = "identity") +
        labs(x = '')+
        theme_classic() +
        scale_y_continuous(
          labels = scales::number_format(accuracy = 0.1))+
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = 'none') 
      plot_list[[i]] = p
    }
    j <- ggplot(data, aes_string(x='cells', 
                                 y=add.backtick(features[i], 
                                                include.backtick = 'all'), 
                                 fill = 'class', color = 'class')) +
      geom_bar(stat = "identity") +
      theme_classic() +
      theme(legend.position = 'bottom') 
    p = get_legend(j)
    
    plot_list[[i+1]] = p
  }else{
    for(i in 1:length(features)){
      p = ggplot(data, aes_string(x='cells',
                                  y=add.backtick(features[i],
                                                 include.backtick = 'all'),
                                  fill = 'class', color = 'class')) + 
        geom_bar(stat = "identity") +
        labs(x = '')+
        theme_classic() +
        scale_color_manual(values = cols)+
        scale_y_continuous(
          labels = scales::number_format(accuracy = 0.1))+
        scale_fill_manual(values = cols)+
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = 'none') 
      plot_list[[i]] = p
    }
    j <- ggplot(data, aes_string(x='cells',
                                 y=add.backtick(features[i],
                                                include.backtick = 'all'), 
                                 fill = 'class', color = 'class')) +
      geom_bar(stat = "identity") +
      theme_classic() +
      scale_color_manual(values = cols)+
      scale_fill_manual(values = cols)+
      theme(legend.position = 'bottom') 
    p = get_legend(j)
    
    plot_list[[i+1]] = p
  }
  
  return(do.call("ggarrange", c(plot_list, ncol=1)))
}




