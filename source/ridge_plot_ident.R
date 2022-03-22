ridge_plot_one_ident <- function(object, features,
                                 split.by = NULL, 
                                 idents = NULL,
                                 colors=NULL,
                                 group.by=NULL) {
  require(ggplot2)
  require(Seurat)
  require(patchwork)
  require(formulaic)
  require(gplots)
  library(scales)
  library(ggridges)
  require(cowplot)
  
  # 
  # checkpoints
  #
  if (sum(features %in% rownames(object)) != length(features)){
    print(paste('Feature(s)', features[!(features %in% rownames(object))], 'not present in the object'))
    features <- features[features %in% rownames(object)]
  }
  if (length(idents) > 1){
    stop('ONLY ONE IDENT WITH THIS FUNCTION')
  }
if (!is.null(x = group.by)) {
  object <- SetIdent(object, value = group.by)
  if(idents %in% levels(object)){
    object <- subset(object, idents = idents)
  }else{
    stop('IDENT NOT PRESENT IN THE GROUP.BY')
  }
  
  
}

#
# select cells!
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
if (is.null(x = colors)) {
  colors <- scales::hue_pal()(length(x = levels(x = idents)))
  colors <- alpha(colors, alpha = 0.5)
} else {
  colors <- Col2Hex(colors)
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
    plot <- ggplot(data, aes_string(x=features, 
                                    y=group.by,
                                    fill = group.by)) +
      geom_density_ridges()+
      theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank())+
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      coord_cartesian(clip = "off") + 
      labs(title = features)
    
  }else{
    plot <- NULL
    list_plot <- vector(mode = "list", length = length(features))
    names(list_plot) <- features
    for(feature in features){
      if(feature != features[length(features)]){
        a <- ggplot(data, aes_string(x=add.backtick(feature), 
                                     y=group.by,
                                     fill = group.by))  + 
          geom_density_ridges() +
          labs(title = feature) +
          theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
          theme(axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                strip.background = element_rect(fill = '#4287f500')) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_discrete(expand = c(0, 0)) +
          coord_cartesian(clip = "off")
        list_plot[[feature]] <- a}
      else{
        a <- ggplot(data, aes_string(x=add.backtick(feature), 
                                     y=group.by,
                                     fill = group.by))  + 
          geom_density_ridges() +
          labs(title = feature)+
          theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
          theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank())+
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_discrete(expand = c(0, 0)) +
          coord_cartesian(clip = "off")
        list_plot[[feature]] <- a
        
      }
    }
    plot <- wrap_plots(list_plot, ncol = length(features),
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
      plot <- ggplot(data, aes_string(x=add.backtick(features), 
                                      y=group.by, fill = group.by)) + 
        geom_density_ridges() +
        labs(title = features)+
        theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())+
        facet_grid(cols = vars(split.by), scales = 'free')
    }
    if(length(features) > 1){
      plot <- NULL
      list_plot <- vector(mode = "list", length = length(features))
      names(list_plot) <- features
      for(feature in features){
        if(feature == features[1]){
          a <-ggplot(data, aes_string(x=add.backtick(feature), 
                                      y=group.by,
                                      fill = group.by)) + 
            geom_density_ridges() +
            labs(title = feature)+
            theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
            facet_grid(cols = vars(split.by), scales = 'free') +
            theme(axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  strip.background = element_rect(fill = '#4287f500'))
          
          list_plot[[feature]] <- a
        }
        if(feature != features[length(features)] &
           feature != features[1]){
          a <-ggplot(data, aes_string(x=add.backtick(feature), 
                                      y=group.by,
                                      fill = group.by)) + 
            geom_density_ridges() +
            labs(title = feature)+
            theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
            facet_grid(cols = vars(split.by), scales = 'free')+
            theme(axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  strip.background = element_blank(),
                  strip.text = element_blank()
            )
          
          list_plot[[feature]] <- a
          
        }
        if(feature == features[length(features)]){
          a <-ggplot(data, aes_string(x=add.backtick(feature), 
                                      y=group.by,
                                      fill = group.by)) + 
            geom_density_ridges() +
            labs(title = feature)+
            theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
            facet_grid(cols = vars(split.by), scales = 'free')+
            theme(axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
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
      plot <- ggplot(data, aes_string(x=add.backtick(features), 
                                      y=group.by, fill = group.by)) + 
        geom_density_ridges() +
        labs(title = feature)+
        theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())+
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
          a <-ggplot(data, aes_string(x=add.backtick(feature), 
                                      y=group.by,
                                      fill = group.by)) + 
            geom_density_ridges() +
            labs(title = feature)+
            theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
            facet_grid(cols = vars(split.by),
                       rows = vars(split.by.1),
                       scales = 'free')+
            theme(axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  strip.background = element_rect(fill = '#4287f500'))
          
          list_plot[[feature]] <- a
        }
        if(feature != features[length(features)] &
           feature != features[1]){
          a <-ggplot(data, aes_string(x=add.backtick(feature), 
                                      y=group.by,
                                      fill = group.by)) + 
            geom_density_ridges() +
            labs(title = feature)+
            theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
            facet_grid(cols = vars(split.by),
                       rows = vars(split.by.1),
                       scales = 'free')+
            theme(axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  strip.background.x = element_blank(),
                  strip.text.x = element_blank()
            )
          
          list_plot[[feature]] <- a
          
        }
        if(feature == features[length(features)]){
          a <-ggplot(data, aes_string(x=add.backtick(feature), 
                                      y=group.by,
                                      fill = group.by)) + 
            geom_density_ridges() +
            labs(title = feature)+
            theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
            facet_grid(cols = vars(split.by),
                       rows = vars(split.by.1),
                       scales = 'free')+
            theme(axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
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

if(group.by == 'object@active.ident'){
  names(colors) <- NULL
  plot <- plot+
    scale_fill_manual(values = colors) + 
    labs(y = 'Idents')+
    NoLegend()
}

plot <- plot &
  scale_x_continuous(expand = c(0, 0)) &
  scale_y_discrete(expand = c(0, 0)) &
  coord_cartesian(clip = "off") &
  theme(strip.background.y = element_rect(fill = '#4287f500'),
        strip.background.x = element_rect(fill = '#4287f500'))

return(plot)

}