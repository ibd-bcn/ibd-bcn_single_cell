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
  require(reshape2)
  
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
  
  if(is.null(x = group.by)){stop('group.by IS MANDATORY')}
  
if (!is.null(x = group.by)) {
  object <- SetIdent(object, value = group.by)
  if(idents %in% levels(object)){
    object <- subset(object, idents = idents)
  }else{
    stop('idents NOT PRESENT IN group.by')
  }}

#
# select cells
#
  cells <- names(x = Idents(object = object)[Idents(object = object) %in% idents])
  idents <-object[[group.by, drop = TRUE]][cells]
  
  if (!is.factor(x = idents)) {
    idents <- factor(x = idents)
  }else{
    idents <- droplevels(idents)
  }
  
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

  data[,group.by] <- factor(data[,group.by] , levels = levels(idents))
  data1 <- melt(data)
  head(data1)

if(is.null(split.by)){
  plot <- ggplot(data1, aes_string(x='value',
                                     y='variable',
                                     fill = group.by)) +
      geom_density_ridges()+
      theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank())+
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      coord_cartesian(clip = "off") +
      labs(title = levels(idents))
  }

  #
  # split by!
  #
# 
  if(!is.null(split.by)){
    colnames(data1)[colnames(data1) %in% split.by] <- make.unique(rep('split.by', length(split.by)))
    print(head(data1))
    #
    # ONE SPLIT.BY
    #
    if(length(split.by) == 1){
      plot <- ggplot(data1, aes(x=value,
                               y=variable,
                               fill = split.by)) +
          geom_density_ridges() +
          labs(title = levels(idents))+
          theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
          theme(axis.title.x = element_blank(),
                axis.title.y = element_blank())
    }
    #
    # MULTIPLE SPLIT.BY
    #
    if(length(split.by) > 1){
      plot <- ggplot(data1, aes_string(x='value',
                                      y='variable',
                                      fill = 'split.by')) +
              geom_density_ridges() +
              labs(title = levels(idents))+
              theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_blank())+
              facet_grid(cols = vars(split.by.1),
                         scales = 'free')
        }

      }

plot <- plot +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(strip.background.y = element_rect(fill = '#4287f500'),
        strip.background.x = element_rect(fill = '#4287f500'))

return(plot)

}

