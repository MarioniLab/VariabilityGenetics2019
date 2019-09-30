theme_mike <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
  + theme(plot.title = element_text(face = "bold",
                                    size = rel(1.2), hjust = 0.5),
          text = element_text(size=16),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold",size = rel(1)),
          axis.title.y = element_text(angle=90,vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(), 
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          panel.grid.major = element_line(colour="#f0f0f0"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.width= unit(1.0, "cm"),
          legend.key.height = unit(0.5, "cm"),
          legend.title = element_text(face="italic"),
          plot.margin=unit(c(10,5,5,5),"mm"),
          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold")
  ))
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",
                 manual_pal(values = c("#386cb0", "#fdb462", "#ef3b2c",
                                       "#fb9a99","#984ea3", "#ffff33",
                                       "#662506", "#b2df8a", "#a6cee3", 
                                       "#CA834E", "#518A87", "#5B113C", 
                                       "#55813B", "#E704C4", "#DDAEA2", 
                                       "#77837F", "#A53327", "#608EFF", 
                                       "#B599D7", "#A50149", "#4E0025", 
                                       "#C9B1A9","#7fc97f", "#3E89BE",
                                       "#33a02c", "#fb9a99", "#e31a1c",
                                       "#fdbf6f", "#ff7f00", "#cab2d6",
                                       "#6a3d9a", "#ffff99", "#b15928",
                                       "#4A1930", "#E8C282", "#E7DBBC", 
                                       "#A68486")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour", "Publication",
                 manual_pal(values = c("#386cb0", "#fdb462",
                                       "#7fc97f","#ef3b2c","#662506",
                                       "#a6cee3","#fb9a99","#984ea3",
                                       "#ffff33", "#b2df8a", "#3E89BE", 
                                       "#CA834E", "#518A87", "#5B113C", 
                                       "#55813B", "#E704C4", "#DDAEA2", 
                                       "#77837F", "#A53327", "#608EFF", 
                                       "#B599D7", "#A50149", "#4E0025", 
                                       "#C9B1A9",
                                       "#33a02c", "#fb9a99", "#e31a1c",
                                       "#fdbf6f", "#ff7f00", "#cab2d6",
                                       "#6a3d9a", "#ffff99", "#b15928", 
                                       "#4A1930", "#E8C282", "#E7DBBC",
                                       "#A68486")), ...)
  
}

# This is from David Robinson on GitHub: https://github.com/dgrtwo/drlib/blob/master/R/reorder_within.R

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}


scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}


scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}