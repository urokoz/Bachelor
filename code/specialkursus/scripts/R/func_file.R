library(tidyverse)
library(ggpubr)

#function for sigmoid transformation
sigmoid <- function(si) {
  if (si >= 3.5) {
    si = 1 / (1 + exp(-(si-3.5)*0.5))
  } 
  else {
    si = 1 / (1 + exp(-(si-3.5)*0.2))
  }
}


#function for scatterplots ----------------------------------------
scatter_plot <- function(df, variable1, variable2, si_filter = TRUE, title, xlab, ylab){
  
  if (si_filter == TRUE) {
    df <- df %>% 
      filter(si >=2)
    plt <- ggplot(df, mapping = aes(x = {{variable1}}, 
                                    y = {{variable2}})) +
      
      geom_point(size = 1) +
      labs(title = title, 
           x = xlab, 
           y = ylab)
      
    plt <- plt + stat_cor(method = "spearman")
  }
  
  else {
    plt <- ggplot(df, mapping = aes(x = {{variable1}}, 
                                    y = {{variable2}})) +
      
      geom_point(size = 1) +
      labs(title = title, 
           x = xlab, 
           y = ylab)
    
    plt <- plt + stat_cor(method = "spearman")
  }
  
  return(plt)
}

#function for boxplots ----------------------------------------
boxplot <- function(df, variable_dis, variable_con, title, xlab, ylab) {
  
  plt <- ggplot(df, mapping = aes(y = {{variable_dis}},
                                  x = {{variable_con}}, 
                                  fill = {{variable_dis}})) + 
    geom_boxplot() + 
    labs(title = title, 
         x = xlab,
         y = ylab) +
    theme_classic()
  
  return(plt)
}



