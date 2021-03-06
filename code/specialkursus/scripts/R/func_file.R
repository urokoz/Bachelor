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
  return(si)
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
           y = ylab) +
      ylim(0, 50)
      
    plt <- plt + stat_cor(method = "spearman", 
                          digits = 3,
                          size = 4.5,
                          hjust = -0.9,
                          vjust = 1)
  }
  
  else {
    plt <- ggplot(df, mapping = aes(x = {{variable1}}, 
                                    y = {{variable2}})) +
      
      geom_point(size = 1, alpha = 0.2) + 
      labs(title = title, 
           x = xlab, 
           y = ylab) +
      #ylim(0, 50) + 
      theme(plot.title = element_text(size=15))
    
    plt <- plt + stat_cor(method = "spearman", 
                          digits = 4,
                          size = 3,
                          hjust = -1.4,
                          vjust = 1)
  }
  
  return(plt)
}

#function for boxplots ----------------------------------------
boxplot <- function(df, variable_dis, variable_con, title, xlab, ylab) {
  
  plt <- ggplot(df, mapping = aes(y = {{variable_con}},
                                  x = {{variable_dis}}, 
                                  fill = {{variable_dis}})) + 
    geom_violin() + 
    geom_boxplot(width = 0.05,
                 color = "black",
                 fill = "white",
                 outlier.shape = NA) +
    geom_signif(comparisons = list(c("Non Reactive", "Reactive"))) +
    labs(title = title, 
         x = xlab,
         y = ylab) +
    theme_classic() +
    theme(legend.position="none")

  return(plt)
}

boxplot_core_sim <- function(df, variable_dis, variable_con, title, xlab, ylab) {
  max_val <- df %>% 
    summarise(max = max({{variable_con}})) %>% 
    unlist()
  plt <- ggplot(df, mapping = aes(y = {{variable_con}},
                                  x = {{variable_dis}}, 
                                  fill = {{variable_dis}})) + 
    geom_violin() + 
    geom_boxplot(width = 0.05,
                 color = "black",
                 fill = "white",
                 outlier.shape = NA) +
    geom_signif(comparisons = list(c("Cross reactive", 
                                     "Non cross reactive"))) +
    labs(title = title, 
         x = xlab,
         y = ylab) +
    theme_classic() +
    theme(legend.position="none")
  
  return(plt)
}


