plot_interval <- function(
    draws,
    true_values,
    naive_values,
    double_bounded,
    max,
    min = 0,
    true_targets = 1,
    naive_targets = 1
){
  
  plot <- ggplot(draws, aes(x = value, y = param)) +
    # theme
    scale_color_manual(values = c("#dedace", "#ccc7ba", "#b3aea1", "#a39d8e")) +
    theme_bw() +
    theme(
      title = element_text(size = 7.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = margin(1, 1, 1, 1),
      axis.text.y = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(vjust = -0.2, size = 8),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white", color = "white")
    ) +
    
    # Lower bound
    annotate(
      xmin = -0.1 * (max - min),
      xmax = 0,
      ymin = -Inf,
      ymax = Inf,
      geom = "rect",
      fill = "#f5f2ed"
    ) +
    scale_x_continuous(
      limits = c(-0.1 * (max - min), max + 0.1 * (max - min)),
      expand = c(0, 0),
      breaks = c(0, max / 2, max)
    ) +
    geom_vline(
      xintercept = min,
      color = "#b8b5ab",
      linewidth = 0.3,
      linetype = "dotted"
    )
  
  
  if (double_bounded == 1) {
    plot <- plot +
      # upper bound
      annotate(
        xmin = max,
        xmax = max + 0.1 * (max - min),
        ymin = -Inf,
        ymax = Inf,
        geom = "rect",
        fill = "#f5f2ed"
      ) +
      geom_vline(
        xintercept = max,
        color = "#b8b5ab",
        linewidth = 0.3,
        linetype = "dotted"
      )
  }
  
  
  plot <- plot + stat_interval(.width = c(.25, .5, .75, .95),
                               size = 2.5,
                               alpha = 0.6)
  
  if(true_targets == 1){
    plot <- plot + geom_point(
      data = true_values,
      pch = 21,
      colour = "white",
      fill = "#ab0f24",
      size = 1.8
    )
  }
  
  if(naive_targets == 1){
    plot <- plot + geom_point(
      data = naive_values,
      shape = 1,
      colour = "#8c2d3a",
      size = 1.8,
      alpha = 0.5
    )
  }
  return(plot) 
}

plot_interval %>% saveRDS("./functions/plot_interval.rds")

