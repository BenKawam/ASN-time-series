post_plot <- function(post_data,
                      target_data,
                      target_values = TRUE,
                      min = 0,
                      max = 2,
                      scale_slab = 0.7,
                      bounds = "none",
                      alpha_filling = 0.05,
                      alpha_outline = 0.15,
                      line_width = 0.1) {
  plot <- post_data %>%
    ggplot(aes(x = value, y = param)) +
    # Theme and layout
    theme_bw() +
    labs(x = "", y = "", color = 'Percentile Interval') +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 6.5),
      axis.text.x = element_text(vjust = -1.5),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(
        color = "#c2c2b0",
        size = 0.3,
        linetype = "dotted"
      ),
      strip.background = element_rect(fill = "white", color = "white")
    ) +
    # filling of the slab
    stat_slab(
      aes(group = iter),
      slab_color = NA,
      fill = "#a19987",
      scale = scale_slab,
      normalize = "groups",
      slab_alpha = alpha_filling,
      trim = FALSE
    ) +
    
    # Contour of slab
    stat_slab(
      aes(group = iter),
      slab_color = "#363533",
      fill = NA,
      slab_linewidth = line_width,
      scale = scale_slab,
      normalize = "groups",
      slab_alpha = alpha_outline,
      trim = FALSE
    ) +
    scale_x_continuous(limits = c(min - 0.1 * (max - min), max + 0.1 * (max - min)),
                       expand = c(0, 0))
  
  if (target_values == TRUE) {
    # Target values
    plot <- plot + geom_point(
      data = target_data,
      pch = 21,
      colour = "white",
      fill = "#ab0f24",
      size = 1.8
    ) +
      geom_text(
        data = target_data,
        aes(label = round(value, digits = 2)),
        nudge_y = -0.21,
        size = 2.6,
        col = "#ab0f24"
      )}
  
  if (bounds == "lower" | bounds == "double") {
    plot <- plot +
      geom_vline(
        xintercept = min,
        color = "#b8b5ab",
        linewidth = 0.3,
        linetype = "dotted"
      ) +
      # Lower bound
      annotate(
        xmin = -0.1 * (max - min),
        xmax = 0,
        ymin = -Inf,
        ymax = Inf,
        geom = "rect",
        fill = "#f5f2ed"
      )
  }
  
  if (bounds == "double") {
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
  
  return(plot)
}

post_plot %>% saveRDS("./functions/post_plot.rds")
