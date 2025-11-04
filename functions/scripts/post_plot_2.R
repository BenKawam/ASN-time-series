post_plot_2 <- function(
    posterior_draws,
    min = -10,
    max = 10,
    y_space = -0.175,
    scale_dens = 1,
    target_values = TRUE,
    target_data,
    n_breaks = 4,
    alpha_filling = 0.1,
    alpha_outline = 0.15,
    line_width = 0.1,
    lower_bound = 0,
    upper_bound = 0) {
  
  plot <- posterior_draws %>%
    ggplot(aes(x = value)) +
    # Theme and layout
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.text.x = element_text(vjust = -0.5, size = 7),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_line(
        color = "#dedbd3",
        linewidth = 0.3,
        linetype = "dotted"),
      strip.background = element_rect(fill = "white", color = "white"),
      plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")
    ) +
    geom_vline(
      xintercept = 0,
      color = "#b8b5ab",
      linewidth = 0.3,
      linetype = "dotted"
    ) +
    geom_hline(
      yintercept = 0,
      color = "#dedbd3",
      linewidth = 0.3,
      linetype = "dotted"
    ) +
    
    scale_x_continuous(limits = c(min - 0.1 * (max - min), max + 0.1 * (max - min)),
                       expand = c(0, 0),
                       breaks = seq(from = min, to = max, by = (max - min) / n_breaks)) +
    
    scale_y_continuous(limits = c(-0.35, 1.15),
                       expand = c(0, 0))  +
    # filling of the slab
    stat_slab(
      aes(group = iter),
      slab_color = NA,
      fill = "#a19987",
      scale = scale_dens,
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
      scale = scale_dens,
      normalize = "groups",
      slab_alpha = alpha_outline,
      trim = FALSE
    )
  
  if (target_values == TRUE) {
    # Target values
    plot <- plot + geom_point(
      data = tibble(value = target_data, iter = 1, y = 0),
      aes(y = y),
      pch = 21,
      colour = "white",
      fill = "#ab0f24",
      size = 1.8
    ) +
      geom_text(
        data = tibble(value = target_data, iter = 1, y = 0),
        aes(label = round(value, digits = 2), y = y),
        nudge_y = y_space,
        size = 1.8,
        col = "#7d0e1d"
      )}
  
  if (lower_bound == 1) {
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
  
  if (upper_bound == 1) {
    plot <- plot +
      geom_vline(
        xintercept = max,
        color = "#b8b5ab",
        linewidth = 0.3,
        linetype = "dotted"
      ) +
      # Lower bound
      annotate(
        xmin = max + 0.1 * (max - min),
        xmax = max,
        ymin = -Inf,
        ymax = Inf,
        geom = "rect",
        fill = "#f5f2ed"
      ) +
      geom_vline(
        xintercept = 1,
        color = "#b8b5ab",
        linewidth = 0.3,
        linetype = "dotted"
      )
  }
  
  
  return(plot)
}

post_plot_2 %>% saveRDS("./functions/post_plot_2.rds")

