hist_plot <- function(obs,
                      breaks_hist = seq(from = 0, to = 40, by = 5),
                      colbin,
                      double_bounded = 1,
                      lim_x = c(0, 40),
                      break_points = c(0, 20, 40),
                      empty = 0,
                      posterior = 0,
                      post_data,
                      alpha_point = 0.3,
                      alpha_line = 0.3,
                      alpha_hist = 1,
                      col_post = "#a19987") {
  plot <- ggplot(obs, aes(x = x)) +
    
    theme_bw() +
    theme(
      title = element_text(size = 7.5),
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white", color = "white"),
      plot.margin = margin(0.1, 5, 0.1, 0.1, "cm")
    ) +
    annotate(
      xmin = -0.1 * (lim_x[2] - lim_x[1]),
      xmax = 0,
      ymin = -Inf,
      ymax = Inf,
      geom = "rect",
      fill = "#f5f2ed"
    ) +
    geom_vline(
      xintercept = 0,
      color = "#b8b5ab",
      linewidth = 0.3,
      linetype = "dotted"
    )
  
  if (double_bounded == 1){
    plot <- plot + 
      annotate(
        xmin = lim_x[2],
        xmax = lim_x[2] + 0.1 * (lim_x[2] - lim_x[1]),
        ymin = -Inf,
        ymax = Inf,
        geom = "rect",
        fill = "#f5f2ed"
      ) +
      geom_vline(
        xintercept = lim_x[2],
        color = "#b8b5ab",
        linewidth = 0.3,
        linetype = "dotted"
      )
  }
  
  plot <- plot + 
    scale_y_continuous(expand = expansion(mult = 0.1)) +
    
    scale_x_continuous(
      limits = c(-0.1 * (lim_x[2] - lim_x[1]),
                 lim_x[2] + 0.1 * (lim_x[2] - lim_x[1])),
      expand = c(0, 0),
      breaks = break_points
    ) +
    
    geom_histogram(
      breaks = breaks_hist,
      color = ifelse(empty == 0, "white", NA),
      linewidth = 0.3,
      alpha = alpha_hist,
      fill = ifelse(empty == 0, colbin, NA)
    ) +
    geom_hline(
      yintercept = 0,
      color = "#b8b5ab",
      linewidth = 0.3,
      linetype = "dotted"
    )
  
  if (posterior == 1) {
    plot <- plot +
      geom_point(
        data = post_data,
        aes(bin, y = count, group = sample),
        pch = 21,
        colour = "white",
        fill = col_post,
        size = 2,
        alpha = alpha_point
      ) +
      geom_line(
        data = post_data,
        aes(bin, y = count, group = sample),
        col = col_post,
        size = 0.5,
        alpha = alpha_line
      )
  }
  
  plot <- plot + 
    geom_histogram(
      breaks = breaks_hist,
      color = ifelse(empty == 0, "#3d3d3d", NA),
      linewidth = 0.5,
      fill = ifelse(empty == 0, "transparent", NA)
    )
  
  return(plot)
}

hist_plot %>% saveRDS("./functions/hist_plot.rds")

