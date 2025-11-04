binom_pmf <- function(x_txt = 1,
                      data,
                      col_data = "#363533",
                      col_post = "#363533",
                      empty = 0,
                      ln_width = 0.4,
                      otl_col = "white",
                      posterior = 0,
                      post_data,
                      post_transparent = 0) {
  plot <- ggplot(data, aes(x = s_to, y = count)) +
    theme_bw( ) +
    theme(
      title = element_text(size = 7.5),
      legend.position = "none",
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white", color = "white"),
      plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
    ) +
    scale_y_continuous(expand = expansion(mult = 0.1)) +
    geom_hline(
      yintercept = 0,
      color = "#b8b5ab",
      linewidth = 0.3,
      linetype = "dotted"
    )
  
  if (empty == 0) {
    plot <- plot + geom_bar(
      stat = "identity",
      fill = col_data,
      width = 0.4,
      linewidth = ln_width,
      col = otl_col
    )
  }
  
  if (empty == 1) {
    plot <- plot + geom_bar(
      stat = "identity",
      fill = NA,
      width = 0.4,
      col = NA
    )
  }
  
  if (x_txt == 0) {
    plot <- plot +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  
  if (posterior == 1) {
    plot <- plot +
      geom_point(
        data = post_data,
        aes(group = as.factor(sample), fill = s_to),
        pch = 21,
        colour = "white",
        size = 2,
        alpha = ifelse(post_transparent == 0, 0.3, 0)
      ) +
      scale_fill_manual(values = col_post) +
      geom_line(
        data = post_data,
        aes(group = as.factor(sample)),
        col = "#a1a1a1",
        linewidth = 0.5,
        alpha = ifelse(post_transparent == 0, 0.15, 0)
      ) +
      geom_bar(
        stat = "identity",
        linewidth = 0.5,
        width = 0.4,
        fill = "transparent",
        color = "#3d3d3d"
      )
  }
  
  return(plot)
}

binom_pmf %>% saveRDS("./functions/binom_pmf.rds")

