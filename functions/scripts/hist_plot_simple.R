hist_plot_simple <- function(vec, # Vector of numbers
                             fill = "#27495c",
                             col = "#f7f7f0",
                             min = 0,
                             min_count = 5,
                             breaks = 20,
                             max,
                             linewidth = 0.5) {
  # Reformatting
  vec %>%
    tibble() %>%
    ggplot(aes(., after_stat(ncount)))  +
    # Histogram
    geom_histogram(
      color = col,
      linewidth = linewidth,
      fill = fill,
      breaks = seq(from = min, to = max, by = ((max - min) / breaks))
    ) +
    stat_bin(
      breaks = seq(from = min, to = max, by = ((max - min) / breaks)),
      geom = "text",
      aes(label = after_stat(
        if_else (condition = count > 0 & count < min_count,
                 as.character(count), "")
      )),
      vjust = -1,
      size = 2,
      col = "#4d4d4d"
    ) +
    # Layour
    theme_bw() +
    theme(
      axis.text.x = element_text(vjust = 0),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(
        color = "#ccccc0",
        linewidth = 0.3,
        linetype = "dotted"
      ),
      panel.border = element_blank(),
      strip.background = element_rect(fill = "white", color = "white"),
      plot.margin = margin(0, 2, 0, 0, "cm"),
    ) +
    labs(y = "", x = "")
}

hist_plot_simple %>% saveRDS("./functions/hist_plot_simple.rds")
