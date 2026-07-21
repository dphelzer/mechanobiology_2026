library(dplyr)
library(purrr)
library(spatstat)
library(ggplot2)
library(ggrastr)
library(patchwork)

results_dir <- "patch_analysis_results"

tissue_windows <- readRDS(file.path(results_dir, "tissue_windows.rds"))
coords <- readRDS(file.path(results_dir, "coords.rds"))
sample_ids <- names(tissue_windows)

base_theme <- theme_minimal(base_size = 8) +
  theme(
    text = element_text(color = "black", size = 8),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )

class_colors <- c(
  "Negative" = "blue",
  "Positive" = "red"
)

boundary_plot_df <- imap_dfr(tissue_windows, function(win, s) {
  if (is.null(win)) return(NULL)
  imap_dfr(win$bdry, function(seg, i) {
    data.frame(x = seg$x, y = seg$y, segment = i, Sample = s)
  })
})


# Bin rectangles, filled by positive-cell count -----------------------

build_bin_rects <- function(sample_id, bin_size) {
  
  win <- tissue_windows[[sample_id]]
  
  df_s <- coords %>%
    filter(Sample == sample_id, Class == "Positive")
  
  nx <- ceiling(diff(win$xrange) / bin_size)
  ny <- ceiling(diff(win$yrange) / bin_size)
  
  xbreaks <- seq(win$xrange[1], win$xrange[2], length.out = nx + 1)
  ybreaks <- seq(win$yrange[1], win$yrange[2], length.out = ny + 1)
  
  xi <- cut(df_s$CentroidX_um, breaks = xbreaks, include.lowest = TRUE, labels = FALSE)
  yi <- cut(df_s$CentroidY_um, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
  
  counts <- as.vector(table(
    factor(xi, levels = 1:nx),
    factor(yi, levels = 1:ny)
  ))
  
  xcenters <- (xbreaks[-length(xbreaks)] + xbreaks[-1]) / 2
  ycenters <- (ybreaks[-length(ybreaks)] + ybreaks[-1]) / 2
  
  bin_grid <- expand.grid(
    xcenter = xcenters,
    ycenter = ycenters,
    KEEP.OUT.ATTRS = FALSE
  ) %>%
    mutate(
      count = counts,
      xmin = xcenter - bin_size / 2,
      xmax = xcenter + bin_size / 2,
      ymin = ycenter - bin_size / 2,
      ymax = ycenter + bin_size / 2,
      inside = inside.owin(xcenter, ycenter, win)
    ) %>%
    filter(inside)
  
  bin_grid
}

plot_bin_grid_sample <- function(sample_id, bin_size) {
  
  bin_grid <- build_bin_rects(sample_id, bin_size)
  
  ggplot() +
    geom_polygon(
      data = boundary_plot_df %>% filter(Sample == sample_id),
      aes(x = x, y = y, group = interaction(Sample, segment)),
      fill = NA,
      color = "black",
      linewidth = 0.3
    ) +
    rasterise(
      geom_point(
        data = coords %>% filter(Sample == sample_id),
        aes(x = CentroidX_um, y = CentroidY_um, color = Class),
        size = 0.01,
        alpha = 0.3
      ),
      dpi = 600
    ) +
    geom_rect(
      data = bin_grid,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = count),
      color = "black",
      linewidth = 0.1,
      alpha = 0.5
    ) +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "Reds")) +
    scale_color_manual(values = class_colors) +
    coord_equal() +
    labs(title = paste(sample_id, "| bin =", bin_size, "um")) +
    base_theme +
    theme(legend.position = "none")
}

sample_ids_to_plot <- sample_ids
bin_size_to_check <- 75

bin_grid_plots <- map(sample_ids_to_plot, plot_bin_grid_sample, bin_size = bin_size_to_check)

wrap_plots(bin_grid_plots, ncol = 4)

ggsave("Figures/grid_overlay.svg",
       height = 15, width = 15, units = "cm")
