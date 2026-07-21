library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(spatstat)
library(ggrastr)
library(patchwork)
library(RColorBrewer)


# Setup: load boundaries, coordinates, tissue windows --------------------

## Load boundary files per sample -----
tissue_windows <- readRDS("patch_analysis_results/tissue_windows.rds")
coords <- readRDS("patch_analysis_results/coords.rds")
sample_ids <- names(tissue_windows)

## Shared plot styling -----
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

boundary_plot_df <- imap_dfr(tissue_windows, function(win, s) {
  if (is.null(win)) return(NULL)
  imap_dfr(win$bdry, function(seg, i) {
    data.frame(x = seg$x, y = seg$y, segment = i, Sample = s)
  })
})


# Density visualization: per-sample bandwidth (diagnostic use) ----------

## NOTE: this version uses each sample's own bw.diggle() bandwidth, so
## values are NOT directly comparable across samples. Useful for
## inspecting one sample in isolation; use the normalized (shared
## bandwidth) contour version below for any cross-sample comparison.

## Build ppp object for CD68+ cells, one sample -----
build_positive_ppp <- function(sample_id) {
  
  win <- tissue_windows[[sample_id]]
  
  df_s <- coords %>%
    filter(
      Sample == sample_id,
      Class == "Positive"
    )
  
  ppp(
    x = df_s$CentroidX_um,
    y = df_s$CentroidY_um,
    window = win
  )
}

## Compute density surface for one sample -----
compute_density_df <- function(sample_id) {
  
  pp_s <- build_positive_ppp(sample_id)
  
  dens <- density(
    pp_s,
    sigma = bw.diggle(pp_s)
  )
  
  as.data.frame(dens) %>%
    mutate(Sample = sample_id)
}

density_df <- map_dfr(sample_ids, compute_density_df)
saveRDS(density_df, "cell_density_dataframe.rds")

## Publication raster plot, one sample -----
plot_density_sample <- function(sample_id) {
  
  ggplot() +
    rasterise(
      geom_tile(
        data = density_df %>% filter(Sample == sample_id),
        aes(x = x, y = y, fill = value)
      ),
      dpi = 600
    ) +
    geom_polygon(
      data = boundary_plot_df %>% filter(Sample == sample_id),
      aes(
        x = x,
        y = y,
        group = interaction(Sample, segment)
      ),
      fill = NA,
      color = "black",
      linewidth = 0.35
    ) +
    scale_fill_gradientn(
      colours = brewer.pal(9, "Reds"),
      name = "CD68+ intensity\n(cells / um^2)",
      limits = range(density_df$value)
    ) +
    coord_equal() +
    labs(title = sample_id) +
    base_theme +
    theme(legend.position = "none")
}

density_plots <- map(sample_ids, plot_density_sample)

wrap_plots(density_plots, ncol = 4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


# Normalized contour plot: shared bandwidth across all samples ----------

## Derive one shared sigma from the denser samples only
## bw.diggle is unreliable at low n, so don't let it influence the
## shared value; the point of a fixed sigma is comparability, not
## re-deriving per-sample
reference_sigma <- coords %>%
  filter(Class == "Positive") %>%
  group_by(Sample) %>%
  filter(n() >= 20) %>%
  group_map(function(df_s, key) {
    pp_s <- ppp(
      x = df_s$CentroidX_um,
      y = df_s$CentroidY_um,
      window = tissue_windows[[key$Sample]]
    )
    bw.diggle(pp_s)
  }) %>%
  unlist() %>%
  median()

reference_sigma
saveRDS(reference_sigma, "shared_bw_sigma.rds")

## Recompute density surfaces with the SAME sigma for every sample -----
compute_density_df_fixed <- function(sample_id) {
  
  df_s <- coords %>%
    filter(
      Sample == sample_id,
      Class == "Positive"
    )
  
  win <- tissue_windows[[sample_id]]
  
  if (nrow(df_s) == 0) return(NULL)
  
  pp_s <- ppp(
    x = df_s$CentroidX_um,
    y = df_s$CentroidY_um,
    window = win
  )
  
  dens <- density(
    pp_s,
    sigma = reference_sigma
  )
  
  as.data.frame(dens) %>%
    mutate(Sample = sample_id)
}

density_df_fixed <- map_dfr(sample_ids, compute_density_df_fixed)

shared_breaks <- pretty(range(density_df_fixed$value), n = 10)

shared_limits <- range(density_df_fixed$value)

plot_density_sample_fixed <- function(sample_id) {
  
  ggplot() +
    rasterise(
      geom_tile(
        data = density_df_fixed %>% filter(Sample == sample_id),
        aes(x = x, y = y, fill = value)
      ),
      dpi = 600
    ) +
    geom_polygon(
      data = boundary_plot_df %>% filter(Sample == sample_id),
      aes(
        x = x,
        y = y,
        group = interaction(Sample, segment)
      ),
      fill = NA,
      color = "black",
      linewidth = 0.3
    ) +
    
    # scale_fill_gradientn(
    #   colours = brewer.pal(9, "Reds"),
    #   name = "Cell\ndensity",
    #   limits = shared_limits
    # ) +
    
    scale_fill_viridis_c(
      option = "turbo",
      name = "Cell\ndensity",
      limits = shared_limits
    ) +
    
    coord_equal() +
    labs(title = sample_id) +
    base_theme +
    theme(legend.position = "none")
}

density_plots_fixed <- map(sample_ids, plot_density_sample_fixed)

wrap_plots(density_plots_fixed, ncol = 4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave("Figures/density_plot.svg",
       height = 15, width = 15, units = "cm")