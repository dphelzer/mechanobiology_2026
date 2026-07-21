library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(spatstat)
library(patchwork)
library(ggrastr)


# Tissue boundaries and coordinate plots -------------------------------

## Load boundary files per sample -----
## (handles multiple annotation files per sample)
boundary_files <- list.files(
  "tissue_boundaries",
  pattern = "_boundary_ann[0-9]+\\.csv$",
  full.names = TRUE
)

boundary_sample_ids <- str_remove(
  basename(boundary_files),
  "_counts_boundary_ann[0-9]+\\.csv$"
)

boundary_lookup <- tibble(
  file = boundary_files,
  Sample = boundary_sample_ids
)

## Load cell detection coordinates -----
files <- list.files(
  "coordinate_exports",
  pattern = "_detections.csv$",
  full.names = TRUE
)

coords <- map_dfr(
  files,
  read_tsv,
  show_col_types = FALSE
)

coords <- coords %>%
  dplyr::mutate(
    Sample = str_remove(Image, "_counts.*$"),
    CentroidX_um = as.numeric(CentroidX_um),
    CentroidY_um = as.numeric(CentroidY_um)
  )

## Build tissue_window for one sample -----
## (combines all boundary files for that sample; falls back to reversed
## winding order if the initial polygon traversal direction is invalid)
build_sample_window <- function(sample_id) {
  
  files_s <- boundary_lookup %>%
    dplyr::filter(Sample == sample_id) %>%
    pull(file)
  
  b_combined <- map_dfr(files_s, function(f) {
    read_csv(f, show_col_types = FALSE) %>%
      mutate(SourceFile = f)
  })
  
  poly_list_asis <- b_combined %>%
    dplyr::group_by(SourceFile, PolygonIndex, RingType, RingIndex) %>%
    dplyr::group_split() %>%
    lapply(function(ring) {
      list(x = ring$x_um, y = ring$y_um)
    })
  
  result <- tryCatch(
    owin(poly = poly_list_asis),
    error = function(e) NULL
  )
  
  if (is.null(result)) {
    poly_list_rev <- b_combined %>%
      dplyr::group_by(SourceFile, PolygonIndex, RingType, RingIndex) %>%
      dplyr::group_split() %>%
      lapply(function(ring) {
        list(x = rev(ring$x_um), y = rev(ring$y_um))
      })
    
    result <- tryCatch(
      owin(poly = poly_list_rev),
      error = function(e) NULL
    )
  }
  
  result
}

sample_ids <- unique(boundary_lookup$Sample)

tissue_windows <- set_names(
  map(sample_ids, build_sample_window),
  sample_ids
)

## flag any samples where window construction failed entirely
failed_samples <- sample_ids[map_lgl(tissue_windows, is.null)]
failed_samples


## QC: percent of cells falling outside their own tissue window -----
qc_results <- map_dfr(sample_ids, function(s) {
  
  win <- tissue_windows[[s]]
  
  if (is.null(win)) {
    return(tibble(
      Sample = s,
      n_cells = NA,
      n_outside = NA,
      pct_outside = NA,
      status = "WINDOW FAILED"
    ))
  }
  
  df_s <- coords %>% filter(Sample == s)
  n_total <- nrow(df_s)
  n_outside <- sum(!inside.owin(
    df_s$CentroidX_um,
    df_s$CentroidY_um,
    win
  ))
  
  tibble(
    Sample = s,
    n_cells = n_total,
    n_outside = n_outside,
    pct_outside = round(100 * n_outside / n_total, 2),
    status = "ok"
  )
})

print(qc_results, n = Inf)

## samples needing a closer look
qc_results %>%
  filter(status != "ok" | pct_outside > 5)


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


## Publication-quality coordinate plot -----
## (patchwork version -- correct 1:1 aspect ratio preserved per sample)
plot_coords_sample <- function(sample_id) {
  
  ggplot() +
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
    rasterise(
      geom_point(
        data = coords %>% filter(Sample == sample_id),
        aes(x = CentroidX_um, y = CentroidY_um, color = Class),
        size = 0.01,
        alpha = 0.4
      ),
      dpi = 600
    ) +
    scale_color_manual(values = class_colors) +
    coord_equal() +
    labs(title = sample_id) +
    base_theme +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    )
}

coord_plots <- map(sample_ids, plot_coords_sample)

wrap_plots(coord_plots, ncol = 4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "none")

ggsave("Figures/coordinate_plot.svg",
       height = 15, width = 15, units = "cm")
