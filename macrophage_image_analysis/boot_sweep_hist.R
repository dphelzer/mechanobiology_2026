library(dplyr)
library(purrr)
library(spatstat)
library(ggplot2)
library(patchwork)


# Load data ----------------------------------------------------------------

results_dir <- "patch_analysis_results"

tissue_windows <- readRDS(file.path(results_dir, "tissue_windows.rds"))
coords <- readRDS(file.path(results_dir, "coords.rds"))
sample_ids <- names(tissue_windows)


# Index of dispersion (D), parametrized by bin size -------------------------

compute_D <- function(
    sample_id,
    bin_size
) {
  
  win <- tissue_windows[[sample_id]]
  
  df_s <- coords %>%
    filter(
      Sample == sample_id,
      Class == "Positive"
    )
  
  if (nrow(df_s) < 2) {
    return(tibble(
      Sample = sample_id,
      bin_size = bin_size,
      D = NA_real_,
      n_bins = NA_integer_,
      mean_expected_count = NA_real_
    ))
  }
  
  nx <- ceiling(diff(win$xrange) / bin_size)
  ny <- ceiling(diff(win$yrange) / bin_size)
  
  xbreaks <- seq(win$xrange[1], win$xrange[2], length.out = nx + 1)
  ybreaks <- seq(win$yrange[1], win$yrange[2], length.out = ny + 1)
  
  xi <- cut(
    df_s$CentroidX_um,
    breaks = xbreaks,
    include.lowest = TRUE,
    labels = FALSE
  )
  yi <- cut(
    df_s$CentroidY_um,
    breaks = ybreaks,
    include.lowest = TRUE,
    labels = FALSE
  )
  
  counts <- as.vector(table(
    factor(xi, levels = 1:nx),
    factor(yi, levels = 1:ny)
  ))
  
  xcenters <- (xbreaks[-length(xbreaks)] + xbreaks[-1]) / 2
  ycenters <- (ybreaks[-length(ybreaks)] + ybreaks[-1]) / 2
  bin_centers <- expand.grid(x = xcenters, y = ycenters)
  bin_inside <- inside.owin(bin_centers$x, bin_centers$y, win)
  counts <- counts[bin_inside]
  
  tibble(
    Sample = sample_id,
    bin_size = bin_size,
    D = var(counts) / mean(counts),
    n_bins = length(counts),
    mean_expected_count = sum(counts) / length(counts)
  )
}


# Sweep ----------------------------------------------------------------

bin_size_grid <- seq(25, 300, by = 25)

qc_grid <- expand.grid(
  Sample = sample_ids,
  bin_size = bin_size_grid,
  stringsAsFactors = FALSE
)

quadrat_sweep <- map2_dfr(
  qc_grid$Sample,
  qc_grid$bin_size,
  compute_D
)

print(quadrat_sweep, n = Inf)
saveRDS(quadrat_sweep, file.path(results_dir, "quadrat_sweep.rds"))


# QC plots ----------------------------------------------------------------

p_D <- ggplot(
  quadrat_sweep,
  aes(bin_size, D, color = Sample)
) +
  geom_line(linewidth = 0.35) +
  geom_point(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3) +
  labs(
    x = "Bin size (um)",
    y = "Index of dispersion (D = var/mean)"
  ) +
  theme_minimal(base_size = 8)

p_nbins <- ggplot(
  quadrat_sweep,
  aes(bin_size, n_bins, color = Sample)
) +
  geom_line(linewidth = 0.35) +
  geom_point(size = 1) +
  labs(
    x = "Bin size (um)",
    y = "Number of bins"
  ) +
  theme_minimal(base_size = 8)

print(
  (p_D / p_nbins) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
)


# Final sweep: 50-100 by 10, averaged per sample -----------------------------

bin_size_final <- seq(50, 100, by = 10)

qc_grid_final <- expand.grid(
  Sample = sample_ids,
  bin_size = bin_size_final,
  stringsAsFactors = FALSE
)

quadrat_sweep_final <- map2_dfr(
  qc_grid_final$Sample,
  qc_grid_final$bin_size,
  compute_D
)

D_summary <- quadrat_sweep_final %>%
  group_by(Sample) %>%
  summarize(mean_D = mean(D, na.rm = TRUE), .groups = "drop")

print(D_summary)
saveRDS(D_summary, file.path(results_dir, "D_summary.rds"))