library(dplyr)
library(purrr)
library(spatstat)
library(furrr)
library(progressr)

results_dir <- "patch_analysis_results"

tissue_windows <- readRDS(file.path(results_dir, "tissue_windows.rds"))
coords <- readRDS(file.path(results_dir, "coords.rds"))
sample_ids <- names(tissue_windows)


# Bin counts for one sample at one bin size ----------------------------------

bin_counts_for_sample <- function(
    sample_id,
    bin_size
) {
  
  win <- tissue_windows[[sample_id]]
  
  df_s <- coords %>%
    filter(
      Sample == sample_id,
      Class == "Positive"
    )
  
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
  
  counts[bin_inside]
}


# Mean D across the 50-100 sweep, with bootstrap CI ---------------------------

compute_mean_D_with_bootstrap_CI <- function(
    sample_id,
    bin_size_grid = seq(50, 100, by = 10),
    n_boot = 10000,
    parallel_boot = FALSE
) {
  
  counts_by_binsize <- map(bin_size_grid, function(bs) {
    bin_counts_for_sample(sample_id, bs)
  })
  
  D_obs_vec <- map_dbl(counts_by_binsize, function(cnts) var(cnts) / mean(cnts))
  mean_D_obs <- mean(D_obs_vec)
  
  one_boot <- function(i) {
    D_vals <- map_dbl(counts_by_binsize, function(cnts) {
      resampled <- sample(cnts, length(cnts), replace = TRUE)
      var(resampled) / mean(resampled)
    })
    mean(D_vals)
  }
  
  boot_mean_D <- if (parallel_boot) {
    future_map_dbl(
      1:n_boot,
      one_boot,
      .options = furrr_options(seed = TRUE)
    )
  } else {
    map_dbl(1:n_boot, one_boot)
  }
  
  n_nan <- sum(is.nan(boot_mean_D))
  if (n_nan > 0) {
    message(
      sample_id, ": ", n_nan, " of ", n_boot,
      " bootstrap iterations produced NaN, excluded from summary."
    )
  }
  
  tibble(
    Sample = sample_id,
    mean_D_obs = mean_D_obs,
    D_boot_mean = mean(boot_mean_D, na.rm = TRUE),
    D_ci_lower = quantile(boot_mean_D, 0.025, na.rm = TRUE),
    D_ci_upper = quantile(boot_mean_D, 0.975, na.rm = TRUE),
    D_boot_se = sd(boot_mean_D, na.rm = TRUE),
    n_nan_dropped = n_nan
  )
}


# Run, parallelized -----------------------------------------------------

n_cores <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = n_cores)
handlers(global = TRUE)
handlers("txtprogressbar")

use_inner_parallel <- length(sample_ids) < n_cores

if (use_inner_parallel) {
  
  D_bootstrap_results <- with_progress({
    p <- progressor(along = sample_ids)
    map_dfr(sample_ids, function(s) {
      result <- compute_mean_D_with_bootstrap_CI(s, parallel_boot = TRUE)
      p(message = s)
      result
    })
  })
  
} else {
  
  D_bootstrap_results <- with_progress({
    p <- progressor(along = sample_ids)
    future_map_dfr(
      sample_ids,
      function(s) {
        result <- compute_mean_D_with_bootstrap_CI(s, parallel_boot = FALSE)
        p(message = s)
        result
      },
      .options = furrr_options(seed = TRUE)
    )
  })
  
}

print(D_bootstrap_results, n = Inf, width = Inf)
saveRDS(D_bootstrap_results, file.path(results_dir, "D_bootstrap_results.rds"))