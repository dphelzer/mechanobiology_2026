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
  
  summary_row <- tibble(
    Sample = sample_id,
    mean_D_obs = mean_D_obs,
    D_boot_mean = mean(boot_mean_D, na.rm = TRUE),
    D_ci_lower = quantile(boot_mean_D, 0.025, na.rm = TRUE),
    D_ci_upper = quantile(boot_mean_D, 0.975, na.rm = TRUE),
    D_boot_se = sd(boot_mean_D, na.rm = TRUE),
    n_nan_dropped = n_nan
  )
  
  list(summary = summary_row, boot_dist = boot_mean_D)
}


# Run, parallelized -----------------------------------------------------

n_cores <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = n_cores)
handlers(global = TRUE)
handlers("txtprogressbar")

use_inner_parallel <- length(sample_ids) < n_cores

if (use_inner_parallel) {
  
  results_list <- with_progress({
    p <- progressor(along = sample_ids)
    map(sample_ids, function(s) {
      result <- compute_mean_D_with_bootstrap_CI(s, parallel_boot = TRUE)
      p(message = s)
      result
    })
  })
  
} else {
  
  results_list <- with_progress({
    p <- progressor(along = sample_ids)
    future_map(
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

names(results_list) <- sample_ids

D_bootstrap_results <- map_dfr(results_list, ~ .x$summary)

boot_dist_df <- imap_dfr(results_list, function(res, s) {
  tibble(Sample = s, D_boot = res$boot_dist)
})

print(D_bootstrap_results, n = Inf, width = Inf)
saveRDS(D_bootstrap_results, file.path(results_dir, "D_bootstrap_results.rds"))
saveRDS(boot_dist_df, file.path(results_dir, "boot_dist_df.rds"))


# Bootstrap distribution overlay, colored by genotype ------------------------

library(ggridges)
library(stringr)
library(patchwork)

boot_dist_df <- boot_dist_df %>%
  mutate(Group = case_when(
    str_starts(Sample, "WT")   ~ "WT",
    str_starts(Sample, "MDX")  ~ "mdx",
    str_starts(Sample, "SSPN") ~ "TG",
    TRUE ~ NA_character_
  ))

stopifnot(!anyNA(boot_dist_df$Group))

summary_df <- boot_dist_df %>%
  group_by(Sample, Group) %>%
  drop_na() %>%
  summarize(
    median = median(D_boot),
    q025 = quantile(D_boot, 0.025),
    q25  = quantile(D_boot, 0.25),
    q75  = quantile(D_boot, 0.75),
    q975 = quantile(D_boot, 0.975),
    .groups = "drop"
  )

genotype_colors <- c("WT" = "grey10", "mdx" = "red", "TG" = "blue")

base_theme <- theme_minimal(base_size = 8) +
  theme(
    text = element_text(color = "black", size = 8),
    axis.text = element_text(color = "black", size = 8),
    axis.title = element_text(color = "black", size = 8),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.ticks.length = unit(-0.1, 'cm'),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )

x_range <- range(c(boot_dist_df$D_boot, summary_df$q025, summary_df$q975), na.rm = TRUE)

## bottom-to-top display order: group by genotype, TG-mdx-WT (matching
## the reference script's convention), samples ordered within group
sample_order <- summary_df %>%
  mutate(Group = factor(Group, levels = c("TG", "mdx", "WT"))) %>%
  arrange(Group, Sample) %>%
  pull(Sample)

boot_dist_df$Sample <- factor(boot_dist_df$Sample, levels = sample_order)
summary_df$Sample <- factor(summary_df$Sample, levels = sample_order)

p_ridge <- ggplot(
  boot_dist_df,
  aes(x = D_boot, y = Sample, fill = Group, color = Group)
) +
  geom_density_ridges(
    alpha = 0.4,
    linewidth = 0.35,
    scale = 1.5,
    rel_min_height = 0.001
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "black",
    linewidth = 0.35
  ) +
  scale_fill_manual(values = genotype_colors) +
  scale_color_manual(values = genotype_colors) +
  scale_x_continuous(limits = x_range) +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(legend.position = "none")

p_box <- ggplot(summary_df, aes(y = Sample, color = Group)) +
  geom_segment(aes(x = q025, xend = q975, yend = Sample), linewidth = 0.35) +
  geom_crossbar(
    aes(xmin = q25, x = median, xmax = q75),
    width = 0.4,
    fill = "white",
    linewidth = 0.35
  ) +
  scale_color_manual(values = genotype_colors) +
  scale_x_continuous(limits = x_range) +
  labs(x = "Index of dispersion (D)", y = NULL) +
  base_theme +
  theme(legend.position = "none")

p_ridge / p_box + plot_layout(heights = c(1, 1))

ggsave(filename = "Figures/bootstrap_quadrant_per_sample.svg",
       height = 10, 
       width = 8, 
       units = "cm")
