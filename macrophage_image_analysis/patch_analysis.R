library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(spatstat)
library(dbscan)
library(sf)
library(furrr)
library(progressr)
library(ggplot2)
library(ggrastr)
library(patchwork)


# Setup --------------------------------------------------------------------

tissue_windows <- readRDS("patch_analysis_results/tissue_windows.rds")
coords <- readRDS("patch_analysis_results/coords.rds")
sample_ids <- names(tissue_windows)

results_dir <- "patch_analysis_results"
if (!dir.exists(results_dir)) dir.create(results_dir)

## Hole-only boundary geometry, re-derived from raw boundary CSVs (not
## preserved in the cached tissue_windows.rds, which has no ring-type
## metadata). Used so eligibility exclusion is based on distance to
## HOLES specifically, not distance to the true outer tissue edge.
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

build_hole_psp <- function(sample_id) {
  
  files_s <- boundary_lookup %>%
    filter(Sample == sample_id) %>%
    pull(file)
  
  b_combined <- map_dfr(files_s, function(f) {
    read_csv(f, show_col_types = FALSE)
  })
  
  hole_rings <- b_combined %>%
    filter(RingType == "hole")
  
  if (nrow(hole_rings) == 0) return(NULL)
  
  segs <- hole_rings %>%
    group_by(PolygonIndex, RingIndex) %>%
    group_split() %>%
    map_dfr(function(ring) {
      tibble(
        x0 = ring$x_um,
        y0 = ring$y_um,
        x1 = c(ring$x_um[-1], ring$x_um[1]),
        y1 = c(ring$y_um[-1], ring$y_um[1])
      )
    })
  
  psp(
    segs$x0, segs$y0, segs$x1, segs$y1,
    window = owin(
      xrange = range(c(segs$x0, segs$x1)),
      yrange = range(c(segs$y0, segs$y1))
    ),
    check = FALSE
  )
}

hole_psp_by_sample <- set_names(
  map(sample_ids, build_hole_psp),
  sample_ids
)

saveRDS(hole_psp_by_sample, file.path(results_dir, "hole_psp_by_sample.rds"))

## Distance from a set of points to the nearest hole boundary. If a
## sample has no holes, everything is fully eligible regardless of
## radius (Inf distance).
dist_to_nearest_hole <- function(x, y, hole_psp) {
  if (is.null(hole_psp)) return(rep(Inf, length(x)))
  pts <- ppp(x, y, window = hole_psp$window, check = FALSE)
  nncross(pts, hole_psp)$dist
}


# Test mode toggle -----------------------------------------------------

TEST_MODE <- FALSE   # FALSE = full run, all samples

if (TEST_MODE) {
  n_sim_L <- 99
  radii <- c(50, 75, 100, 200)
  test_sample_ids <- c("WT_C", "MDX_D", "SSPN_C")
  sample_ids_run <- test_sample_ids
} else {
  n_sim_L <- 999
  radii <- c(50, 60, 70, 80, 90, 100)
  sample_ids_run <- sample_ids
}

run_suffix <- if (TEST_MODE) "test3" else "full"

SKIP_EXPENSIVE_QC <- FALSE  # TRUE skips QC 1/3/4 (radius + envelope sweeps)

save_result <- function(obj, name) {
  path <- file.path(results_dir, paste0(name, "_", run_suffix, ".rds"))
  saveRDS(obj, path)
  message("Saved: ", path)
}


# Core functions -------------------------------------------------------

## Custom replacement for localL(): border correction via closepairs()
## + hole-distance, verified to match localL(correction="none") exactly.
## closepairs(twice=TRUE) already lists each pair from both points'
## perspectives -- tabulate cp$i alone, do NOT combine i and j.
compute_L_excess_per_point <- function(
    x, 
    y, 
    win, 
    radii, 
    hole_psp = NULL) {
  
  n <- length(x)
  if (n < 2) return(rep(NA_real_, n))
  
  X <- ppp(x, y, window = win, check = FALSE)
  area_win <- spatstat.geom::area(win)
  rmax <- max(radii)
  
  cp <- closepairs(X, rmax, what = "ijd")
  hole_dist <- dist_to_nearest_hole(x, y, hole_psp)
  
  excess_mat <- matrix(
    NA_real_,
    nrow = n,
    ncol = length(radii)
  )
  
  for (k in seq_along(radii)) {
    r <- radii[k]
    within_r <- cp$d <= r
    counts <- tabulate(cp$i[within_r], nbins = n)
    K_i <- (area_win / (n - 1)) * counts
    L_i <- sqrt(K_i / pi)
    excess <- L_i - r
    eligible <- hole_dist >= r
    excess[!eligible] <- NA_real_
    excess_mat[, k] <- excess
  }
  
  apply(
    excess_mat, 1,
    function(row) {
      if (all(is.na(row))) return(NA_real_)
      max(row, na.rm = TRUE)
    }
  )
}

## Returns the RAW simulated maxima (not just the 99th percentile), so
## downstream code can bootstrap the threshold estimate itself.
compute_local_L_sim_maxima <- function(
    pp_obs, win, radii, n_sim,
    hole_psp = NULL, parallel_sims = FALSE
) {
  
  one_sim <- function(i) {
    pp_sim <- rpoint(pp_obs$n, win = win)
    excess <- compute_L_excess_per_point(
      pp_sim$x, pp_sim$y, win, radii,
      hole_psp = hole_psp
    )
    if (all(is.na(excess))) return(-Inf)
    max(excess, na.rm = TRUE)
  }
  
  if (parallel_sims) {
    future_map_dbl(
      1:n_sim, one_sim,
      .options = furrr_options(seed = TRUE)
    )
  } else {
    map_dbl(1:n_sim, one_sim)
  }
}

compute_significant_points <- function(sample_id, radii, n_sim, parallel_sims = FALSE) {
  
  win <- tissue_windows[[sample_id]]
  hole_psp <- hole_psp_by_sample[[sample_id]]
  
  df_s <- coords %>%
    filter(
      Sample == sample_id,
      Class == "Positive"
    )
  
  if (nrow(df_s) < 2) {
    return(list(
      sig_coords = df_s[0, c("CentroidX_um", "CentroidY_um")],
      threshold_val = NA_real_,
      sim_maxima = numeric(0),
      point_max_excess = numeric(0)
    ))
  }
  
  pp_obs <- ppp(
    x = df_s$CentroidX_um,
    y = df_s$CentroidY_um,
    window = win
  )
  
  sim_maxima <- compute_local_L_sim_maxima(
    pp_obs, win, radii, n_sim,
    hole_psp = hole_psp,
    parallel_sims = parallel_sims
  )
  
  threshold_val <- quantile(sim_maxima, probs = 0.99)
  
  point_max_excess <- compute_L_excess_per_point(
    pp_obs$x, pp_obs$y, win, radii,
    hole_psp = hole_psp
  )
  
  is_significant <- !is.na(point_max_excess) & (point_max_excess > threshold_val)
  
  if (pp_obs$n != nrow(df_s)) {
    warning(
      sample_id, ": ppp() dropped ", nrow(df_s) - pp_obs$n,
      " point(s) at/outside window boundary"
    )
  }
  
  list(
    sig_coords = tibble(
      CentroidX_um = pp_obs$x[is_significant],
      CentroidY_um = pp_obs$y[is_significant]
    ),
    threshold_val = threshold_val,
    sim_maxima = sim_maxima,
    point_max_excess = point_max_excess
  )
}

## DBSCAN minPts=1 = connected-components on the eps-radius graph
## (pure single-linkage chaining -- long strings of marginal points can
## merge into one elongated patch).
compute_patches_from_significant <- function(sig_coords, linking_radius) {
  
  if (nrow(sig_coords) == 0) {
    return(tibble(
      total_inflamed_area = 0,
      n_patches = 0,
      mean_patch_size = 0
    ))
  }
  
  clusters <- dbscan(
    sig_coords,
    eps = linking_radius,
    minPts = 1
  )
  
  sig_coords$cluster <- clusters$cluster
  
  patch_areas <- sig_coords %>%
    group_by(cluster) %>%
    group_map(function(pts, key) {
      pts_sf <- st_as_sf(pts, coords = c("CentroidX_um", "CentroidY_um"))
      st_area(st_union(st_buffer(pts_sf, dist = linking_radius / 2)))
    }) %>%
    unlist()
  
  tibble(
    total_inflamed_area = sum(patch_areas),
    n_patches = length(patch_areas),
    mean_patch_size = mean(patch_areas)
  )
}

compute_patch_geometries <- function(sig_coords, linking_radius) {
  
  if (nrow(sig_coords) == 0) {
    return(tibble(
      X = numeric(0),
      Y = numeric(0),
      cluster = integer(0),
      group_id = character(0)
    ))
  }
  
  clusters <- dbscan(
    sig_coords,
    eps = linking_radius,
    minPts = 1
  )
  
  sig_coords$cluster <- clusters$cluster
  
  sig_coords %>%
    group_by(cluster) %>%
    group_split() %>%
    map_dfr(function(pts) {
      cl <- unique(pts$cluster)
      pts_sf <- st_as_sf(pts, coords = c("CentroidX_um", "CentroidY_um"))
      unioned <- st_union(st_buffer(pts_sf, dist = linking_radius / 2))
      coords_mat <- as.data.frame(st_coordinates(unioned))
      part_id <- if ("L2" %in% names(coords_mat)) coords_mat$L2 else 1
      
      coords_mat %>%
        mutate(
          cluster = cl,
          group_id = paste(cl, L1, part_id, sep = "_")
        )
    })
}

## Bootstraps the ALREADY-COMPUTED sim_maxima (no resimulation) to get
## an SE on pct_significant driven purely by Monte Carlo noise in the
## threshold estimate.
compute_pct_significant_se <- function(sim_maxima, point_max_excess, n_total, n_boot = 200) {
  
  if (length(sim_maxima) == 0 || length(point_max_excess) == 0 || n_total == 0) {
    return(tibble(se_pct = NA_real_))
  }
  
  boot_pct <- map_dbl(1:n_boot, function(i) {
    boot_sims <- sample(sim_maxima, length(sim_maxima), replace = TRUE)
    boot_threshold <- quantile(boot_sims, 0.99)
    100 * sum(point_max_excess > boot_threshold, na.rm = TRUE) / n_total
  })
  
  tibble(se_pct = sd(boot_pct))
}


# Run: expensive step once per sample -----------------------------------

n_cores <- max(1, parallel::detectCores() - 1)

plan(multisession, workers = n_cores)
handlers(global = TRUE)
handlers("txtprogressbar")

## Fewer samples than cores -> parallelize simulations, not samples,
## so idle cores aren't wasted waiting on the last sample.
use_inner_parallel <- length(sample_ids_run) < n_cores

if (use_inner_parallel) {
  
  significant_points_by_sample <- with_progress({
    p <- progressor(along = sample_ids_run)
    map(sample_ids_run, function(s) {
      result <- compute_significant_points(
        s, radii, n_sim_L,
        parallel_sims = TRUE
      )
      p(message = s)
      result
    })
  }) %>%
    set_names(sample_ids_run)
  
} else {
  
  significant_points_by_sample <- with_progress({
    p <- progressor(along = sample_ids_run)
    future_map(
      sample_ids_run,
      function(s) {
        result <- compute_significant_points(
          s, radii, n_sim_L,
          parallel_sims = FALSE
        )
        p(message = s)
        result
      },
      .options = furrr_options(seed = TRUE)
    )
  }) %>%
    set_names(sample_ids_run)
  
}

save_result(significant_points_by_sample, "significant_points_by_sample")


# Patch sweep: cheap, no resimulation -----------------------------------

linking_radius_candidates <- seq(50, 100, by = 10)

patch_sweep <- map_dfr(linking_radius_candidates, function(lr) {
  map_dfr(sample_ids_run, function(s) {
    compute_patches_from_significant(
      significant_points_by_sample[[s]]$sig_coords,
      lr
    ) %>%
      mutate(
        Sample = s,
        linking_radius = lr
      )
  })
})

save_result(patch_sweep, "patch_sweep")

# Plot styling -----------------------------------------------------------

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
    data.frame(
      x = seg$x,
      y = seg$y,
      segment = i,
      Sample = s
    )
  })
})

sig_points_df <- imap_dfr(significant_points_by_sample, function(res, s) {
  if (nrow(res$sig_coords) == 0) return(NULL)
  res$sig_coords %>% mutate(Sample = s)
})


# QC: overlay significant points + patches on coordinate plots ----------

## Linking radius = max(radii): the main analysis combines multiple
## radii into one test, so a fixed small linking radius would fragment
## patches found at the larger end of that range.

# overlay_linking_radius_main <- max(radii)
overlay_linking_radius_main <- (135)

plot_sig_overlay_sample <- function(sample_id) {
  
  sig_s <- sig_points_df %>%
    filter(Sample == sample_id)
  
  patch_geom_df <- compute_patch_geometries(
    sig_s %>% select(-Sample),
    overlay_linking_radius_main
  )
  
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
      linewidth = 0.35
    ) +
    rasterise(
      geom_point(
        data = coords %>% filter(Sample == sample_id),
        aes(
          x = CentroidX_um,
          y = CentroidY_um,
          color = Class
        ),
        size = 0.01,
        alpha = 0.3
      ),
      dpi = 600
    ) +
    geom_polygon(
      data = patch_geom_df,
      aes(
        x = X,
        y = Y,
        group = group_id,
        fill = factor(cluster)
      ),
      alpha = 0.85,
      color = "black",
      linewidth = 0.35
    ) +
    scale_fill_viridis_d(option = "turbo") +
    scale_color_manual(values = class_colors) +
    coord_equal() +
    labs(title = sample_id) +
    base_theme +
    theme(legend.position = "none")
}

sig_overlay_plots <- map(sample_ids_run, plot_sig_overlay_sample)

wrap_plots(sig_overlay_plots, ncol = 4)

## save plot
ggsave("Figures/clusters.svg",
       height = 15, width = 15, units = "cm")

# QC 1 + QC 4: significant points across a (sample x radius) grid -------
# GATED: validates radii choice on the test set.

if (TEST_MODE && !SKIP_EXPENSIVE_QC) {
  
  
  radius_grid <- seq(50, 100, by = 10)
  
  qc_grid <- expand.grid(
    Sample = sample_ids_run,
    radius = radius_grid,
    stringsAsFactors = FALSE
  )
  
  sig_by_sample_radius <- with_progress({
    p <- progressor(nrow(qc_grid))
    future_map2(
      qc_grid$Sample, qc_grid$radius,
      function(s, r) {
        result <- compute_significant_points(s, radii = r, n_sim = n_sim_L)
        p()
        result
      },
      .options = furrr_options(seed = TRUE)
    )
  })
  
  radius_stability <- qc_grid %>%
    mutate(
      n_total = map_dbl(Sample, function(s) {
        nrow(coords %>% filter(Sample == s, Class == "Positive"))
      }),
      n_sig = map_dbl(sig_by_sample_radius, function(res) nrow(res$sig_coords)),
      n_eligible = map2_dbl(Sample, radius, function(s, r) {
        hole_psp <- hole_psp_by_sample[[s]]
        df_s <- coords %>% filter(Sample == s, Class == "Positive")
        if (nrow(df_s) < 2) return(0)
        sum(dist_to_nearest_hole(df_s$CentroidX_um, df_s$CentroidY_um, hole_psp) >= r)
      }),
      pct_significant = ifelse(n_total > 0, 100 * n_sig / n_total, 0),
      pct_eligible = ifelse(n_total > 0, 100 * n_eligible / n_total, 0),
      se_pct = map2_dbl(sig_by_sample_radius, n_total, function(res, nt) {
        compute_pct_significant_se(res$sim_maxima, res$point_max_excess, nt)$se_pct
      }),
      threshold_val = map_dbl(sig_by_sample_radius, function(res) res$threshold_val)
    ) %>%
    select(Sample, radius, pct_significant, se_pct, pct_eligible, threshold_val)
  
  save_result(sig_by_sample_radius, "sig_by_sample_radius")
  save_result(radius_stability, "radius_stability")
  
  radius_stability %>%
    group_by(Sample) %>%
    summarize(
      min_eligible = min(pct_eligible),
      max_eligible = max(pct_eligible),
      any_na_ntotal = anyNA(pct_eligible)
    ) %>%
    print()
  
  p_significant <- ggplot(
    radius_stability,
    aes(radius, pct_significant, color = Sample, fill = Sample)
  ) +
    geom_ribbon(
      aes(
        ymin = pmax(0, pct_significant - se_pct),
        ymax = pmin(100, pct_significant + se_pct)
      ),
      alpha = 0.2,
      color = NA
    ) +
    geom_line(linewidth = 0.35) +
    geom_point(size = 1) +
    labs(
      x = "Radius (um)",
      y = "% flagged significant (+/- SE)"
    ) +
    base_theme +
    theme(
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.position = "right"
    )
  
  p_eligible <- ggplot(
    radius_stability,
    aes(radius, pct_eligible, color = Sample)
  ) +
    geom_line(linewidth = 0.35) +
    geom_point(size = 1) +
    labs(
      x = "Radius (um)",
      y = "% eligible (bdist >= r)"
    ) +
    base_theme +
    theme(
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.position = "right"
    )
  

  p_threshold <- ggplot(
    radius_stability,
    aes(radius, threshold_val, color = Sample)
  ) +
    geom_line(linewidth = 0.35) +
    geom_point(size = 1) +
    labs(
      x = "Radius (um)",
      y = "Significance threshold (L excess)"
    ) +
    base_theme +
    theme(
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.position = "right"
    )
  
  print(
    (p_significant / p_eligible / p_threshold) +
      plot_layout(guides = "collect") &
      theme(legend.position = "right")
  )
  
  select_radius_criteria <- function(radius_stability, eligible_floor = 50, min_consecutive = 3) {
    
    eligible_summary <- radius_stability %>%
      filter(pct_eligible >= eligible_floor) %>%
      group_by(Sample) %>%
      summarize(max_eligible_radius = max(radius), .groups = "drop")
    
    eligible_summary <- tibble(Sample = unique(radius_stability$Sample)) %>%
      left_join(eligible_summary, by = "Sample")
    
    plateau_summary <- radius_stability %>%
      group_by(Sample) %>%
      arrange(radius, .by_group = TRUE) %>%
      mutate(
        delta_sig = pct_significant - lag(pct_significant),
        quiet_step = !is.na(delta_sig) & abs(delta_sig) <= se_pct,
        quiet_run = ave(quiet_step, cumsum(!quiet_step), FUN = cumsum)
      ) %>%
      summarize(
        plateau_radius = radius[which(quiet_run >= min_consecutive)[1]],
        .groups = "drop"
      )
    
    eligible_summary %>%
      left_join(plateau_summary, by = "Sample")
  }
  
  radius_criteria <- select_radius_criteria(radius_stability)
  print(radius_criteria)
  save_result(radius_criteria, "radius_criteria")
  
} else {
  message("Skipping QC 1+4 -- TEST_MODE is FALSE.")
}


# QC 2: DBSCAN linking-radius stability -----------------------------------
## Cheap, no resimulation. Narrowed to 50-150, finer resolution.

linking_radius_grid <- seq(50, 200, by = 10)

linking_stability <- map_dfr(linking_radius_grid, function(lr) {
  map_dfr(sample_ids_run, function(s) {
    compute_patches_from_significant(
      significant_points_by_sample[[s]]$sig_coords,
      lr
    ) %>%
      mutate(
        Sample = s,
        linking_radius = lr
      )
  })
})

## Smallest linking radius at which n_patches has stopped decreasing
## meaningfully -- n_patches stabilizing is a cleaner signal than
## mean_patch_size, which keeps growing by construction even after
## patch count stops changing.
select_linking_radius_criteria <- function(linking_stability) {
  linking_stability %>%
    group_by(Sample) %>%
    arrange(linking_radius, .by_group = TRUE) %>%
    summarize(
      stable_linking_radius = linking_radius[which(n_patches == min(n_patches))[1]],
      .groups = "drop"
    )
}

linking_radius_criteria <- select_linking_radius_criteria(linking_stability)
print(linking_radius_criteria)
save_result(linking_radius_criteria, "linking_radius_criteria")

p_npatches <- ggplot(
  linking_stability,
  aes(linking_radius, n_patches, color = Sample)
) +
  geom_line(linewidth = 0.35) +
  geom_point(size = 1) +
  labs(
    x = "Linking radius (um)",
    y = "Number of patches"
  ) +
  base_theme +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  )

p_meansize <- ggplot(
  linking_stability,
  aes(linking_radius, mean_patch_size, color = Sample)
) +
  geom_line(linewidth = 0.35) +
  geom_point(size = 1) +
  labs(
    x = "Linking radius (um)",
    y = "Mean patch size (um^2)"
  ) +
  base_theme +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  )

(p_npatches / p_meansize) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

save_result(linking_stability, "linking_stability")


# QC 3: global L(r) - r with CSR envelope, per sample --------------------
# GATED, same reasoning as QC 1+4.

compute_L_envelope_df <- function(sample_id, nsim = 499) {
  
  win <- tissue_windows[[sample_id]]
  
  df_s <- coords %>%
    filter(
      Sample == sample_id,
      Class == "Positive"
    )
  
  if (nrow(df_s) < 2) return(NULL)
  
  pp_obs <- ppp(
    x = df_s$CentroidX_um,
    y = df_s$CentroidY_um,
    window = win
  )
  
  env <- envelope(
    pp_obs, Lest,
    nsim = nsim,
    correction = "translate",
    fix.n = TRUE,
    verbose = FALSE
  )
  
  tibble(
    Sample = sample_id,
    r = env$r,
    obs_excess = env$obs - env$r,
    lo_excess = env$lo - env$r,
    hi_excess = env$hi - env$r
  )
}

if (TEST_MODE && !SKIP_EXPENSIVE_QC) {
  
  L_envelope_df <- with_progress({
    p <- progressor(along = sample_ids_run)
    future_map_dfr(
      sample_ids_run,
      function(s) {
        result <- compute_L_envelope_df(s)
        p(message = s)
        result
      },
      .options = furrr_options(seed = TRUE)
    )
  })
  
  save_result(L_envelope_df, "L_envelope_df")
  
  print(
    ggplot(L_envelope_df, aes(x = r)) +
      geom_ribbon(
        aes(ymin = lo_excess, ymax = hi_excess),
        fill = "grey80",
        alpha = 0.6
      ) +
      geom_line(
        aes(y = obs_excess),
        color = "red",
        linewidth = 0.4
      ) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.3
      ) +
      facet_wrap(~ Sample, scales = "free") +
      labs(
        x = "r (um)",
        y = "L(r) - r"
      ) +
      base_theme +
      theme(
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8)
      )
  )
  
  ## Same data, overlaid on one axis instead of faceted -- direct
  ## visual comparison of L(r)-r across samples. Envelope ribbon
  ## dropped here (would overlap messily with 3+ samples); the
  ## faceted version above still shows each sample's own envelope.
  print(
    ggplot(L_envelope_df, aes(x = r, y = obs_excess, color = Sample)) +
      geom_line(linewidth = 0.5) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.3
      ) +
      labs(
        x = "r (um)",
        y = "L(r) - r"
      ) +
      base_theme +
      theme(
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = "right"
      )
  )
  
} else {
  message("Skipping QC 3 -- TEST_MODE is FALSE.")
}


# QC 4: overlay grid across radius sweep ----------------------------------
# GATED -- depends on sig_by_sample_radius from QC 1.

if (TEST_MODE && !SKIP_EXPENSIVE_QC) {
  
  ## Small representative subset spanning 50-100. Must use values that
  ## exist EXACTLY in radius_grid (subsetted, not recomputed).
  radius_grid_qc4 <- c(50, 70, 90, 100)
  
  missing_r <- setdiff(radius_grid_qc4, radius_grid)
  if (length(missing_r) > 0) {
    stop(
      "radius_grid_qc4 contains values not in radius_grid: ",
      paste(missing_r, collapse = ", ")
    )
  }
  
  sig_by_radius <- map(radius_grid_qc4, function(r) {
    idx <- which(qc_grid$radius == r)
    set_names(
      map(idx, function(i) sig_by_sample_radius[[i]]$sig_coords),
      qc_grid$Sample[idx]
    )
  }) %>%
    set_names(as.character(radius_grid_qc4))
  
  plot_sig_overlay_at_radius <- function(sample_id, r) {
    
    sig_r <- sig_by_radius[[as.character(r)]][[sample_id]]
    
    patch_geom_df <- compute_patch_geometries(sig_r, linking_radius = r)
    
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
          aes(
            x = CentroidX_um,
            y = CentroidY_um,
            color = Class
          ),
          size = 0.15,
          alpha = 0.3
        ),
        dpi = 600
      ) +
      geom_polygon(
        data = patch_geom_df,
        aes(
          x = X,
          y = Y,
          group = group_id,
          fill = factor(cluster)
        ),
        alpha = 0.45,
        color = "black",
        linewidth = 0.6
      ) +
      scale_fill_viridis_d(option = "turbo") +
      scale_color_manual(values = class_colors) +
      coord_equal() +
      labs(title = paste(sample_id, "| r =", r, "| n_sig =", nrow(sig_r))) +
      base_theme +
      theme(legend.position = "none")
  }
  
  overlay_grid_plots <- map(sample_ids_run, function(s) {
    map(radius_grid_qc4, function(r) plot_sig_overlay_at_radius(s, r))
  }) %>%
    flatten()
  
  print(wrap_plots(overlay_grid_plots, ncol = length(radius_grid_qc4)))
  
} else {
  message("Skipping QC 4 -- TEST_MODE is FALSE.")
}


# Summary stats ---------------------------------------------------------
compute_auc <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

auc_summary <- patch_sweep %>%
  group_by(Sample) %>%
  arrange(linking_radius, .by_group = TRUE) %>%
  summarize(
    auc_n_patches = compute_auc(linking_radius, n_patches),
    auc_mean_patch_size = compute_auc(linking_radius, mean_patch_size),
    range_width = max(linking_radius) - min(linking_radius),
    mean_n_patches = auc_n_patches / range_width,
    mean_patch_size = auc_mean_patch_size / range_width,
    .groups = "drop"
  )

auc_summary <- auc_summary %>%
  mutate(Group = case_when(
    str_starts(Sample, "WT")   ~ "WT",
    str_starts(Sample, "MDX")  ~ "mdx",
    str_starts(Sample, "SSPN") ~ "TG",
    TRUE ~ NA_character_
  ))

# safety check: confirm every sample got labeled
stopifnot(!anyNA(auc_summary$Group))

auc_summary$Group <- factor(
  auc_summary$Group,
  levels = c("WT", "mdx", "TG")
  )

auc_summary %>% 
  ggplot(aes(x=Group, y=mean_n_patches)) +
  geom_boxplot() +
  geom_jitter() +
  theme_minimal()

auc_summary %>% 
  ggplot(aes(x=Group, y=mean_patch_size)) +
  geom_boxplot() +
  geom_jitter() +
  theme_minimal()
openxlsx::write.xlsx(auc_summary, "patch_analysis.xlsx")
save_result(auc_summary, "auc_summary")


# simple average

patch_sweep <- readRDS("patch_analysis_results/patch_sweep_full.rds")

simple_summary <- patch_sweep %>%
  group_by(Sample) %>%
  summarize(
    mean_n_patches = mean(n_patches),
    mean_patch_size = mean(mean_patch_size),
    .groups = "drop"
  )

simple_summary <- simple_summary %>%
  mutate(Group = case_when(
    str_starts(Sample, "WT")   ~ "WT",
    str_starts(Sample, "MDX")  ~ "mdx",
    str_starts(Sample, "SSPN") ~ "TG",
    TRUE ~ NA_character_
  ))

# safety check: confirm every sample got labeled
stopifnot(!anyNA(simple_summary$Group))

simple_summary$Group <- factor(
  simple_summary$Group,
  levels = c("WT", "mdx", "TG")
)

simple_summary %>% 
  ggplot(aes(x=Group, y=mean_n_patches)) +
  geom_boxplot() +
  geom_jitter() +
  theme_minimal()

simple_summary %>% 
  ggplot(aes(x=Group, y=mean_patch_size)) +
  geom_boxplot() +
  geom_jitter() +
  theme_minimal()
openxlsx::write.xlsx(simple_summary, "mean_results.xlsx")

