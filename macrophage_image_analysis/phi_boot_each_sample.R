library(dplyr)
library(purrr)
library(furrr)
library(progressr)
library(metafor)
library(emmeans)
library(stringr)

CD68 <- readxl::read_excel(
  "Sarcospan Project CD68 Macrophage Densities.xlsx",
  sheet = 14
)

CD68$Genotype <- factor(CD68$Genotype, levels = c("WT", "mdx", "TG"))

area <- 0.39811499
CD68 <- CD68 %>%
  mutate(Count = round(Value * area, 0))

sample_counts <- split(CD68$Count, CD68$Sample)
animal_geno <- CD68 %>% distinct(Sample, Genotype)
sample_ids <- unique(CD68$Sample)


# Per-animal phi, with bootstrap SE ------------------------------------

compute_phi_animal <- function(vals) {
  m <- mean(vals)
  sum((vals - m)^2 / m) / (length(vals) - 1)
}

compute_phi_with_bootstrap_CI <- function(
    sample_id,
    n_boot = 10000,
    parallel_boot = FALSE
) {
  
  vals <- sample_counts[[sample_id]]
  n_images <- length(vals)
  
  phi_obs <- compute_phi_animal(vals)
  
  one_boot <- function(i) {
    resampled <- sample(vals, n_images, replace = TRUE)
    compute_phi_animal(resampled)
  }
  
  phi_boot <- if (parallel_boot) {
    future_map_dbl(
      1:n_boot,
      one_boot,
      .options = furrr_options(seed = TRUE)
    )
  } else {
    map_dbl(1:n_boot, one_boot)
  }
  
  n_nan <- sum(is.nan(phi_boot))
  if (n_nan > 0) {
    message(
      sample_id, ": ", n_nan, " of ", n_boot,
      " bootstrap iterations produced NaN, excluded from summary."
    )
  }
  
  ## Log-scale SE, computed directly from log-transformed bootstrap
  ## draws (not a delta-method approximation) -- phi's variance
  ## scales with its own magnitude, which badly distorts inverse-
  ## variance weighting on the raw scale (confirmed directly: implied
  ## weights varied ~2400-fold across animals, driven almost entirely
  ## by phi's own size rather than genuine measurement precision).
  ## Log transform is standard for ratio-type statistics (phi/D are
  ## fold-deviations from CSR=1) for exactly this reason.
  log_phi_boot <- log(phi_boot[phi_boot > 0])
  
  tibble(
    Sample = sample_id,
    n_images = n_images,
    phi_obs = phi_obs,
    log_phi_obs = log(phi_obs),
    phi_boot_se = sd(phi_boot, na.rm = TRUE),
    log_phi_boot_se = sd(log_phi_boot, na.rm = TRUE),
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
  
  phi_results <- with_progress({
    p <- progressor(along = sample_ids)
    map_dfr(sample_ids, function(s) {
      result <- compute_phi_with_bootstrap_CI(s, parallel_boot = TRUE)
      p(message = s)
      result
    })
  })
  
} else {
  
  phi_results <- with_progress({
    p <- progressor(along = sample_ids)
    future_map_dfr(
      sample_ids,
      function(s) {
        result <- compute_phi_with_bootstrap_CI(s, parallel_boot = FALSE)
        p(message = s)
        result
      },
      .options = furrr_options(seed = TRUE)
    )
  })
  
}

phi_final <- phi_results %>%
  left_join(animal_geno, by = "Sample") %>%
  rename(Group = Genotype) %>%
  mutate(vi_log = log_phi_boot_se^2)

print(phi_final, n = Inf, width = Inf)


# Random-effects meta-regression on LOG scale, genotype as moderator --------

phi_rma_model <- rma(
  yi = log_phi_obs,
  vi = vi_log,
  mods = ~ Group,
  data = phi_final,
  method = "REML",
  test = "knha"
)

print(phi_rma_model)

ref_grid <- emmprep(phi_rma_model)
phi_emm <- emmeans(ref_grid, ~ Group)
phi_emm

## Manual back-transformation

phi_emm_df <- as.data.frame(phi_emm)
phi_emm_response <- phi_emm_df %>%
  mutate(
    emmean_phi = exp(emmean),
    lower.CL_phi = exp(lower.CL),
    upper.CL_phi = exp(upper.CL)
  )
print(phi_emm_response)

phi_pairwise <- pairs(phi_emm, adjust = "holm")
print(phi_pairwise)

## phi_pairwise's 'estimate' is a difference in log(phi) = log fold-
## change; exponentiate for an interpretable fold-change per contrast
phi_pairwise_df <- as.data.frame(phi_pairwise) %>%
  mutate(fold_change = exp(estimate))
print(phi_pairwise_df)

saveRDS(phi_rma_model, "phi_rma_model.rds")
saveRDS(phi_pairwise, "phi_pairwise_results.rds")

library(ggplot2)

genotype_colors <- c("WT" = "grey10", "mdx" = "red", "TG" = "blue")

base_theme <- theme_minimal(base_size = 8) +
  theme(
    text = element_text(color = "black", size = 8),
    axis.text = element_text(color = "black", size = 8),
    axis.title = element_text(color = "black", size = 8),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.ticks.length.y = unit(-0.1, 'cm'),
    axis.ticks.length.x = unit(0, "cm"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )

phi_emm_response$Group <- factor(phi_emm_response$Group, levels = c("WT", "mdx", "TG"))
phi_final$Group <- factor(phi_final$Group, levels = c("WT", "mdx", "TG"))

ggplot() +
  geom_segment(
    data = phi_emm_response,
    aes(
      y = lower.CL_phi,
      yend = upper.CL_phi,
      x = Group,
      xend = Group
    ),
    linewidth = 0.35
  ) +
  geom_point(
    data = phi_emm_response,
    aes(y = emmean_phi, x = Group),
    size = 3,
    shape = 21,
    fill = "white"
  ) +
  geom_jitter(
    data = phi_final,
    aes(y = phi_obs, x = Group, color = Group),
    width = 0.1,
    alpha = 0.6,
    size = 2,
    shape = 16
  ) +
  scale_color_manual(values = genotype_colors) +
  scale_y_log10() +
  labs(y = "Pearson dispersion (phi)", x = NULL) +
  base_theme +
  theme(legend.position = "none")

ggsave(filename = "Figures/phi_CI_dotplot_log.svg",
       height = 3, width = 3, units = "cm")