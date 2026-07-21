library(tidyverse)
library(boot)
library(Seurat)

FAPs <- readRDS('seurat_FAPs_single_cell.rds')

# 1. per-condition Pearson r, CI, and lm fit --------------------------------

lm_cor <- function(
    seurat_obj, 
    score_x, 
    score_y, 
    condition_col = "condition"
    ) {
  
  data <- FetchData(
    seurat_obj, 
    vars = c(score_x, score_y, condition_col)
    )
  
  colnames(data) <- c("x", "y", "condition")
  
  results <- do.call(
    rbind, 
    lapply(split(data, data$condition), function(df) {
      
    ct <- cor.test(df$x, df$y, method = "pearson")
    fit <- lm(y ~ x, data = df)
    
    data.frame(
      condition = unique(df$condition),
      r = unname(ct$estimate), ci_low = ct$conf.int[1], ci_high = ct$conf.int[2],
      p_value = ct$p.value,
      slope = coef(fit)[2], intercept = coef(fit)[1],
      n = nrow(df)
    )
  }))
  
  rownames(results) <- NULL
  results
}

lm_res <- lm_cor(FAPs, "YAP", "Fibrotic", condition_col = "condition")
print(lm_res)

# scatter with r + CI annotation, lm fit line
plot_lm_scatter <- function(seurat_obj, score_x, score_y, lm_res, condition_col = "condition") {
  
  data <- FetchData(seurat_obj, vars = c(score_x, score_y, condition_col))
  colnames(data) <- c("x", "y", "condition")
  
  labels <- lm_res %>%
    mutate(label = sprintf("r = %.2f [%.2f, %.2f]", r, ci_low, ci_high))
  
  ggplot(data, aes(x = x, y = y)) +
    ggpointdensity::geom_pointdensity(size = 0.5, adjust = 1) +
    scale_color_viridis_c(option = "inferno") +
    geom_smooth(method = "lm", se = TRUE, color = "white", linewidth = 0.5) +
    facet_wrap(~ condition) +
    geom_text(data = labels, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.1, vjust = 1.5, inherit.aes = FALSE, size = 3.5, color = "black") +
    labs(x = score_x, y = score_y, color = "local\ndensity") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())
}

plot_lm_scatter(
  FAPs, 
  "YAP", 
  "Fibrotic", 
  lm_res, 
  condition_col = "condition"
  )

# 2. residuals: flag cells deviating from the fitted relationship ----------

lm_residuals <- function(seurat_obj, score_x, score_y, condition_col = "condition",
                         resid_threshold = 2) {
  
  data <- FetchData(seurat_obj, vars = c(score_x, score_y, condition_col))
  colnames(data) <- c("x", "y", "condition")
  data$cell <- rownames(data)
  
  results <- do.call(rbind, lapply(split(data, data$condition), function(df) {
    fit <- lm(y ~ x, data = df)
    df$fitted <- fitted(fit)
    df$resid_std <- rstandard(fit)
    df$extreme <- abs(df$resid_std) > resid_threshold
    df$direction <- ifelse(df$resid_std > 0, "higher than expected", "lower than expected")
    df
  }))
  rownames(results) <- NULL
  results
}

resid_res <- lm_residuals(FAPs, "YAP", "Fibrotic", condition_col = "condition", resid_threshold = 2)

resid_res %>%
  group_by(condition) %>%
  summarise(n = n(), n_extreme = sum(extreme), frac_extreme = n_extreme / n)

extreme_cells <- resid_res %>% filter(extreme)
print(head(extreme_cells))

plot_lm_residuals <- function(resid_res) {
  ggplot(resid_res, aes(x = x, y = y)) +
    geom_point(data = subset(resid_res, !extreme), color = "grey70",
               size = 0.5, alpha = 0.4) +
    geom_point(data = subset(resid_res, extreme), aes(color = direction),
               size = 1.2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
    facet_wrap(~ condition) +
    scale_color_manual(values = c("higher than expected" = "firebrick",
                                  "lower than expected" = "steelblue")) +
    labs(x = "YAP", y = "Fibrotic", color = NULL) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())
}

plot_lm_residuals(resid_res)

# 3. closed-form eccentricity and mutual information ------------------------
# both monotonic in |r|, so CI bounds on r transform directly -- no
# bootstrap or eigendecomposition needed under the validated Gaussian fit

eccentricity_from_r <- function(r) sqrt(1 - (1 - abs(r)) / (1 + abs(r)))
mi_from_r <- function(r) -0.5 * log(1 - r^2)

gaussian_metrics <- function(lm_res) {
  lm_res %>%
    mutate(
      eccentricity = eccentricity_from_r(r),
      ecc_low  = eccentricity_from_r(ci_low),
      ecc_high = eccentricity_from_r(ci_high),
      mi = mi_from_r(r),
      mi_low  = mi_from_r(ci_low),
      mi_high = mi_from_r(ci_high)
    )
}

metrics_res <- gaussian_metrics(lm_res)
print(metrics_res)

plot_gaussian_metrics <- function(metrics_res) {
  long_df <- rbind(
    data.frame(condition = metrics_res$condition, metric = "Pearson r",
               est = metrics_res$r, low = metrics_res$ci_low, high = metrics_res$ci_high),
    data.frame(condition = metrics_res$condition, metric = "Eccentricity",
               est = metrics_res$eccentricity, low = metrics_res$ecc_low, high = metrics_res$ecc_high),
    data.frame(condition = metrics_res$condition, metric = "Mutual information",
               est = metrics_res$mi, low = metrics_res$mi_low, high = metrics_res$mi_high)
  )
  
  ggplot(long_df, aes(x = est, y = condition)) +
    geom_errorbar(aes(xmin = low, xmax = high), height = 0.15, orientation = "y") +
    geom_point(size = 3) +
    facet_wrap(~ metric, scales = "free_x") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 13)
}

plot_gaussian_metrics(metrics_res)

# 4. pairwise genotype comparison: Fisher z-test on r -----------------------

fisher_z_test <- function(r1, n1, r2, n2) {
  z1 <- atanh(r1); z2 <- atanh(r2)
  se <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  z_stat <- (z1 - z2) / se
  p_value <- 2 * pnorm(-abs(z_stat))
  data.frame(z_stat = z_stat, p_value = p_value)
}

pairwise_condition_comparison <- function(lm_res) {
  pairs <- combn(lm_res$condition, 2, simplify = FALSE)
  do.call(rbind, lapply(pairs, function(p) {
    row1 <- lm_res[lm_res$condition == p[1], ]
    row2 <- lm_res[lm_res$condition == p[2], ]
    ft <- fisher_z_test(row1$r, row1$n, row2$r, row2$n)
    data.frame(condition1 = p[1], condition2 = p[2],
               r1 = row1$r, r2 = row2$r, z_stat = ft$z_stat, p_value = ft$p_value)
  }))
}

condition_comparison <- pairwise_condition_comparison(lm_res)
condition_comparison$p_adj <- p.adjust(condition_comparison$p_value, method = "holm")
print(condition_comparison)

# 5. optional: bootstrap CI on r as a robustness check against Fisher z ----
# run once; if bootstrap and Fisher z CIs agree closely, the Gaussian
# assumption is well-supported and eccentricity/mi CIs above can be trusted
# as-is (they're monotonic transforms of r, so don't need separate bootstrap)

bootstrap_r_check <- function(seurat_obj, score_x, score_y, condition_col = "condition",
                              n_boot = 2000) {
  
  data <- FetchData(seurat_obj, vars = c(score_x, score_y, condition_col))
  colnames(data) <- c("x", "y", "condition")
  
  results <- do.call(rbind, lapply(split(data, data$condition), function(df) {
    cor_fn <- function(d, i) cor(d$x[i], d$y[i], method = "pearson")
    b <- boot(df, cor_fn, R = n_boot)
    ci <- boot.ci(b, type = "perc")
    data.frame(condition = unique(df$condition),
               r_boot = b$t0, boot_ci_low = ci$percent[4], boot_ci_high = ci$percent[5])
  }))
  rownames(results) <- NULL
  results
}

boot_check <- bootstrap_r_check(FAPs, "YAP", "Fibrotic", condition_col = "condition")
comparison_check <- merge(lm_res[, c("condition", "r", "ci_low", "ci_high")], boot_check, by = "condition")
print(comparison_check)
