library(tidyverse)
library(SpatialFeatureExperiment)

sfe_list        <- readRDS("spatial_corr_SFE.rds")
sfe_scores_list <- readRDS("sfe_scroes_list.rds")
score_cols <- c("YAP_UCell", "Inflammatory_UCell", "Fibrotic_UCell")
pairs <- list(
  c("YAP_UCell", "Inflammatory_UCell"),
  c("YAP_UCell", "Fibrotic_UCell"),
  c("Inflammatory_UCell", "Fibrotic_UCell")
)
BPPARAM <- BiocParallel::SnowParam(workers = parallel::detectCores() - 1)
# local Moran's I -----------------------------------------------------------
for (nm in names(sfe_list)) {
  sfe_list[[nm]] <- colDataUnivariate(
    sfe_list[[nm]], type = "localmoran_perm", features = score_cols,
    colGraphName = "visium", nsim = 99999, zero.policy = TRUE,
    BPPARAM = BPPARAM
  )
}
saveRDS(sfe_list, "SFE_Moran_I_local.rds")
# check column names once before trusting the extraction below
# colnames(localResult(sfe_list[["WT"]], "localmoran_perm", "Inflammatory_UCell"))
local_moran_df <- do.call(rbind, lapply(names(sfe_list), function(nm) {
  do.call(rbind, lapply(score_cols, function(s) {
    lm_res <- as.data.frame(localResult(sfe_list[[nm]], "localmoran_perm", s))
    data.frame(
      genotype  = nm,
      signature = s,
      spot      = rownames(lm_res),
      local_I   = lm_res$Ii,
      p_value   = lm_res[["Pr(z != E(Ii))"]],
      p_sim     = lm_res[["Pr(z != E(Ii)) Sim"]]
    )
  }))
}))
saveRDS(local_moran_df, "local_moran_df.rds")
# openxlsx::write.xlsx(local_moran_df, "local_moran_I.xlsx")
# local Lee's L ---------------------------------------------------------
local_lee_results <- list()
for (nm in names(sfe_scores_list)) {
  sfe <- sfe_scores_list[[nm]]
  for (pr in pairs) {
    f1 <- pr[1]; f2 <- pr[2]
    key <- paste(nm, f1, f2, sep = "__")
    
    ll <- calculateBivariate(
      sfe, type = "locallee", feature1 = f1, feature2 = f2,
      colGraphName = "visium", exprs_values = "scores",
      zero.policy = TRUE
    )
    
    ll_df <- as.data.frame(ll)
    
    local_lee_results[[key]] <- data.frame(
      genotype = nm, feature1 = f1, feature2 = f2,
      spot = colnames(sfe), local_L = ll_df[[1]]
    )
  }
}
local_lee_df <- do.call(rbind, local_lee_results)
rownames(local_lee_df) <- NULL
saveRDS(local_lee_df, "local_lee_df.rds")
# openxlsx::write.xlsx(local_lee_df, "local_lee_L.xlsx")

# hotspot fraction (raw p, fixed threshold across all groups) -----------
# fixed threshold rather than per-group FDR: BH-FDR's rejection cutoff is
# data-adaptive, so it can shift more leniently for genotypes with more
# true signal, inflating apparent between-genotype differences in hotspot
# fraction beyond the real effect size. spot counts are similar across
# genotypes here, so a fixed raw threshold keeps the comparison fair.

local_moran_df <- local_moran_df %>%
  group_by(signature) %>%              # pool across genotype, keep score separate
  mutate(p_adj = p.adjust(p_sim, method = "fdr")) %>%
  ungroup()

# p_threshold <- 0.05
# local_moran_df$significant <- local_moran_df$p_value < p_threshold

local_moran_df$significant <- local_moran_df$p_adj < 0.05

hotspot_fraction <- local_moran_df %>%
  group_by(genotype, signature) %>%
  summarise(n_spots = n(),
            n_na = sum(is.na(significant)),
            n_significant = sum(significant, na.rm = TRUE),
            frac_hotspot = n_significant / (n_spots - n_na),
            .groups = "drop")

print(hotspot_fraction)
openxlsx::write.xlsx(hotspot_fraction, "hotspot_fraction.xlsx")

# patch count + area percentage -----------------------------------------

count_patches <- function(sfe, genotype_name, score, local_df, graph_name = "visium") {
  g  <- colGraph(sfe, graph_name)
  nb <- g$neighbours
  n_total <- ncol(sfe)
  
  sig_spots <- local_df$spot[local_df$genotype == genotype_name &
                               local_df$signature == score &
                               !is.na(local_df$significant) & local_df$significant]
  spot_idx <- match(sig_spots, colnames(sfe))
  spot_idx <- spot_idx[!is.na(spot_idx)]
  
  if (length(spot_idx) == 0) {
    return(data.frame(genotype = genotype_name, signature = score, n_patches = 0,
                      mean_patch_size = NA, max_patch_size = NA,
                      total_area_pct = 0))
  }
  
  visited <- setNames(rep(FALSE, length(spot_idx)), spot_idx)
  n_patches <- 0L
  patch_sizes <- c()
  
  for (idx in spot_idx) {
    if (visited[as.character(idx)]) next
    n_patches <- n_patches + 1L
    queue <- idx
    members <- c()
    while (length(queue) > 0) {
      node <- queue[1]; queue <- queue[-1]
      if (visited[as.character(node)]) next
      visited[as.character(node)] <- TRUE
      members <- c(members, node)
      neigh <- nb[[node]]
      neigh <- neigh[neigh %in% spot_idx]
      neigh <- neigh[!visited[as.character(neigh)]]
      queue <- c(queue, neigh)
    }
    patch_sizes <- c(patch_sizes, length(members))
  }
  
  data.frame(genotype = genotype_name, signature = score, n_patches = n_patches,
             mean_patch_size = mean(patch_sizes), max_patch_size = max(patch_sizes),
             total_area_pct = 100 * sum(patch_sizes) / n_total)
}

patch_summary <- do.call(rbind, lapply(names(sfe_list), function(nm) {
  do.call(rbind, lapply(score_cols, function(s) {
    count_patches(sfe_list[[nm]], nm, s, local_moran_df)
  }))
}))
print(patch_summary)
openxlsx::write.xlsx(patch_summary, "patch_summary.xlsx")


# spatial plots ----------------------------------------------------------
local_moran_df <- readRDS("local_moran_df.rds")
merged <- readRDS("merged_scored.rds")

plot_local_I <- function(merged, local_moran_df, score,
                         clip_low = "q02", clip_high = "q98") {
  sub <- local_moran_df[local_moran_df$signature == score, ]
  lookup <- setNames(sub$local_I, sub$spot)
  merged$local_I_plot <- lookup[colnames(merged)]
  
  lims <- quantile(merged$local_I_plot, probs = c(0.02, 0.98), na.rm = TRUE)
  bound <- max(abs(lims))
  
  SpatialFeaturePlot(merged, features = "local_I_plot", image.alpha = 0,
                     alpha = c(1,1), pt.size.factor = 3, stroke = 0.1,
                     min.cutoff = clip_low, max.cutoff = clip_high) &
    
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-bound, bound),
                         oob = scales::squish)
}

plot_local_I(merged, local_moran_df, "Inflammatory_UCell")

plot_local_L <- function(merged, local_lee_df, feature1, feature2,
                         clip_low = "q02", clip_high = "q98") {
  sub <- local_lee_df[local_lee_df$feature1 == feature1 &
                        local_lee_df$feature2 == feature2, ]
  lookup <- setNames(sub$local_L, sub$spot)
  merged$local_L_plot <- lookup[colnames(merged)]
  
  lims <- quantile(merged$local_L_plot, probs = c(0.02, 0.98), na.rm = TRUE)
  
  SpatialFeaturePlot(merged, features = "local_L_plot", image.alpha = 0,
                     alpha = c(1,1), pt.size.factor = 3, stroke = 0.1,
                     min.cutoff = clip_low, max.cutoff = clip_high) &
    
    # scale_fill_gradientn(values = RColorBrewer::brewer.pal(9, 'Reds'),
    #                      limits = c(lims[1], lims[2]))
  
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(lims[1], lims[2]),
                         oob = scales::squish)
}

plot_local_L(merged, local_lee_df, 
             "YAP_UCell", 
             "Fibrotic_UCell")

sfe_list        <- readRDS("spatial_corr_SFE.rds")
sfe_scores_list <- readRDS("sfe_scroes_list.rds")
score_cols <- c("YAP_UCell", "Inflammatory_UCell", "Fibrotic_UCell")
pairs <- list(
  c("YAP_UCell", "Inflammatory_UCell"),
  c("YAP_UCell", "Fibrotic_UCell"),
  c("Inflammatory_UCell", "Fibrotic_UCell")
)
BPPARAM <- BiocParallel::SnowParam(workers = parallel::detectCores() - 1)
# local Moran's I -----------------------------------------------------------
for (nm in names(sfe_list)) {
  sfe_list[[nm]] <- colDataUnivariate(
    sfe_list[[nm]], type = "localmoran_perm", features = score_cols,
    colGraphName = "visium", nsim = 999, zero.policy = TRUE,
    BPPARAM = BPPARAM
  )
}
saveRDS(sfe_list, "SFE_Moran_I_local.rds")
# check column names once before trusting the extraction below
# colnames(localResult(sfe_list[["WT"]], "localmoran_perm", "Inflammatory_UCell"))
local_moran_df <- do.call(rbind, lapply(names(sfe_list), function(nm) {
  do.call(rbind, lapply(score_cols, function(s) {
    lm_res <- as.data.frame(localResult(sfe_list[[nm]], "localmoran_perm", s))
    data.frame(
      genotype  = nm,
      signature = s,
      spot      = rownames(lm_res),
      local_I   = lm_res$Ii,
      p_value   = lm_res[["Pr(z != E(Ii)) Sim"]]
    )
  }))
}))
saveRDS(local_moran_df, "local_moran_df.rds")
# openxlsx::write.xlsx(local_moran_df, "local_moran_I.xlsx")
# local Lee's L ---------------------------------------------------------
local_lee_results <- list()
for (nm in names(sfe_scores_list)) {
  sfe <- sfe_scores_list[[nm]]
  for (pr in pairs) {
    f1 <- pr[1]; f2 <- pr[2]
    key <- paste(nm, f1, f2, sep = "__")
    
    ll <- calculateBivariate(
      sfe, type = "locallee", feature1 = f1, feature2 = f2,
      colGraphName = "visium", exprs_values = "scores",
      zero.policy = TRUE
    )
    
    ll_df <- as.data.frame(ll)
    
    local_lee_results[[key]] <- data.frame(
      genotype = nm, feature1 = f1, feature2 = f2,
      spot = colnames(sfe), local_L = ll_df[[1]]
    )
  }
}
local_lee_df <- do.call(rbind, local_lee_results)
rownames(local_lee_df) <- NULL
saveRDS(local_lee_df, "local_lee_df.rds")
# openxlsx::write.xlsx(local_lee_df, "local_lee_L.xlsx")

# hotspot fraction (permutation p, fixed raw threshold across all groups)
# permutation p-values don't assume normality of the underlying score
# (unlike the analytic column), and a fixed raw threshold avoids the
# FDR-scoping/adaptivity issues discussed earlier. uncorrected, so treat
# this as descriptive/exploratory -- at p<0.05 with no correction, expect
# ~5% of spots to appear "significant" by chance alone under a null.

p_threshold <- 0.05

local_moran_df$significant <- local_moran_df$p_value < p_threshold

hotspot_fraction <- local_moran_df %>%
  group_by(genotype, signature) %>%
  summarise(n_spots = n(),
            n_na = sum(is.na(significant)),
            n_significant = sum(significant, na.rm = TRUE),
            frac_hotspot = n_significant / (n_spots - n_na),
            .groups = "drop")

print(hotspot_fraction)
openxlsx::write.xlsx(hotspot_fraction, "hotspot_fraction.xlsx")

# patch count + area percentage -----------------------------------------

count_patches <- function(sfe, genotype_name, score, local_df, graph_name = "visium") {
  g  <- colGraph(sfe, graph_name)
  nb <- g$neighbours
  n_total <- ncol(sfe)
  
  sig_spots <- local_df$spot[local_df$genotype == genotype_name &
                               local_df$signature == score &
                               !is.na(local_df$significant) & local_df$significant]
  spot_idx <- match(sig_spots, colnames(sfe))
  spot_idx <- spot_idx[!is.na(spot_idx)]
  
  if (length(spot_idx) == 0) {
    return(data.frame(genotype = genotype_name, signature = score, n_patches = 0,
                      mean_patch_size = NA, max_patch_size = NA,
                      total_area_pct = 0))
  }
  
  visited <- setNames(rep(FALSE, length(spot_idx)), spot_idx)
  n_patches <- 0L
  patch_sizes <- c()
  
  for (idx in spot_idx) {
    if (visited[as.character(idx)]) next
    n_patches <- n_patches + 1L
    queue <- idx
    members <- c()
    while (length(queue) > 0) {
      node <- queue[1]; queue <- queue[-1]
      if (visited[as.character(node)]) next
      visited[as.character(node)] <- TRUE
      members <- c(members, node)
      neigh <- nb[[node]]
      neigh <- neigh[neigh %in% spot_idx]
      neigh <- neigh[!visited[as.character(neigh)]]
      queue <- c(queue, neigh)
    }
    patch_sizes <- c(patch_sizes, length(members))
  }
  
  data.frame(genotype = genotype_name, signature = score, n_patches = n_patches,
             mean_patch_size = mean(patch_sizes), max_patch_size = max(patch_sizes),
             total_area_pct = 100 * sum(patch_sizes) / n_total)
}

patch_summary <- do.call(rbind, lapply(names(sfe_list), function(nm) {
  do.call(rbind, lapply(score_cols, function(s) {
    count_patches(sfe_list[[nm]], nm, s, local_moran_df)
  }))
}))
print(patch_summary)
openxlsx::write.xlsx(patch_summary, "patch_summary.xlsx")

# spatial plot: local Lee's L only ---------------------------------------
# local Moran's I tracks the raw UCell score's own spatial pattern closely
# (it's essentially deviation-from-mean x neighbor-averaged-deviation), so
# it's largely redundant with the SpatialFeaturePlot(score) already run
# earlier in the pipeline. Lee's L is the one that adds new information,
# since it shows the joint local relationship between two scores.
local_lee_df <- readRDS("local_lee_df.rds")
merged <- readRDS("merged_scored.rds")

plot_local_L <- function(merged, local_lee_df, feature1, feature2,
                         clip_low = "q02", clip_high = "q98") {
  sub <- local_lee_df[local_lee_df$feature1 == feature1 &
                        local_lee_df$feature2 == feature2, ]
  lookup <- setNames(sub$local_L, sub$spot)
  merged$local_L_plot <- lookup[colnames(merged)]
  
  lims <- quantile(merged$local_L_plot, probs = c(0.02, 0.98), na.rm = TRUE)
  bound <- max(abs(lims))
  
  SpatialFeaturePlot(merged, features = "local_L_plot", image.alpha = 0,
                     alpha = c(1,1), pt.size.factor = 3, stroke = 0.1,
                     min.cutoff = clip_low, max.cutoff = clip_high) &
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-bound, bound),
                         oob = scales::squish)
}

plot_local_L(merged, local_lee_df, "YAP_UCell", "Fibrotic_UCell")

# patch summary dotplot ---------------------------------------------------

patch_summary$genotype <- factor(patch_summary$genotype, levels = c("WT", "mdx", "TG"))

ggplot(patch_summary, aes(x = genotype, y = signature)) +
  geom_point(aes(size = n_patches, color = mean_patch_size)) +
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(9, 'Reds')) +
  # scale_color_viridis_c(option = "magma", direction = -1) +
  scale_size_continuous(range = c(1, 3)) +
  labs(x = NULL, y = NULL, size = "n patches", color = "mean patch\nsize (spots)") +
  theme_minimal(base_size = 13) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.35, colour = 'black'),
    axis.ticks = element_line(linewidth = 0.35, colour = 'black'),
    axis.text = element_text(colour = 'black', size = 8),
    axis.title = element_text(colour = 'black', size = 8),
    title = element_text(colour = 'black', size = 8),
    legend.title = element_text(colour = 'black', size = 8),
    legend.text = element_text(colour = 'black', size = 8),
    legend.ticks = element_line(linewidth = 0.35, colour = 'black'),
    axis.ticks.length = unit(-0.1, 'cm')
  )

ggsave("Figures/spatial_score_size.svg",
       height = 3, width = 8, units = 'cm')
