BiocManager::install("Voyager")
library(Seurat)
library(Voyager)
library(SingleCellExperiment)
library(tidyverse)
library(ggplot2)
library(spdep)

# make gene sets ----------------------------------------------------------
merged <- readRDS("merged_clusters.rds")

# YAP_targets <- readxl::read_xlsx('YAP target genes.xlsx', sheet = 2)
gene_sets <- list(
  
  YAP = c(
    'Ankrd1', 'Ccn1', 'Ccn2', 'Cdc42ep3', 'Parva', 'Nedd9',
    'Csrp1', 'Fmn2', 'Thbs1', 'Dcn', 'Sparc', 'Col6a2', 'Mfap5',
    'Ltbp1', 'Pxdn', 'Inhba', 'Tgfb2', 'Cdh13', 'Dab2', 'Myof',
    'Asap1', 'Ccdc80', 'Dock5', 'Arhgef17', 'Nuak2', 'Lats2', 
    'Amotl2'
    ),
  
  Inflammatory = c(
    'Lgals3', 'Fabp5', 'Spp1', 'Ctss',
    'Trem2', 'Thbs1', 'Gpnmb', 'Cd68',
    'Mmp19', 'Adgre1', 'Stat3', 'Lgals1', 'Il2rb',
    'Ifngr1', 'S100a6', 'Ptprc'
  ),
  
  Fibrotic = c(
    'Col1a1', 'Col1a2', 'Col3a1', 'Fn1', 
    'Tgfb1', 'Ccn2', 'Mmp3', 'Tnc', 'Bgn',
    'Timp1', 'Timp2', 'Pdgfra',
    'Pdgfrb', 'Postn', 'Cd44', 'Thbs4'
  )
)

merged <- AddModuleScore(object = merged,
                         features = gene_sets,
                         name = names(gene_sets))

merged <- UCell::AddModuleScore_UCell(
  merged,
  features = gene_sets
  )


p1 <- SpatialFeaturePlot(
  merged, features = "Fibrotic_UCell",
  image.alpha = 0,
  alpha = c(1,1),
  pt.size.factor = 3,
  stroke = 0.1
  ) &
  NoLegend()

p1_raster <- list()
for (i in seq(p1)) {
  p1_raster[[i]] <- ggrastr::rasterise(p1[[i]], 
                                       layers = "Spatial", 
                                       dpi = 800)
}

patchwork::wrap_plots(p1_raster, ncol = 3)
ggsave(file = "Figures/FibroticUCell_raster.svg",
       height = 10, width = 10,
       units = "cm")

merged$Genotype <- factor(merged$Genotype,
                          levels = c("WT", "mdx", "mdx_TG"))

## for addmodulescore
VlnPlot(merged,
        features = c("YAP1", "Inflammatory2", "Fibrotic3"),
        group.by = "Genotype",
        pt.size = 0) &
  geom_boxplot(fill="white", outliers = FALSE)

## for UCell Scores
VlnPlot(merged, features = "Fibrotic_UCell", 
        group.by = 'Genotype',
        pt.size = 0) +
  scale_fill_manual(values = c('black', 'red', 'blue')) +
  NoLegend() +
  labs(x = NULL,
       y = 'Score',
       title = NULL) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.35, colour = 'black'),
    axis.ticks = element_line(linewidth = 0.35, colour = 'black'),
    axis.text.y = element_text(colour = 'black', size = 8),
    axis.text.x = element_text(colour = 'black', size = 8, 
                               angle = 0, hjust = 0.5),
    axis.title = element_text(colour = 'black', size = 8),
    title = element_text(colour = 'black', size = 8),
    legend.title = element_text(colour = 'black', size = 8),
    legend.text = element_text(colour = 'black', size = 8),
    legend.ticks = element_line(linewidth = 0.35, colour = 'black'),
    axis.ticks.length.y = unit(-0.1, 'cm'),
    axis.ticks.length.x = unit(0, 'cm')
  )

ggsave("Figures/Fibrotic_score_vln.svg",
       height = 3.5, width = 4.5, units = 'cm')

FeatureScatter(merged, 'YAP_UCell', 'Fibrotic_UCell', split.by = 'Genotype')

RidgePlot(merged, features = c("YAP_UCell", "Inflammatory_UCell", "Fibrotic_UCell"),
          group.by = "Genotype")

merged <- UpdateSeuratObject(merged)
saveRDS(merged, "merged_scored.rds")

merged <- readRDS("merged_scored.rds")
# UCell separates signal so much cleaner


# run Moran's I -----------------------------------------------------------
merged[["Spatial"]]$scale.data <- NULL
merged[["SCT"]] <- NULL

# convert to SFE
WT  <- subset(merged, cells = Cells(merged@images[["slice1"]]))
mdx <- subset(merged, cells = Cells(merged@images[["slice1.2"]]))
TG  <- subset(merged, cells = Cells(merged@images[["slice1.3"]]))

sfe_WT  <- toSpatialFeatureExperiment(
  WT,  
  image_scalefactors = "lowres", 
  unit = "micron"
  )

sfe_mdx <- toSpatialFeatureExperiment(
  mdx,
  image_scalefactors = "lowres", 
  unit = "micron"
  )

sfe_TG  <- toSpatialFeatureExperiment(
  TG,  
  image_scalefactors = "lowres", 
  unit = "micron"
  )


score_cols <- c(
  "YAP_UCell",
  "Inflammatory_UCell",
  "Fibrotic_UCell"
)

sfe_list <- list(WT = sfe_WT, mdx = sfe_mdx, TG = sfe_TG)

for (nm in names(sfe_list)) {
  g <- findVisiumGraph(sfe_list[[nm]], zero.policy = TRUE)
  colGraph(sfe_list[[nm]], "visium") <- g
}

## windows
BPPARAM <- BiocParallel::SnowParam(workers = parallel::detectCores()-1)

## mac
# BPPARAM <- BiocParallel::MulticoreParam(parallel::detectCores()-1)

for (nm in names(sfe_list)) {
  sfe_list[[nm]] <- colDataMoransI(
    sfe_list[[nm]], features = score_cols, colGraphName = "visium",
    zero.policy = TRUE
  )
}

for (nm in names(sfe_list)) {
  sfe_list[[nm]] <- colDataUnivariate(
    sfe_list[[nm]], "moran.mc", features = score_cols,
    colGraphName = "visium", nsim = 99999, zero.policy = TRUE,
    BPPARAM = BPPARAM
  )
}

morans_results <- list()

for (nm in names(sfe_list)) {
  g <- colGraph(sfe_list[[nm]], "visium")
  
  res <- as.data.frame(colFeatureData(sfe_list[[nm]])[score_cols, ])
  colnames(res)[grep("^moran_", colnames(res))]                <- "moran"
  colnames(res)[grep("^K_", colnames(res))]                    <- "K"
  colnames(res)[grep("^moran\\.mc_statistic_", colnames(res))] <- "mc_statistic"
  colnames(res)[grep("^moran\\.mc_parameter_", colnames(res))] <- "mc_rank"
  colnames(res)[grep("^moran\\.mc_p\\.value_", colnames(res))] <- "mc_p_value"
  
  res$signature <- score_cols
  res$genotype  <- nm
  
  # analytic test, per signature
  analytic <- lapply(score_cols, function(s) {
    v <- colData(sfe_list[[nm]])[[s]]
    mt <- spdep::moran.test(v, g, zero.policy = TRUE)
    se <- sqrt(mt$estimate["Variance"])
    data.frame(
      analytic_I       = unname(mt$estimate["Moran I statistic"]),
      analytic_se      = unname(se),
      analytic_ci_low  = unname(mt$estimate["Moran I statistic"] - 1.96 * se),
      analytic_ci_high = unname(mt$estimate["Moran I statistic"] + 1.96 * se),
      analytic_p_value = mt$p.value
    )
  })
  analytic_df <- do.call(rbind, analytic)
  
  res <- cbind(res[, c("genotype", "signature", "moran", "K",
                       "mc_statistic", "mc_rank", "mc_p_value")],
               analytic_df)
  morans_results[[nm]] <- res
}


morans_df <- do.call(rbind, morans_results)
rownames(morans_df) <- NULL
print(morans_df)

openxlsx::write.xlsx(morans_df, "moran_univariate.xlsx")

morans_df$genotype <- factor(morans_df$genotype, 
                             levels = c("WT", "mdx", "TG"))

ggplot(morans_df, aes(x = signature, y = moran, fill = genotype)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(y = "Moran's I (observed)", x = NULL) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# two features ------------------------------------------------------------
score_cols <- c("YAP_UCell", "Inflammatory_UCell", "Fibrotic_UCell")

pairs <- list(
  c("YAP_UCell", "Inflammatory_UCell"),
  c("YAP_UCell", "Fibrotic_UCell"),
  c("Inflammatory_UCell", "Fibrotic_UCell")
)

# build the small "scores" assay SFE objects (as before), one per genotype
sfe_scores_list <- list()
for (nm in names(sfe_list)) {
  sfe <- sfe_list[[nm]]
  score_mat <- t(as.matrix(colData(sfe)[, score_cols]))
  rownames(score_mat) <- score_cols
  colnames(score_mat) <- colnames(sfe)
  
  sfe_scores <- SpatialFeatureExperiment(
    assays = list(scores = score_mat),
    colData = colData(sfe),
    spatialCoords = spatialCoords(sfe),
    spatialCoordsNames = NULL
  )
  colGraph(sfe_scores, "visium") <- colGraph(sfe, "visium")
  sfe_scores_list[[nm]] <- sfe_scores
}

saveRDS(sfe_scores_list, "sfe_scroes_list.rds")

# loop genotype x pair, store everything in one flat list
all_results <- list()

for (nm in names(sfe_scores_list)) {
  sfe <- sfe_scores_list[[nm]]
  
  for (pr in pairs) {
    f1 <- pr[1]
    f2 <- pr[2]
    key <- paste(nm, f1, f2, sep = "__")
    
    mc <- calculateBivariate(
      sfe, type = "lee.mc", feature1 = f1, feature2 = f2,
      colGraphName = "visium", exprs_values = "scores",
      zero.policy = TRUE, nsim = 99999,
      BPPARAM = BPPARAM
    )
    
    all_results[[key]] <- list(genotype = nm, feature1 = f1, feature2 = f2, mc = mc)
  }
}


summary_df <- do.call(rbind, lapply(names(all_results), function(key) {
  r <- all_results[[key]]
  mc <- r$mc
  
  data.frame(
    genotype  = r$genotype,
    feature1  = r$feature1,
    feature2  = r$feature2,
    lee_L     = as.numeric(mc$lee.mc_statistic),
    rank      = as.numeric(mc$lee.mc_parameter),
    p_value   = as.numeric(mc$lee.mc_p.value),
    nsim      = length(mc$lee.mc_res[[1]]),
    stringsAsFactors = FALSE
  )
}))

rownames(summary_df) <- NULL
print(summary_df)

all_results_analytic <- list()

for (nm in names(sfe_scores_list)) {
  sfe <- sfe_scores_list[[nm]]
  g <- colGraph(sfe, "visium")
  
  for (pr in pairs) {
    f1 <- pr[1]; f2 <- pr[2]
    key <- paste(nm, f1, f2, sep = "__")
    
    x <- assay(sfe, "scores")[f1, ]
    y <- assay(sfe, "scores")[f2, ]
    
    lt <- spdep::lee.test(x, y, g, zero.policy = TRUE)
    se <- unname(sqrt(lt$estimate["Variance"]))
    L  <- unname(lt$estimate["Lee's L statistic"])
    
    all_results_analytic[[key]] <- data.frame(
      genotype = nm, feature1 = f1, feature2 = f2,
      analytic_L = L,
      analytic_se = se,
      analytic_ci_low = L - 1.96 * se,
      analytic_ci_high = L + 1.96 * se,
      analytic_p_value = lt$p.value
    )
  }
}

analytic_bivariate_df <- do.call(rbind, all_results_analytic)
rownames(analytic_bivariate_df) <- NULL

# merge with existing permutation-based summary_df (from the lee.mc loop)
summary_df_full <- merge(
  summary_df, analytic_bivariate_df,
  by = c("genotype", "feature1", "feature2")
)
print(summary_df_full)

openxlsx::write.xlsx(summary_df_full, "bivariate_spatial_corr.xlsx")

saveRDS(sfe_list, "spatial_corr_SFE.rds")

# plotting ----------------------------------------------------------------
sfe_list <- readRDS("spatial_corr_SFE.rds")

theme_custom <- 
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.35, colour = 'black'),
    axis.ticks = element_line(linewidth = 0.35, colour = 'black'),
    axis.text = element_text(colour = 'black', size = 8),
    axis.title = element_text(colour = 'black', size = 8),
    title = element_text(colour = 'black', size = 8)
  )

## single feature ----
library(ggrastr)
library(svglite)
moran_scatter_gg <- function(
    sfe, 
    feature, 
    graph_name = "visium",
    colData_source = TRUE, 
    assay_name = "scores",
    zero.policy = TRUE, 
    label = NULL
    ) {
  
  g <- colGraph(sfe, graph_name)
  v <- if (colData_source) colData(sfe)[[feature]] else assay(sfe, assay_name)[feature, ]
  v_lag <- lag.listw(g, v, zero.policy = zero.policy)
  df <- data.frame(value = v, lag = v_lag)
  fit <- lm(lag ~ value, data = df)
  slope <- coef(fit)[2]
  
  resid_std <- (residuals(fit) - mean(residuals(fit))) / sd(residuals(fit))
  infl <- which(abs(resid_std) > 2)
  df$influential <- FALSE
  df$influential[infl] <- TRUE
  
  title_txt <- if (is.null(label)) feature else label
  
  ggplot(df, aes(x = value, y = lag)) +
    
    geom_point_rast(aes(color = influential), 
                    shape = 16, alpha = 0.6, size = 0.8,
                    raster.dpi = 1200) +
    
    # geom_point(aes(color = influential), 
    #            shape = 16, alpha = 0.6, size = 0.8) +
    
    geom_smooth(method = "lm", 
                se = TRUE, 
                color = "black", 
                linewidth = 0.35) +
    
    geom_vline(xintercept = mean(v), 
               linetype = "dashed", 
               color = "grey50",
               linewidth = 0.35) +
    
    geom_hline(yintercept = mean(v_lag), 
               linetype = "dashed", 
               color = "grey50",
               linewidth = 0.35) +
    
    scale_color_manual(values = c(`FALSE` = "steelblue", `TRUE` = "firebrick")) +
    
    labs(title = title_txt, 
         # subtitle = paste0("Moran's I = ", round(slope, 4)),
         x = paste(feature, "(score)"), 
         y = paste(feature, "(lag)"),
         color = "influential") +
    theme_custom +
    theme(legend.position = 'none',
          aspect.ratio = 1)
}

moran_scatter_gg(sfe_list[["WT"]], "YAP_UCell", label = "WT: YAP score")

## bivariate Lee's L ----
lee_scatter_gg <- function(
    sfe, 
    feature1, 
    feature2,
    graph_name = "visium",
    assay_name = "scores",
    zero.policy = TRUE,
    label = NULL
    ) {
  
  g <- colGraph(sfe, graph_name)
  
  x <- assay(sfe, assay_name)[feature1, ]
  y <- assay(sfe, assay_name)[feature2, ]
  
  x_lag <- lag.listw(g, x, zero.policy = zero.policy)
  y_lag <- lag.listw(g, y, zero.policy = zero.policy)
  
  df <- data.frame(x_lag = x_lag, y_lag = y_lag)
  fit <- lm(y_lag ~ x_lag, data = df)
  
  resid_std <- (residuals(fit) - mean(residuals(fit))) / sd(residuals(fit))
  infl <- which(abs(resid_std) > 2)
  df$influential <- FALSE
  df$influential[infl] <- TRUE
  
  title_txt <- if (is.null(label)) paste(feature1, "vs", feature2) else label
  
  ggplot(df, aes(x = x_lag, y = y_lag)) +
    
    geom_point_rast(aes(color = influential), 
                    shape = 16, alpha = 0.6, size = 0.8,
                    raster.dpi = 1200) +
    
    # geom_point(aes(color = influential), 
    #            shape = 16,
    #            alpha = 0.6, size = 0.8) +
    
    geom_smooth(method = "lm", 
                se = TRUE, 
                color = "black", 
                linewidth = 0.35) +
    
    scale_color_manual(values = c(`FALSE` = "steelblue", `TRUE` = "firebrick")) +
    
    labs(
      title = title_txt,
      x = paste(feature1, "(spatial lag)"),
      y = paste(feature2, "(spatial lag)"),
      color = "influential"
    ) +
    theme_custom +
    theme(legend.position = "none",
          aspect.ratio = 1)
}

lee_scatter_gg(sfe_scores_list[["mdx"]], "Inflammatory_UCell", "Fibrotic_UCell", label = "mdx")

plot_list <- list()
for (i in seq(sfe_list)) {
  plot_list[[i]] <- moran_scatter_gg(sfe_list[[i]], 
                   "YAP_UCell", 
                   label = paste(names(sfe_list[i]), 
                                 ": YAP score"))
}

patchwork::wrap_plots(plot_list, ncol = 3)
ggsave(filename = "Figures/YAP_score_scatter_updated.svg",
       height = 5.2, width = 5.2*3, units = "cm")

plot_list <- list()
for (i in seq(sfe_scores_list)) {
  plot_list[[i]] <- lee_scatter_gg(
    sfe_scores_list[[i]], 
    "YAP_UCell", "Fibrotic_UCell",
    label = names(sfe_scores_list[i]))
}

patchwork::wrap_plots(plot_list, ncol = 3)
ggsave(filename = "Figures/YAP_Fibrotic_scatter_updated.svg",
       height = 5.2, width = 5.2*3, units = "cm")

