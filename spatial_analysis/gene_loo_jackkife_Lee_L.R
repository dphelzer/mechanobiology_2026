library(UCell)
library(spdep)
library(SpatialFeatureExperiment)
library(tidyverse)

sfe_list  <- readRDS("spatial_corr_SFE.rds")
ranks_list <- readRDS("ucell_ranks_by_genotype.rds")

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

pairs_lee <- list(
  c("YAP", "Inflammatory"),
  c("YAP", "Fibrotic"),
  c("Inflammatory", "Fibrotic")
)

lodo_lee_score <- function(sfe, ranks, genes1, genes2, graph_name = "visium", zero.policy = TRUE) {
  g <- colGraph(sfe, graph_name)
  n <- length(g$neighbours)
  S2 <- sum(sapply(g$weights, function(w) sum(w)^2))
  k1 <- length(genes1)
  k2 <- length(genes2)
  
  full1 <- UCell::ScoreSignatures_UCell(features = list(sig = genes1), precalc.ranks = ranks)[, 1]
  full2 <- UCell::ScoreSignatures_UCell(features = list(sig = genes2), precalc.ranks = ranks)[, 1]
  theta_full <- spdep::lee(full1, full2, g, n, S2, zero.policy = zero.policy)$L
  
  loo1 <- vapply(seq_len(k1), function(i) {
    s1 <- UCell::ScoreSignatures_UCell(features = list(sig = genes1[-i]), precalc.ranks = ranks)[, 1]
    spdep::lee(s1, full2, g, n, S2, zero.policy = zero.policy)$L
  }, numeric(1))
  
  loo2 <- vapply(seq_len(k2), function(i) {
    s2 <- UCell::ScoreSignatures_UCell(features = list(sig = genes2[-i]), precalc.ranks = ranks)[, 1]
    spdep::lee(full1, s2, g, n, S2, zero.policy = zero.policy)$L
  }, numeric(1))
  
  loo_all <- c(loo1, loo2)
  genes_all <- c(genes1, genes2)
  k <- length(loo_all)
  
  se <- sqrt((k - 1) / k * sum((loo_all - mean(loo_all))^2))
  
  infl <- abs(loo_all - mean(loo_all))
  top2 <- order(infl, decreasing = TRUE)[1:min(2, k)]
  
  data.frame(
    mean_L = theta_full,
    se_jackknife = se,
    k_genes = k,
    top_influential_genes = paste(genes_all[top2], collapse = ", "),
    top_influential_deltas = paste(round(infl[top2], 4), collapse = ", ")
  )
}

lee_jk_results <- list()

for (nm in names(sfe_list)) {
  for (pr in pairs_lee) {
    s1 <- pr[1]; s2 <- pr[2]
    t0 <- Sys.time()
    res <- lodo_lee_score(sfe_list[[nm]], ranks_list[[nm]], gene_sets[[s1]], gene_sets[[s2]])
    cat(nm, s1, s2, ":", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "sec\n")
    
    res$genotype <- nm
    res$set1 <- s1
    res$set2 <- s2
    lee_jk_results[[paste(nm, s1, s2, sep = "__")]] <- res
  }
}

lee_jk_df <- do.call(rbind, lee_jk_results)
rownames(lee_jk_df) <- NULL
lee_jk_df <- lee_jk_df[, c("genotype", "set1", "set2", "mean_L", "se_jackknife", "k_genes",
                           "top_influential_genes", "top_influential_deltas")]
print(lee_jk_df)

openxlsx::write.xlsx(lee_jk_df, "lee_jackknife_score_level.xlsx")
