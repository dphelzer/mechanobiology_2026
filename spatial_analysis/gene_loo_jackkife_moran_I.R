library(UCell)
library(spdep)
library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(tidyverse)

sfe_list  <- readRDS("spatial_corr_SFE.rds")
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

lodo_moran_score <- function(sfe, ranks, genes, graph_name = "visium", zero.policy = TRUE) {
  g <- colGraph(sfe, graph_name)
  k <- length(genes)
  
  full_score <- UCell::ScoreSignatures_UCell(features = list(sig = genes), precalc.ranks = ranks)[, 1]
  theta_full <- unname(spdep::moran.test(full_score, g, zero.policy = zero.policy)$estimate["Moran I statistic"])
  
  theta_loo <- vapply(seq_len(k), function(i) {
    score <- UCell::ScoreSignatures_UCell(features = list(sig = genes[-i]), precalc.ranks = ranks)[, 1]
    unname(spdep::moran.test(score, g, zero.policy = zero.policy)$estimate["Moran I statistic"])
  }, numeric(1))
  
  se <- sqrt((k - 1) / k * sum((theta_loo - mean(theta_loo))^2))
  
  infl <- abs(theta_loo - mean(theta_loo))
  top2 <- order(infl, decreasing = TRUE)[1:min(2, k)]
  
  data.frame(
    mean_I = theta_full,
    se_jackknife = se,
    k_genes = k,
    top_influential_genes = paste(genes[top2], collapse = ", "),
    top_influential_deltas = paste(round(infl[top2], 4), collapse = ", ")
  )
}

ranks_list <- list()
for (nm in names(sfe_list)) {
  t0 <- Sys.time()
  ranks_list[[nm]] <- UCell::StoreRankings_UCell(assay(sfe_list[[nm]], "logcounts"))
  cat(nm, "StoreRankings_UCell:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "sec\n")
}

moran_jk_results <- list()

for (nm in names(sfe_list)) {
  for (set_name in names(gene_sets)) {
    t0 <- Sys.time()
    res <- lodo_moran_score(sfe_list[[nm]], ranks_list[[nm]], gene_sets[[set_name]])
    cat(nm, set_name, ":", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "sec\n")
    
    res$genotype  <- nm
    res$signature <- set_name
    moran_jk_results[[paste(nm, set_name, sep = "__")]] <- res
  }
}

moran_jk_df <- do.call(rbind, moran_jk_results)
rownames(moran_jk_df) <- NULL
moran_jk_df <- moran_jk_df[, c("genotype", "signature", "mean_I", "se_jackknife", "k_genes",
                               "top_influential_genes", "top_influential_deltas")]
print(moran_jk_df)

openxlsx::write.xlsx(moran_jk_df, "moran_jackknife_score_level.xlsx")
saveRDS(ranks_list, "ucell_ranks_by_genotype.rds")
