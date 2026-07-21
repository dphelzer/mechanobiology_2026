loo_gene_contribution <- function(
    seurat_obj, yap_genes, fibrotic_col = "Fibrotic",
                                  cell_subset = NULL, 
    method = "spearman") {
  
  obj <- seurat_obj
  if (!is.null(cell_subset)) obj <- subset(obj, cells = cell_subset)
  
  fibrotic_score <- FetchData(obj, vars = fibrotic_col)[[1]]
  
  obj <- UCell::AddModuleScore_UCell(obj, features = list(YAP_full = yap_genes), name = "_full")
  baseline_score <- FetchData(obj, vars = "YAP_full_full")[[1]]
  baseline_cor <- cor(baseline_score, fibrotic_score, method = method, use = "complete.obs")
  
  loo_results <- lapply(yap_genes, function(g) {
    genes_sub <- setdiff(yap_genes, g)
    obj_tmp <- UCell::AddModuleScore_UCell(obj, features = list(YAP_loo = genes_sub), name = "_loo")
    score_loo <- FetchData(obj_tmp, vars = "YAP_loo_loo")[[1]]
    cor_loo <- cor(score_loo, fibrotic_score, method = method, use = "complete.obs")
    data.frame(gene_removed = g, cor_without = cor_loo,
               delta = baseline_cor - cor_loo)
  })
  
  loo_df <- do.call(rbind, loo_results)
  loo_df <- loo_df[order(-loo_df$delta), ]
  loo_df$gene_removed <- factor(loo_df$gene_removed, levels = rev(loo_df$gene_removed))
  
  p <- ggplot(loo_df, aes(x = gene_removed, y = delta)) +
    geom_col(fill = "steelblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    coord_flip() +
    labs(x = NULL, y = "drop in correlation when gene removed\n(baseline_cor - cor_without)",
         title = paste0("baseline correlation (full YAP set) = ", round(baseline_cor, 3))) +
    theme_minimal(base_size = 12)
  
  list(plot = p, data = loo_df, baseline_cor = baseline_cor)
}

yap_genes <- gene_sets$YAP
res <- loo_gene_contribution(
  FAPs, yap_genes, fibrotic_col = "Fibrotic")
res$plot
res$data
saveRDS(res, 'res_jackknife_loo_cor.rds')
