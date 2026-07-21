library(tidyverse)

genotype_pairs <- list(c("WT", "mdx"), c("WT", "TG"), c("mdx", "TG"))

welch_t_jackknife <- function(est1, se1, k1, est2, se2, k2) {
  diff <- est1 - est2
  se_diff <- sqrt(se1^2 + se2^2)
  t_stat <- diff / se_diff
  df <- (se1^2 + se2^2)^2 /
    ((se1^4 / (k1 - 1)) + (se2^4 / (k2 - 1)))
  p_value <- 2 * pt(-abs(t_stat), df = df)
  data.frame(diff = diff, se_diff = se_diff, stat = t_stat, df = df,
             p_value = p_value, method = "Welch-t")
}

z_test_raw_se <- function(est1, se1, est2, se2) {
  diff <- est1 - est2
  se_diff <- sqrt(se1^2 + se2^2)
  z_stat <- diff / se_diff
  p_value <- 2 * (1 - pnorm(abs(z_stat)))
  data.frame(diff = diff, se_diff = se_diff, stat = z_stat, df = NA,
             p_value = p_value, method = "z-test")
}

## Moran's I ----
moran_jk_df <- readxl::read_excel("moran_jackknife_score_level.xlsx")
signatures <- unique(moran_jk_df$signature)
all_comparisons_I <- list()

for (gp in genotype_pairs) {
  g1 <- gp[1]; g2 <- gp[2]
  for (s in signatures) {
    row1 <- moran_jk_df[moran_jk_df$genotype == g1 & moran_jk_df$signature == s, ]
    row2 <- moran_jk_df[moran_jk_df$genotype == g2 & moran_jk_df$signature == s, ]
    if (nrow(row1) == 0 || nrow(row2) == 0) next
    
    res_t <- welch_t_jackknife(row1$mean_I, row1$se_jackknife, row1$k_genes,
                               row2$mean_I, row2$se_jackknife, row2$k_genes)
    res_z <- z_test_raw_se(row1$mean_I, row1$se_jackknife, row2$mean_I, row2$se_jackknife)
    
    res_t$genotype1 <- g1; res_t$genotype2 <- g2; res_t$signature <- s
    res_z$genotype1 <- g1; res_z$genotype2 <- g2; res_z$signature <- s
    
    all_comparisons_I[[paste(g1, g2, s, "t", sep = "__")]] <- res_t
    all_comparisons_I[[paste(g1, g2, s, "z", sep = "__")]] <- res_z
  }
}

comparison_df_I <- do.call(rbind, all_comparisons_I)
rownames(comparison_df_I) <- NULL
comparison_df_I <- comparison_df_I[, c("genotype1", "genotype2", "signature", "method",
                                       "diff", "se_diff", "stat", "df", "p_value")]

comparison_df_I <- comparison_df_I %>%
  group_by(method, signature) %>%
  mutate(padj = p.adjust(p_value, method = "holm")) %>%
  ungroup()

print(comparison_df_I)
openxlsx::write.xlsx(comparison_df_I, "Moran_I_sig_test_score_jackknife.xlsx")

## Lee's L ----
lee_jk_df <- readxl::read_excel("lee_jackknife_score_level.xlsx")

pairs_lee <- list(
  c("YAP", "Inflammatory"),
  c("YAP", "Fibrotic"),
  c("Inflammatory", "Fibrotic")
)

all_comparisons_L <- list()

for (gp in genotype_pairs) {
  g1 <- gp[1]; g2 <- gp[2]
  for (pr in pairs_lee) {
    s1 <- pr[1]; s2 <- pr[2]
    row1 <- lee_jk_df[lee_jk_df$genotype == g1 & lee_jk_df$set1 == s1 & lee_jk_df$set2 == s2, ]
    row2 <- lee_jk_df[lee_jk_df$genotype == g2 & lee_jk_df$set1 == s1 & lee_jk_df$set2 == s2, ]
    if (nrow(row1) == 0 || nrow(row2) == 0) next
    
    res_t <- welch_t_jackknife(row1$mean_L, row1$se_jackknife, row1$k_genes,
                               row2$mean_L, row2$se_jackknife, row2$k_genes)
    res_z <- z_test_raw_se(row1$mean_L, row1$se_jackknife, row2$mean_L, row2$se_jackknife)
    
    res_t$genotype1 <- g1; res_t$genotype2 <- g2; res_t$set1 <- s1; res_t$set2 <- s2
    res_z$genotype1 <- g1; res_z$genotype2 <- g2; res_z$set1 <- s1; res_z$set2 <- s2
    
    all_comparisons_L[[paste(g1, g2, s1, s2, "t", sep = "__")]] <- res_t
    all_comparisons_L[[paste(g1, g2, s1, s2, "z", sep = "__")]] <- res_z
  }
}

comparison_df_L <- do.call(rbind, all_comparisons_L)
rownames(comparison_df_L) <- NULL
comparison_df_L <- comparison_df_L[, c("genotype1", "genotype2", "set1", "set2", "method",
                                       "diff", "se_diff", "stat", "df", "p_value")]

comparison_df_L <- comparison_df_L %>%
  mutate(score_contrast = paste(set1, set2, sep = "_")) %>%
  group_by(method, score_contrast) %>%
  mutate(padj = p.adjust(p_value, method = "holm")) %>%
  ungroup()

print(comparison_df_L)
openxlsx::write.xlsx(comparison_df_L, "Lee_L_sig_test_score_jackknife.xlsx")
