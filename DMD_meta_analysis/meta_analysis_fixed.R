library(org.Hs.eg.db) 
library(HGNChelper)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(stringr)
library(readxl)
library(mashr)
library(ashr)

# load in datasets --------------------------------------------------------
datasets <- list(
  GSE6011      = read_xlsx("GSE6011_DEGs.xlsx"),
  GSE3307_FSH  = read_xlsx("GSE3307_FSH_set_DEGs.xlsx"),
  GSE3307_HUD  = read_xlsx("GSE3307_HUD_set_DEGs.xlsx"),
  GSE109178    = read_xlsx("GSE109178_DEGs.xlsx"),
  GSE86356     = read_xlsx("GSE86356_DEGs.xlsx"),
  sothers_2026 = read_xlsx("Sothers_2025_PMID_40641092.xlsx"),
  GSE38417     = read_xlsx("human microarray GSE38417.xlsx")
)

symbol_col <- c(GSE6011 = "Gene.symbol", GSE3307_FSH = "Gene.symbol", GSE3307_HUD = "Gene.symbol",
                GSE109178 = "Gene.symbol", GSE86356 = "Symbol", sothers_2026 = "hgnc_symbol",
                GSE38417 = "Gene.symbol")
logfc_col  <- c(GSE6011 = "logFC", GSE3307_FSH = "logFC", GSE3307_HUD = "logFC", GSE109178 = "logFC",
                GSE86356 = "log2FoldChange", sothers_2026 = "log2FoldChange", GSE38417 = "logFC")
stat_col   <- c(GSE6011 = "t", GSE3307_FSH = "t", GSE3307_HUD = "t", GSE109178 = "t",
                GSE86356 = "stat", sothers_2026 = "stat", GSE38417 = "t")
se_col     <- c(GSE6011 = NA, GSE3307_FSH = NA, GSE3307_HUD = NA, GSE109178 = NA,
                GSE86356 = "lfcSE", sothers_2026 = "lfcSE", GSE38417 = NA)
flip       <- c(GSE6011 = TRUE, GSE3307_FSH = FALSE, GSE3307_HUD = TRUE, GSE109178 = FALSE,
                GSE86356 = TRUE, sothers_2026 = FALSE, GSE38417 = TRUE)

# symbol mapping ----------------------------------------------------------
# previous symbols/aliases -> current HGNC
# LOC<entrez> -> symbol; ambiguous aliases kept as-is
hgnc_map <- getCurrentHumanMap()

map_symbols <- function(x) {
  u  <- unique(x)
  hg <- checkGeneSymbols(u, unmapped.as.na = FALSE, map = hgnc_map)
  new <- ifelse(!hg$Approved & !str_detect(hg$Suggested.Symbol, "///"), hg$Suggested.Symbol, u)
  loc <- str_match(new, "^LOC(\\d+)$")[, 2]
  if (any(!is.na(loc))) {
    eg  <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(na.omit(loc)),
                                 keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")
    hit <- !is.na(loc) & !is.na(eg[loc])
    new[hit] <- eg[loc[hit]]
  }
  setNames(new, u)[x]
}

# pre-process -------------------------------------------------------------
# make sure all gene symbols match
harmonize <- function(df, name) {
  df <- df %>%
    rename(symbol = all_of(symbol_col[[name]]),
           logFC  = all_of(logfc_col[[name]]),
           stat   = all_of(stat_col[[name]]))
  df$SE <- if (!is.na(se_col[[name]])) df[[se_col[[name]]]] else abs(df$logFC / df$stat)
  
  df <- df %>%
    mutate(symbol = trimws(sub("\\s*//+\\s*.*$", "", symbol))) %>%
    filter(!is.na(symbol), symbol != "", is.finite(logFC), is.finite(stat), is.finite(SE), SE > 0) %>%
    mutate(symbol = unname(map_symbols(symbol)))
  
  if (flip[[name]]) { df$logFC <- -df$logFC; df$stat <- -df$stat }
  
  # one row per gene (probes and merged aliases): keep strongest
  df %>%
    group_by(symbol) %>%
    slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(symbol, logFC, SE, stat)
}

harmonized <- imap(datasets, harmonize)
map_int(harmonized, nrow)

saveRDS(harmonized, "shared_mapping_all_data.rds")

# directional check: DMD down, rest up
dmd_up <- c("DMD", "SPP1", "POSTN", "COL1A1", "COL3A1", "COL1A2", "MYH3", "MYH8", "TNC")
imap_dfr(harmonized, ~ filter(.x, symbol %in% dmd_up) %>% mutate(dataset = .y)) %>%
  select(dataset, symbol, logFC) %>%
  pivot_wider(names_from = dataset, values_from = logFC) %>%
  print(width = Inf)

# DMD down in all, consistent sign

# matrices ----------------------------------------------------------------
to_wide <- function(col) {
  imap(harmonized, ~ select(.x, symbol, all_of(col)) %>% rename(!!.y := all_of(col))) %>%
    reduce(full_join, by = "symbol") %>%
    column_to_rownames("symbol") %>%
    as.matrix()
}
Bhat <- to_wide("logFC")
Shat <- to_wide("SE")[rownames(Bhat), colnames(Bhat)]

# fit on complete cases ---------------------------------------------------
complete_idx <- complete.cases(Bhat)
sum(complete_idx); nrow(Bhat)

data_fit <- mash_set_data(Bhat[complete_idx, ], Shat[complete_idx, ])
strong   <- get_significant_results(mash_1by1(data_fit), thresh = 0.05)
length(strong)

U_pca <- cov_pca(data_fit, npc = 5, subset = strong)
U_ed  <- cov_ed(data_fit, U_pca, subset = strong)
U_c   <- cov_canonical(data_fit)
m_fit <- mash(data_fit, c(U_c, U_ed))

# apply to all genes (missing = uninformative) ----------------------------
Bhat_full <- Bhat; Shat_full <- Shat
Shat_full[is.na(Bhat)] <- 1e6
Bhat_full[is.na(Bhat)] <- 0

data_full <- mash_set_data(Bhat_full, Shat_full)
m_all <- mash(data_full, g = get_fitted_g(m_fit), fixg = TRUE)

# per-dataset results -----------------------------------------------------
results <- list(pm = get_pm(m_all), psd = get_psd(m_all), lfsr = get_lfsr(m_all)) %>%
  imap(~ as_tibble(.x, rownames = "symbol") %>%
         pivot_longer(-symbol, names_to = "dataset", values_to = .y)) %>%
  reduce(left_join, by = c("symbol", "dataset")) %>%
  mutate(observed = !is.na(Bhat[cbind(symbol, dataset)])) %>%
  filter(observed) %>%
  select(-observed)

# overall effect across each gene's observed datasets ----------------------
# FE: mean of the observed datasets' true effects
# RE: mean effect allowing between-dataset variance tau^2
obs  <- !is.na(Bhat)
post <- mash_compute_posterior_matrices(m_all, data_full, output_posterior_cov = TRUE)
mu   <- post$PosteriorMean
Sig  <- post$PosteriorCov  # R x R x genes

re_stats <- t(vapply(seq_len(nrow(Bhat)), function(j) {
  o <- which(obs[j, ]); n <- length(o)
  m <- mu[j, o]
  S <- matrix(Sig[o, o, j], n, n)
  w <- rep(1 / n, n)
  var_fe <- drop(t(w) %*% S %*% w)
  # tau^2 = E[ sum_r (beta_r - mean)^2 / (n - 1) | data ]
  tau2 <- if (n > 1) {
    C <- diag(n) - 1 / n
    (sum(diag(C %*% S)) + drop(t(m) %*% C %*% m)) / (n - 1)
  } else NA_real_
  c(n, sum(w * m), var_fe, tau2)
}, numeric(4)))

# exact FE lfsr from the mash mixture: one A per missingness pattern
pattern <- apply(obs, 1, paste, collapse = "")
lfsr_fe <- numeric(nrow(Bhat))
for (p in unique(pattern)) {
  j <- which(pattern == p)
  A <- matrix(obs[j[1], ] / sum(obs[j[1], ]), nrow = 1)
  d <- mash_set_data(Bhat_full[j, , drop = FALSE], Shat_full[j, , drop = FALSE])
  lfsr_fe[j] <- mash_compute_posterior_matrices(m_all, d, A = A)$lfsr[, 1]
}

overall <- tibble(
  symbol  = rownames(Bhat),
  n_obs   = re_stats[, 1],
  pm_mean = re_stats[, 2],
  psd_fe  = sqrt(re_stats[, 3]),
  tau     = sqrt(re_stats[, 4]),
  psd_re  = sqrt(re_stats[, 3] + re_stats[, 4] / re_stats[, 1]),  # RE mean: var_fe + tau^2 / n
  psd_pred = sqrt(re_stats[, 3] + re_stats[, 4]),                # prediction (new dataset)
  lfsr_fe = lfsr_fe
) %>%
  # Knapp-Hartung: t with n-1 df; floored at exact FE (RE can't be more certain)
  mutate(df_re   = ifelse(n_obs > 1, n_obs - 1, NA),
         crit_re = qt(0.975, df_re),
         lfsr_re = pmax(pt(-abs(pm_mean) / psd_re, df_re), lfsr_fe),
         p_re    = pmin(1, 2 * lfsr_re)) %>%  # two-sided; p_re < 0.05 <=> 95% CI excludes 0
  arrange(lfsr_re)

# QC
sum(overall$lfsr_fe < 0.05, na.rm = TRUE)
sum(overall$p_re < 0.05, na.rm = TRUE)
table(overall$n_obs[overall$p_re < 0.05])

# save --------------------------------------------------------------------
results_wide <- results %>% pivot_wider(names_from = dataset, values_from = c(pm, psd, lfsr))

openxlsx::write.xlsx(results_wide, "prelim_res.xlsx")
openxlsx::write.xlsx(results, "long_res.xlsx")
openxlsx::write.xlsx(overall, "overall_res.xlsx")
saveRDS(m_all, "full_mash_model.rds")
saveRDS(m_fit, "subset_mash_fit.rds")

