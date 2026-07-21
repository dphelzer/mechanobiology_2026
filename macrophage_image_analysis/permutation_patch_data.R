library(dplyr)
library(readxl)
library(openxlsx)
library(purrr)

patch_data <- read_xlsx("patch_summary_results.xlsx")

patch_data$Group <- factor(
  patch_data$Group, 
  levels = c('WT', 'mdx', 'TG')
)

metrics <- c("mean_n_patches", "mean_patch_size", "pct_tissue_inflamed")

perm_results <- map(metrics, function(m) {
  f <- as.formula(paste(m, "~ Group"))
  perm_maxT(
    data = patch_data,
    formula = f,
    stat = "welch",
    resampling = "permutation"
  )
}) %>%
  set_names(metrics)

# optional but recommended: stamp each result with its metric name,
# useful if you ever combine these into one long table later
perm_results <- imap(perm_results, function(df, name) {
  df %>% mutate(Metric = name, .before = 1)
})

write.xlsx(x = perm_results, file = "permutation_results.xlsx")
