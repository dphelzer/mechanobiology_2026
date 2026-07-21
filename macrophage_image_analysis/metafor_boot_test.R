library(dplyr)
library(metafor)
library(emmeans)

results_dir <- "patch_analysis_results"
D_final <- readxl::read_excel("D_final.xlsx")

D_final <- D_final %>%
  mutate(vi = D_boot_se^2)

rma_model <- rma(
  yi = mean_D_obs,
  vi = vi,
  mods = ~ Group,
  data = D_final,
  method = "REML",
  test = "knha"
)

print(rma_model)

ref_grid <- emmprep(rma_model)
emm <- emmeans(ref_grid, ~ Group)

pairwise_results <- pairs(emm, adjust = "holm")
print(pairwise_results)

res <- as.data.frame(pairwise_results)

saveRDS(rma_model, file.path(results_dir, "D_rma_model.rds"))
saveRDS(pairwise_results, file.path(results_dir, "D_pairwise_results.rds"))

openxlsx::write.xlsx(res, "dispersion_res_boot.xlsx")