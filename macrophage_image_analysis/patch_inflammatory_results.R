library(dplyr)
library(purrr)
library(stringr)
library(spatstat)
library(ggplot2)
library(patchwork)


# Load what's needed --------------------------------------------------

results_dir <- "patch_analysis_results"
run_suffix <- "full"   # change to "test3" if working from the test run

tissue_windows <- readRDS(file.path(results_dir, "tissue_windows.rds"))
patch_sweep <- readRDS(file.path(results_dir, paste0("patch_sweep_", run_suffix, ".rds")))

sample_ids_run <- unique(patch_sweep$Sample)


# Patch summary ---------------------------------------------------------

tissue_area_by_sample <- map_dfr(sample_ids_run, function(s) {
  tibble(
    Sample = s,
    tissue_area = spatstat.geom::area(tissue_windows[[s]])
  )
})

patch_summary <- patch_sweep %>%
  group_by(Sample) %>%
  summarize(
    mean_n_patches = mean(n_patches),
    mean_patch_size = mean(mean_patch_size),
    mean_total_inflamed_area = mean(total_inflamed_area),
    .groups = "drop"
  ) %>%
  left_join(tissue_area_by_sample, by = "Sample") %>%
  mutate(pct_tissue_inflamed = 100 * mean_total_inflamed_area / tissue_area)

print(patch_summary)
saveRDS(patch_summary, file.path(results_dir, paste0("patch_summary_", run_suffix, ".rds")))


# Genotype labeling by sample NAME, not row position ---------------------

patch_summary <- patch_summary %>%
  mutate(Group = case_when(
    str_starts(Sample, "WT")   ~ "WT",
    str_starts(Sample, "MDX")  ~ "mdx",
    str_starts(Sample, "SSPN") ~ "TG",
    TRUE ~ NA_character_
  ))

stopifnot(!anyNA(patch_summary$Group))

patch_summary$Group <- factor(
  patch_summary$Group,
  levels = c("WT", "mdx", "TG")
)


# Plots -------------------------------------------------------------------

p_group_npatches <- ggplot(patch_summary, aes(x = Group, y = mean_n_patches)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  labs(y = "Mean n_patches (50-100um)") +
  theme_minimal(base_size = 8)

p_group_patchsize <- ggplot(patch_summary, aes(x = Group, y = mean_patch_size)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  labs(y = "Mean patch size (um^2)") +
  theme_minimal(base_size = 8)

p_group_pctinflamed <- ggplot(patch_summary, aes(x = Group, y = pct_tissue_inflamed)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  labs(y = "% tissue inflamed") +
  theme_minimal(base_size = 8)

p_group_npatches + p_group_patchsize + p_group_pctinflamed

openxlsx::write.xlsx(patch_summary, "patch_summary_results.xlsx")
