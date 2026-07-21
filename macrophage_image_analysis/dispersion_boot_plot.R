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


# Group-level estimate + 95% CI, with individual samples overlaid -----------

library(ggplot2)

genotype_colors <- c("WT" = "grey10", "mdx" = "red", "TG" = "blue")

base_theme <- theme_minimal(base_size = 8) +
  theme(
    text = element_text(color = "black", size = 8),
    axis.text = element_text(color = "black", size = 8),
    axis.title = element_text(color = "black", size = 8),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.ticks.length.y = unit(-0.1, 'cm'),
    axis.ticks.length.x = unit(0, "cm"), 
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )

group_est <- as.data.frame(emm)
group_est$Group <- factor(group_est$Group, levels = c("WT", "mdx", "TG"))

D_final$Group <- factor(D_final$Group, levels = c("WT", "mdx", "TG"))

ggplot() +
  geom_segment(
    data = group_est,
    aes(
      y = lower.CL, 
      yend = upper.CL, 
      x = Group, 
      xend = Group
      ),
    linewidth = 0.35
  ) +
  geom_point(
    data = group_est,
    aes(y = emmean, x = Group),
    size = 3,
    shape = 21,
    fill = "white"
  ) +
  geom_jitter(
    data = D_final,
    aes(y = mean_D_obs, x = Group, color = Group),
    width = 0.1,
    alpha = 0.6,
    size = 2,
    shape = 16
  ) +
  scale_color_manual(values = genotype_colors) +
  labs(y = "Index of dispersion", x = NULL) +
  base_theme +
  theme(legend.position = "none")

ggsave(filename = "Figures/dispersion_CI_dotplot.svg",
       height = 3, width = 3, units = "cm")

