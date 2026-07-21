library(tidyverse)
library(permuco)
library(patchwork)
library(ggridges)

# CI on disperion ---------------------------------------------------------

CD68 <- readxl::read_excel(
  "Sarcospan Project CD68 Macrophage Densities.xlsx",
  sheet = 14
  )

CD68$Genotype <- factor(
  CD68$Genotype, 
  levels = c("WT", "mdx", "TG")
  )

## global D
area <- 0.39811499
CD68 <- CD68 %>% 
  mutate(Count = Value*area) %>% 
  dplyr::group_by(Sample) %>% 
  mutate(D = var(Count)/mean(Count)) %>% 
  ungroup()

CD68$Count <- round(CD68$Count, 0)

dis_dat <- CD68 %>% 
  group_by(Genotype, Sample) %>% 
  summarize(Index_of_disperion = mean(D)) %>% 
  select(Genotype, Index_of_disperion)

openxlsx::write.xlsx(
  dis_dat,
  "index_of_dispersion.xlsx"
  )

## counts
fit_dispersion <- function(df_geno) {
  fit <- glm(Count ~ factor(Sample), family = poisson(), data = df_geno)
  phi <- sum(residuals(fit, type = "pearson")^2) / df.residual(fit)
  list(phi = phi, df = df.residual(fit))
}

WT_fit  <- fit_dispersion(CD68 %>% filter(Genotype == "WT"))
mdx_fit <- fit_dispersion(CD68 %>% filter(Genotype == "mdx"))
TG_fit  <- fit_dispersion(CD68 %>% filter(Genotype == "TG"))

WT_fit$phi; mdx_fit$phi; TG_fit$phi

sample_counts <- split(CD68$Count, CD68$Sample)

boot_phi <- function(df_geno, n_boot = 10000) {
  animals <- unique(df_geno$Sample)
  replicate(n_boot, {
    boot_data <- map_dfr(animals, function(s) {
      vals <- sample_counts[[s]]
      tibble(Sample = s, Count = sample(vals, length(vals), replace = TRUE))
    })
    fit <- glm(Count ~ factor(Sample), family = poisson(), data = boot_data)
    sum(residuals(fit, type = "pearson")^2) / df.residual(fit)
  })
}

boot_phi_WT <- boot_phi(CD68 %>% filter(Genotype == "WT"))
boot_phi_mdx <- boot_phi(CD68 %>% filter(Genotype == "mdx"))
boot_phi_TG  <- boot_phi(CD68 %>% filter(Genotype == "TG"))

# permutation test --------------------------------------------------------
sample_counts <- split(CD68$Count, CD68$Sample)
animal_geno <- CD68 %>% distinct(Sample, Genotype)

## Closed-form phi: Poisson GLM MLE with one intercept per animal
## reproduces each animal's own sample mean exactly, so this is
## mathematically identical to glm(Count ~ factor(Sample), poisson())
## but avoids IRLS overhead across ~35,000 evaluations.
fit_phi_fast <- function(animal_ids) {
  ss <- 0; n_total <- 0
  for (s in animal_ids) {
    vals <- sample_counts[[s]]
    m <- mean(vals)
    ss <- ss + sum((vals - m)^2 / m)
    n_total <- n_total + length(vals)
  }
  ss / (n_total - length(animal_ids))
}

WT_animals  <- animal_geno$Sample[animal_geno$Genotype == "WT"]
mdx_animals <- animal_geno$Sample[animal_geno$Genotype == "mdx"]
TG_animals  <- animal_geno$Sample[animal_geno$Genotype == "TG"]
all_animals <- c(WT_animals, mdx_animals, TG_animals)
n_WT <- length(WT_animals); n_mdx <- length(mdx_animals); n_TG <- length(TG_animals)

## observed statistics
obs <- c(
  mdx_WT = log(fit_phi_fast(mdx_animals) / fit_phi_fast(WT_animals)),
  mdx_TG = log(fit_phi_fast(mdx_animals) / fit_phi_fast(TG_animals)),
  TG_WT  = log(fit_phi_fast(TG_animals)  / fit_phi_fast(WT_animals))
)

## exact enumeration of all joint 3-group label assignments
WT_sets <- combn(all_animals, n_WT, simplify = FALSE)
all_perms <- list(); idx <- 1
for (wt_set in WT_sets) {
  remaining <- setdiff(all_animals, wt_set)
  mdx_sets <- combn(remaining, n_mdx, simplify = FALSE)
  for (mdx_set in mdx_sets) {
    tg_set <- setdiff(remaining, mdx_set)
    all_perms[[idx]] <- list(WT = wt_set, mdx = mdx_set, TG = tg_set)
    idx <- idx + 1
  }
}
n_perm <- length(all_perms)
n_perm  # confirm 11550

null_mat <- matrix(NA_real_, nrow = n_perm, ncol = 3,
                   dimnames = list(NULL, c("mdx_WT", "mdx_TG", "TG_WT")))
for (i in seq_len(n_perm)) {
  p <- all_perms[[i]]
  phi_WT  <- fit_phi_fast(p$WT)
  phi_mdx <- fit_phi_fast(p$mdx)
  phi_TG  <- fit_phi_fast(p$TG)
  null_mat[i, "mdx_WT"] <- log(phi_mdx / phi_WT)
  null_mat[i, "mdx_TG"] <- log(phi_mdx / phi_TG)
  null_mat[i, "TG_WT"]  <- log(phi_TG  / phi_WT)
}

## use |log-ratio| for a two-sided test, since direction for mdx vs TG isn't assumed a priori
full_mat <- rbind(abs(obs), abs(null_mat))
result <- compute_stepdownmaxT(full_mat, alternative = "greater")
result


# bootstrap phi plot ------------------------------------------------------
phi_df <- bind_rows(
  data.frame(Genotype = "WT",  phi = boot_phi_WT),
  data.frame(Genotype = "mdx", phi = boot_phi_mdx),
  data.frame(Genotype = "TG",  phi = boot_phi_TG)
) %>% mutate(Genotype = factor(Genotype, levels = c("WT", "mdx", "TG")))

saveRDS(phi_df, "dispersion_bootstrap.rds")

phi_df <- readRDS("dispersion_bootstrap.rds")

summary_df <- phi_df %>%
  group_by(Genotype) %>%
  summarize(
    median = median(phi),
    q025 = quantile(phi, 0.025),
    q25  = quantile(phi, 0.25),
    q75  = quantile(phi, 0.75),
    q975 = quantile(phi, 0.975)
  )

genotype_colors <- c("WT" = "grey10", "mdx" = "red", "TG" = "blue")

base_theme <- theme_minimal(base_size = 8) +
  theme(
    text = element_text(color = "black", size = 8),
    axis.text = element_text(color = "black", size = 8),
    axis.title = element_text(color = "black", size = 8),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.ticks.length = unit(-0.1,'cm'),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )

## shared x-range across both panels, computed once from the actual data
x_range <- range(c(phi_df$phi, summary_df$q025, summary_df$q975))

phi_df$Genotype <- factor(phi_df$Genotype, levels = c("TG", "mdx", "WT"))     # bottom-to-top order
summary_df$Genotype <- factor(summary_df$Genotype, levels = c("TG", "mdx", "WT"))

p_ridge <- ggplot(phi_df, 
                  aes(x = phi, y = Genotype, fill = Genotype, 
                      color = Genotype)) +
  
  geom_density_ridges(alpha = 0.4, linewidth = 0.35, 
                      scale = 1.5, rel_min_height = 0.001) +
  
  geom_vline(xintercept = 1, 
             linetype = "dashed", 
             color = "black", 
             linewidth = 0.35) +
  
  scale_fill_manual(values = genotype_colors) +
  scale_color_manual(values = genotype_colors) +
  scale_x_continuous(limits = x_range) +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(legend.position = "none")

p_box <- ggplot(summary_df, aes(y = Genotype, color = Genotype)) +
  geom_segment(aes(x = q025, xend = q975, yend = Genotype), linewidth = 0.35) +
  geom_crossbar(aes(xmin = q25, x = median, xmax = q75), width = 0.4, fill = "white", linewidth = 0.35) +
  scale_color_manual(values = genotype_colors) +
  scale_x_continuous(limits = x_range) +
  labs(x = "Pearson dispersion (phi)", y = NULL) +
  base_theme +
  theme(legend.position = "none")

p_ridge / p_box + plot_layout(heights = c(2, 1))

ggsave("Figures/Disperion_distriubtion.svg",
       height = 6, width = 4, units = "cm")

# null distributions ------------------------------------------------------
null_long <- as.data.frame(null_mat) %>%
  pivot_longer(everything(), names_to = "comparison", values_to = "stat")

obs_df <- data.frame(comparison = names(obs), obs_stat = as.numeric(obs))

ggplot(null_long, aes(x = stat)) +
  geom_histogram(bins = 60, fill = "grey75", color = "white") +
  geom_vline(data = obs_df, aes(xintercept = obs_stat), color = "firebrick", linewidth = 1) +
  facet_wrap(~ comparison, scales = "free", ncol = 1) +
  labs(x = "log(phi ratio)", y = "Count",
       title = "Exact permutation null distributions (observed value in red)") +
  theme_minimal()
