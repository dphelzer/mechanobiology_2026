library(spatstat)
library(terra)
library(tidyverse)

PX  <- 2  # match the Groovy export
DIR <- "spatial_export"


# build lookup table ------------------------------------------------------
genotype <- tribble(
  ~id,    ~Group,
  "2065", "mdx", "2066", "mdx", "2067", "mdx", "2068", "mdx", "2069", "mdx",
  "2585", "TG",  "2708", "TG",  "2717", "TG",  "2742", "TG",
  "2832", "TG",  "3148", "TG",  "3149", "TG"
)

ids <- str_remove(
  basename(list.files(DIR, "_mask\\.tif$")), 
  "_mask\\.tif$"
  )

ids <- ids[str_extract(ids, "^\\d+") %in% genotype$id]
ids


# load samples ------------------------------------------------------------
load_sample <- function(s) {
  m  <- as.matrix(rast(file.path(DIR, paste0(s, "_mask.tif"))), wide = TRUE)
  ny <- nrow(m); nx <- ncol(m)
  xr <- c(0, nx * PX); yr <- c(0, ny * PX)
  
  mac <- read_csv(file.path(DIR, paste0(s, "_macrophages.csv")), show_col_types = FALSE)
  
  list(id  = s,
       W   = owin(xrange = xr, yrange = yr, mask = m >= 1),
       IgG  = owin(xrange = xr, yrange = yr, mask = m == 2),
       Macs = mac)
}

samples <- set_names(map(ids, load_sample), ids)

# this is a large object, save so we don't have to run every time
saveRDS(samples, file = "samples.rds")

samples <- readRDS("samples.rds")

# QC and summary stats ----------------------------------------------------
check <- map_dfr(samples, function(s) {
  tibble(
    id        = s$id,
    nx        = ncol(s$W$m),
    ny        = nrow(s$W$m),
    width_um  = ncol(s$W$m) * PX,
    height_um = nrow(s$W$m) * PX,
    tissue_mm2 = sum(s$W$m) * PX^2 / 1e6,
    igg_pct   = 100 * sum(s$IgG$m) / sum(s$W$m),
    n_mac     = nrow(s$Macs),
    frac_in   = mean(inside.owin(s$Macs$x_um, s$Macs$y_um, s$W))
  )
})
print(check, n = Inf)

# testing on one sample ---------------------------------------------------
## distance to IgG regions
s <- samples[["2067-Stitching-13"]] 

system.time({
  D     <- distmap(s$IgG)
  dvals <- D[s$W, drop = TRUE]
})

length(dvals) # should equal sum(s$W$m)
sum(s$W$m)
summary(dvals) # min 0 (inside IgG), max = farthest tissue, in um

## for macrophages
keep <- inside.owin(s$Macs$x_um, s$Macs$y_um, s$W)
mac  <- ppp(s$Macs$x_um[keep], s$Macs$y_um[keep], window = s$W, check = FALSE)
dmac <- D[mac]

summary(dmac)
mean(dmac); mean(dvals)      # macrophage mean vs tissue mean

## percentile transform
dsort <- sort(dvals)
U <- findInterval(dmac, dsort) / length(dsort)

hist(U, breaks=10)
summary(U)
mean(U) # example calculation of test statistic
median(U)

# most sit very close to IgG region, 
# pulled to farther distances in the mean calculation
# keep both
# otherwise, everything seems ok

# scale to all samples ----------------------------------------------------
sample_U <- function(s) {
  D     <- distmap(s$IgG)
  dvals <- D[s$W, drop = TRUE]
  dsort <- sort(dvals)
  
  keep <- inside.owin(s$Macs$x_um, s$Macs$y_um, s$W)
  mac  <- ppp(
    s$Macs$x_um[keep], 
    s$Macs$y_um[keep],
    window = s$W, 
    check = FALSE
    )
  dmac <- D[mac]
  
  out <- tibble(
    id = s$id,
    x_um = mac$x,
    y_um = mac$y,
    dist_um = dmac,
    U = findInterval(dmac, dsort) / length(dsort),
    tissue_mean_dist = mean(dvals),
    n_outside = sum(!keep)
  )
  
  rm(D, dvals, dsort); gc(verbose = FALSE)
  out
}

U_all <- map_dfr(samples, sample_U)

# attach genotype
U_all <- U_all %>%
  mutate(sample_id = str_extract(id, "^\\d+")) %>%
  left_join(genotype, by = c("sample_id" = "id")) %>%
  mutate(Group = factor(Group, levels = c("mdx", "TG")))

count(U_all, id, Group) |> ggplot(aes(x = Group, y = n)) +
  geom_boxplot(outliers = FALSE) + geom_jitter()

openxlsx::write.xlsx(U_all, file = "distance_to_nearest_IgG.xlsx")

# plot distances ----------------------------------------------------------
## macrophages and distances
per_sample <- U_all |>
  group_by(id, Group) |>
  summarise(
    avg_U = mean(U),
    median_U = median(U), 
    n_mac = n(),
    .groups = "drop"
    )

## IgG area
geom <- map_dfr(samples, function(s) {
  tibble(
    id = s$id,
    tissue_mm2 = sum(s$W$m)   * PX^2 / 1e6,
    igg_mm2 = sum(s$IgG$m) * PX^2 / 1e6,
    igg_pct = 100 * sum(s$IgG$m) / sum(s$W$m)
  )
})

per_sample <- per_sample |>
  left_join(geom, by = "id") |>
  left_join(distinct(U_all, id, tissue_mean_dist), by = "id") |>
  mutate(mac_per_mm2 = n_mac / tissue_mm2)

genotype_colors <- c("mdx" = "red", "TG" = "blue")

base_theme <- theme_minimal(base_size = 8) +
  theme(
    text        = element_text(colour = "black", size = 8),
    axis.text   = element_text(colour = "black", size = 8),
    axis.title  = element_text(colour = "black", size = 8),
    axis.line   = element_line(colour = "black", linewidth = 0.35),
    axis.ticks  = element_line(colour = "black", linewidth = 0.35),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.ticks.length.x = unit(0, "cm"),
    panel.grid       = element_blank(),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    legend.position  = "none"
  )

## plot summary data ----
ggplot(per_sample, aes(Group, median_U, colour = Group)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.35) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.45,
               colour = "black", linewidth = 0.35) +
  geom_jitter(width = 0.12, height = 0, size = 1.5, shape = 16) +
  scale_colour_manual(values = genotype_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "Mean tissue percentile", x = NULL) +
  base_theme

# ggsave("Figures/tstar.svg", height = 4, width = 3.5, units = "cm")

## plot full distributions ----
B  <- 10
bw <- 1 / B

enr <- U_all |>
  mutate(bin = pmin(floor(U / bw) + 1, B)) |>
  count(id, Group, bin) |>
  complete(nesting(id, Group), bin = 1:B, fill = list(n = 0)) |>
  group_by(id) |>
  mutate(E = (n / sum(n)) / bw) |>
  ungroup() |>
  mutate(u_mid = (bin - 0.5) * bw)

enr_grp <- enr |>
  group_by(Group, u_mid) |>
  summarise(m = mean(E), s = sd(E), .groups = "drop")

ggplot(enr_grp, aes(u_mid, m, colour = Group, fill = Group)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35) +
  geom_ribbon(aes(ymin = m - s, ymax = m + s), colour = NA, alpha = 0.18) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = genotype_colors) +
  scale_fill_manual(values = genotype_colors) +
  labs(x = "Tissue percentile of distance to IgG+  (u)",
       y = "Macrophage enrichment (observed / expected)") +
  base_theme

# mac in IgG regions -------------------------------------------------------------
areas <- map_dfr(samples, \(s) tibble(
  id = s$id,
  tissue_um2 = sum(s$W$m) * PX^2,
  igg_um2 = sum(s$IgG$m) * PX^2
))

sweep_r <- function(s, rs = c(0, 10, 20, 30, 40, 50), PX = 2) {
  D <- distmap(s$IgG)
  dvals <- D[s$W, drop = TRUE]
  keep <- inside.owin(s$Macs$x_um, s$Macs$y_um, s$W)
  mac <- ppp(
    s$Macs$x_um[keep],
    s$Macs$y_um[keep], 
    window = s$W, 
    check = FALSE
    )
  dmac  <- D[mac]
  A_tot <- length(dvals) * PX^2
  N <- length(dmac)
  
  out <- map_dfr(rs, function(r) {
    A_in <- sum(dvals <= r) * PX^2
    N_in <- sum(dmac  <= r)
    tibble(id = s$id, r = r, N = N, N_in = N_in,
           area_in_mm2 = A_in / 1e6,
           lam_in = N_in / (A_in / 1e6),
           lam_out = (N - N_in) / ((A_tot - A_in) / 1e6),
           frac_in = N_in / N,
           frac_area = A_in / A_tot)
  })
  rm(D, dvals); gc(verbose = FALSE)
  mutate(out, RR = lam_in / lam_out)
}

rsweep <- map_dfr(samples, sweep_r) |>
  mutate(sample_id = str_extract(id, "^\\d+")) |>
  left_join(genotype, by = c("sample_id" = "id"))

print(rsweep, n = Inf, width = Inf)

openxlsx::write.xlsx(rsweep, file = "mac_near_IgG.xlsx")

rsweep |>
  group_by(Group, r) |>
  summarise(across(c(lam_in, lam_out, RR, frac_in, frac_area), mean), .groups = "drop") |>
  print(n = Inf)

ggplot(rsweep |> filter(r == 20), aes(x=Group, y=lam_in)) +
  geom_boxplot() + geom_jitter()

## stat test ----
wide <- rsweep |>
  filter(r > 0) |>
  select(sample_id, Group, r, RR) |>
  pivot_wider(names_from = r, values_from = RR) # ratio of in to out

Y <- as.matrix(wide[, -(1:2)]) # values only
g <- wide$Group # group labels
N <- nrow(Y); n1 <- sum(g == "mdx") # sample sizes

welch <- function(i) {
  a <- Y[i, , drop = FALSE]; b <- Y[-i, , drop = FALSE]
  abs(colMeans(a) - colMeans(b)) /
    sqrt(apply(a, 2, var) / nrow(a) + apply(b, 2, var) / nrow(b))
}

idx <- combn(N, n1, simplify = FALSE) # arrangements for exact test
obs <- welch(which(g == "mdx")) # obs welch-t at each r

# null at each r
null <- vapply(idx, welch, numeric(ncol(Y)))
# rows contain welch t for a given radius across all permutations
# columns have welch t across all radii for a given permutation

mx <- apply(null, 2, max) # get the max stat per permutation across
                          # all tested r 

# number of times max null stat (maximized over all r) is
# >= the max obs stat
cat(
  sprintf(
    "max|t| = %.3f global p = %.4f \n",
    max(obs), 
    mean(mx >= max(obs))))

# maxT controlled p-val on welch t stat 0.0139

# p val at each, few controlled by maxT
openxlsx::write.xlsx(
  x = tibble(
    r = as.numeric(colnames(Y)), 
    t = obs,
    p_adj = sapply(obs, \(o) mean(mx >= o))
  ),
  file = "maxT_all_radii.xlsx"
)

# alternative to avoid one r giving the max stat
# integrated test stat then global pvalue
S <- function(x) mean(x^2)
obsS <- S(obs)
nullS <- apply(null, 2, S)
cat(
  sprintf(
    "T = %.3f   global p = %.4f", 
    obsS, 
    mean(nullS >= obsS))
  )

# p = 0.0139
# same result, more robust when there are differences across radii

# lesion size & # of patches -------------------------------------------------------------
lesion_stats <- function(s, PX = 2) {
  r <- rast(file.path(DIR, paste0(s$id, "_mask.tif")))
  ig <- ifel(r == 2, 1, NA)
  p <- patches(ig, directions = 8, zeroAsNA = FALSE)
  f <- freq(p)
  a <- f$count * PX^2
  tibble(id = s$id, n_lesions = length(a),
         mean_lesion_um2   = mean(a),
         median_lesion_um2 = median(a),
         max_lesion_um2    = max(a),
         total_igg_um2     = sum(a))
}

lesions <- map_dfr(samples, lesion_stats) |>
  mutate(sample_id = str_extract(id, "^\\d+")) |>
  left_join(genotype, by = c("sample_id" = "id"))

print(lesions, n = Inf, width = Inf)

openxlsx::write.xlsx(lesions, "IgG_summary.xlsx")

lesions |> ggplot(aes(x=Group, y=mean_lesion_um2)) +
  geom_boxplot(outliers = F, width = 0.3) +
  geom_jitter(width = 0.1) +
  base_theme


max_res <- perm_maxT(
  lesions, 
  formula = max_lesion_um2~Group,
  stat = "welch",
  resampling = "exact"
  )

mean_res <- perm_maxT(
  lesions, 
  formula = mean_lesion_um2~Group,
  stat = "welch",
  resampling = "exact"
)

openxlsx::write.xlsx(
  x = list(
    max_lesion_results = max_res,
    mean_lesion_results = mean_res
  ),
  file = "welch_maxT_lesion_stats.xlsx"
)