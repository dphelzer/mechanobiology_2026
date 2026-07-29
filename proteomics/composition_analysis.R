install.packages("compositions")
install.packages("ggfortify")
install.packages("RVAideMemoire")
library(compositions)
library(tidyverse)
library(ggfortify)
library(vegan)
library(RVAideMemoire)
library(patchwork)
library(purrr)
library(broom)

# prepare data ------------------------------------------------------------

data <- openxlsx::read.xlsx(
  "QCONCATS.SKELETAL MUSCLE.ALL GENOTYPES.Skeletal Muscle Absolute Quant. Summary.xlsx",
  sheet = 6)

data <- data |> 
  pivot_longer(
    cols = WT1:TG5,
    values_to = "Fraction",
    names_to = "Samples"
    )

data_list <- data |> split(f = data$Annotation)

data_list <- lapply(
  data_list,
  function(x){
    x |> 
      pivot_wider(
        names_from = Gene,
        values_from = Fraction
      ) |> 
      select(-Annotation) |> 
      column_to_rownames(var = "Samples")
  }
)

clr_data <- lapply(
  data_list,
  function(x){
    clr(x)
  }
)


## pca and omnibus permanova -----------------------------------------------
pca_list <- lapply(
  clr_data,
  function(x){
    prcomp(x, scale. = FALSE)
  }
)

for (i in seq(pca_list)) {
  message("PCA results for ", names(pca_list[i]))
  print(summary(pca_list[[i]]))
}

# sample labels
sample <- unique(data$Samples)

# plot pca for all
plot_list <- list()
for (i in seq(pca_list)) {
  plot_list[[i]] <- autoplot(
    pca_list[[i]],
    data = data.frame(Group = sample),
    colour = "Group",
    loadings = FALSE,
    loadings.label = FALSE,
    shape = 19,
    size = 3
  ) +
    labs(title = names(pca_list)[i]) +
    scale_color_manual(
      values = c(rep("grey10",5), rep("red",4), rep("blue",5))
      ) +
    theme_bw(base_size = 8,
             base_line_size = 0.35) +
      theme(legend.position = "none")
}


wrap_plots(plot_list, ncol = 5)

# stat test global composition
PERMANOVA <- function(dat, samples){
  adonis2(
    as.matrix(dat) ~ samples,
    method = "euclidean",
    permutations = 99999
  )
}

groups <- c(rep("WT", 5), rep("mdx", 4), rep("TG", 5))

set.seed(729)
res_list <- lapply(
  clr_data, PERMANOVA, samples = groups
)

for (i in seq(res_list)) {
  message("PERMANOVA Results: ", names(res_list[i]))
  print(res_list[[i]])
}

## save all omnibus results
res_df <- lapply(res_list, as.data.frame)
openxlsx::write.xlsx(
  x = res_df,
  file = "permanova_results.xlsx"
  )

## pairwise contrasts for remainder ---------------------------------------
## remove non-significant

res_filtered <- keep(
  res_df, 
  ~ any(
    .x$`Pr(>F)` < 0.05, 
        na.rm = TRUE
    )
  )


clr_data_filt <- keep(
  clr_data, 
  names(clr_data) %in% names(res_filtered)
  )

pairwise_permanova_test <- function(dat, samples) {
  pairwise.perm.manova(
    dist(as.matrix(dat),
         method = "euclidean"),
    samples,
    nperm = 99999,
    p.method = "holm"
  )
}

set.seed(757)
pairwise_results <- lapply(
  clr_data_filt,
  pairwise_permanova_test,
  samples = groups
)

for (i in seq(pairwise_results)) {
  message("Pairwise results: ", names(pairwise_results[i]))
  print(pairwise_results[[i]])
}

## save all pairwise results
res_pairwise <- lapply(pairwise_results, tidy)
openxlsx::write.xlsx(
  x = res_pairwise,
  file = "pairwise_permanova_results.xlsx"
)