BiocManager::install("DESeq2")
library(Seurat)
library(patchwork)
library(tidyverse)
library(clustree)

# pre-processing ####
## dba wt and mdx
## load in data

counts <- list() # create list for countts
numsamples <- 6 # number of samples

## names
sam_names <- c("WT7", "WT28", 
               "D2.mdx7", "D2.mdx28", 
               "D2.mdx7_steroid", "D2.mdx28_steroid")

## create all
for (i in 1:numsamples) {
  counts[[i]] <- ReadMtx(mtx = paste(sam_names[i], "/matrix.mtx.gz", sep = ""),
                         features = paste(sam_names[i], "/features.tsv.gz", sep = ""),
                         cells = paste(sam_names[i], "/barcodes.tsv.gz", sep = ""))
}

## seurat object
object.list <- list()
for (i in 1:length(counts)) {
  object.list[[i]] <- CreateSeuratObject(counts = counts[[i]],
                                         min.cells = 3,
                                         min.features = 200,
                                         project = sam_names[i])
}

names(object.list) <- sam_names

## mt percent
for (i in 1:length(object.list)) {
  object.list[[i]][["percent.mt"]] <- PercentageFeatureSet(object.list[[i]], pattern = "^mt")
}
head(object.list[[1]]@meta.data,3) # good

plot.list <- list()
for (i in 1:length(object.list)) {
  plot.list[[i]] <- VlnPlot(object.list[[i]],
                            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                            pt.size = 0) &
    geom_boxplot(width = 0.1, fill = "white")
}

wrap_plots(plot.list)

## don't need to filter for percent mito, these are nuclei

thresholds <- function(seurat_object) {
  low_RNAcount <- quantile(seurat_object[['nCount_RNA']]$nCount_RNA, probs = 0.02)
  high_RNAcount <- quantile(seurat_object[['nCount_RNA']]$nCount_RNA, probs = 0.98)
  output <- c(low_RNAcount, high_RNAcount)
  return(output)
}

QC_thresholds <- sapply(object.list, thresholds)
QC_thresholds <- as.data.frame(QC_thresholds)
rownames(QC_thresholds) <- c("Low UMI", "High UMI")
View(QC_thresholds)

for (i in 1:length(object.list)) {
  plot.list[[i]] <- VlnPlot(object.list[[i]], 
                            features = c("nFeature_RNA", "nCount_RNA"), 
                            ncol = 2, pt.size = 0) &
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.length.y = unit(-0.1, "cm")) &
    geom_boxplot(width = 0.1, fill = "white") &
    geom_hline(yintercept = QC_thresholds[[i]][1]) &
    geom_hline(yintercept = QC_thresholds[[i]][2])
}

wrap_plots(plot.list)

for (i in 1:length(object.list)) {
  object.list[[i]] <- subset(
    x = object.list[[i]],
    subset = nCount_RNA > QC_thresholds[[i]][1] &
      nCount_RNA < QC_thresholds[[i]][2]
  )
}

for (i in 1:length(object.list)) {
  plot.list[[i]] <- FeatureScatter(object.list[[i]], 
                                   feature1 = "nCount_RNA", 
                                   feature2 = "nFeature_RNA") & NoLegend() & 
    theme(axis.ticks.length = unit(-0.1, "cm"))
}

gridExtra::grid.arrange(grobs = plot.list, ncol = 3)

# merge data ####
merged <- merge(x = object.list[[1]], 
                y = c(object.list[[2]], object.list[[3]], 
                      object.list[[4]], object.list[[5]],
                      object.list[[6]]))

## processing steps ####
merged <- SCTransform(merged)
merged <- RunPCA(merged)
merged <- RunUMAP(merged, dims = 1:30, verbose = TRUE)

head(merged@meta.data,3)
unique(merged$orig.ident)
merged@active.ident <- factor(merged$orig.ident, levels = unique(merged$orig.ident))
DimPlot(merged) # cool
DimPlot(merged, label = TRUE) + NoLegend()

## clustering
merged <- FindNeighbors(merged, dims = 1:30)

res <- seq(0.1,1,0.1)
for (i in 1:length(res)) {
  merged <- FindClusters(merged, algorithm = 4, resolution = res[i], random.seed = 135)
}

head(merged@meta.data,3)

clustree(x = merged, prefix = "SCT_snn_res.", node_colour = "sc3_stability") # use 0.2

# 0.2 is fine, remove the rest
for (res in c(0.1, seq(0.3, 1, 0.1))) {
  column <- paste0("SCT_snn_res.", res)
  merged[[column]] <- NULL
}

merged$seurat_clusters <- merged$SCT_snn_res.0.2
merged@active.ident <- merged$seurat_clusters

## annotation
DefaultAssay(merged) <- "RNA"
merged <- NormalizeData(merged)
merged <- JoinLayers(merged)

saveRDS(object = merged, file = "merged.rds")

DE <- FindAllMarkers(merged, assay = "RNA", logfc.threshold = 0.25, only.pos = TRUE)
DE <- DE %>% filter(p_val_adj < 0.05)
DE <- DE %>% relocate(gene)

IDs <- c("FAPs", 2, "FAPs", "Myonuclei", "Myonuclei",
         6, 7, "SCs", 9, 10,
         11, "Tenocytes", 13, 14, 15)

# IDs <- c("FAPs", "FAPs", "Innate immune cells", "Myonuclei", 
#          "Myonuclei", "Innate immune cells", "ECs", "Innate immune cells", 
#          "Innate immune cells", "B cells", "T cells", "SCs", 
#          "Tenocytes", "Smooth muscle", 15, "Adipocytes")
names(IDs) <- levels(merged)
merged <- RenameIdents(merged, IDs)

saveRDS(merged, "merged_annotated.rds")
merged <- readRDS("merged_annotated.rds")

Genotype <- c(
  rep("WT", length(which(merged$orig.ident == "WT7")) + 
        length(which(merged$orig.ident == "WT28"))),
  
  rep("mdx", length(which(merged$orig.ident == "D2.mdx7")) + 
        length(which(merged$orig.ident == "D2.mdx28"))),
  
  rep("mdx_steroid", length(which(merged$orig.ident == "D2.mdx7_steroid")) + 
        length(which(merged$orig.ident == "D2.mdx28_steroid")))
  )

merged$Genotype <- factor(Genotype, levels = c("WT", "mdx", "mdx_steroid"))

## just plotting, skipping annotations
genes <- c("Fkrp", "Dolk", "B3galnt2", "Fktn", "B4gat1", 
           "Pomgnt1", "Pomgnt2", "Pomk", "Pomt1", 
           "Crppa", "Pomt2", "Gmppb",
           "Large1", "Dpm1", "Dpm2", "Dpm3")

merged$orig.ident <- factor(merged$orig.ident, levels = c(unique(merged$orig.ident)))
plot.list <- list()
for (i in 1:length(genes)) {
  plot.list[[i]] <- VlnPlot(merged, features = genes[i], pt.size = 0, split.by = "Genotype"
  ) & 
    #geom_boxplot(width = 0.1, fill = "white", outliers = FALSE) &
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 15, hjust = 0.5, angle = 0),
          axis.ticks.x = element_blank(),
          axis.ticks.length.y = unit(-0.1, "cm"),
          axis.text.y = element_text(size = 15),
          axis.title = element_text(size = 15)) & 
    #NoLegend() &
    ylab("")
}

vlnplots <- gridExtra::grid.arrange(grobs = plot.list, ncol = 4)
ggsave(plot = vlnplots, filename = "Figures/ridge.tiff", height = 15, 
       width = 25, units = "cm", dpi = 600)

# FAPs --------------------------------------------------------------------
FAPs <- subset(merged, idents = "FAPs")
FAPs@active.ident <- FAPs$Genotype

## same plot as human
GROUP_COL    <- "Genotype"
GROUP_LEVELS <- c("WT", "mdx", "mdx_steroid")
COMPARISONS  <- list(c("WT", "mdx"), c("WT", "mdx_steroid"), c("mdx", "mdx_steroid"))
COLOR_LIMIT  <- 0.15

gene_sets <- list(
  "Core" = c("YAP1","WWTR1","TEAD1","TEAD2","TEAD3","TEAD4"),
  "Negative regulators" = c("STK3","STK4","SAV1","MOB1A","MOB1B","NF2","AMOT",
                            "LATS1","AMOTL1","WWC1","WWC2","PTPN14","VGLL4",
                            "FRMD6","MAP4K4"),
  "YAP/TAZ/TEAD targets" = c("ANKRD1","CCN1","CCN2","CDC42EP3","PARVA","NEDD9",
                             "CSRP1","FMN2","THBS1","DCN","SPARC","COL6A2","MFAP5",
                             "LTBP1","PXDN","INHBA","TGFB2","CDH13","DAB2","MYOF",
                             "ASAP1","CCDC80","DOCK5","ARHGEF17","NUAK2","LATS2",
                             "AMOTL2","AXL","CRIM1","F3","FJX1","FOXF2","GADD45B",
                             "IGFBP3","NT5E","TNS1"),
  "Mechanotransduction" = c("PIEZO1","PIEZO2","TRPV4",
                            "ITGB1","ITGA5","ITGAV","ITGA11","PTK2","SRC","TLN1","VCL","ZYX","PXN",
                            "RHOA","ROCK1","ROCK2","MYH9","MYL9","ACTB","ACTA2","FLNA",
                            "ARHGAP29","ARHGAP35","GNA12","GNA13","LPAR1",
                            "MRTFA","MRTFB","SRF",
                            "CAV1","LMNA","SYNE1","SYNE2","EMD")
)

to_mouse <- function(x) paste0(substr(x, 1, 1), tolower(substring(x, 2)))
gene_sets <- lapply(gene_sets, to_mouse)

# U effect = U/(n1*n2) - 0.5; positive = higher in second group
u_effect_table <- function(expr, groups, comparisons, genes) {
  genes <- intersect(genes, rownames(expr))
  bind_rows(lapply(comparisons, function(cmp) {
    a <- which(groups == cmp[1]); b <- which(groups == cmp[2])
    n1 <- length(a); n2 <- length(b); N <- n1 + n2
    sub <- expr[genes, c(a, b), drop = FALSE]
    df <- bind_rows(lapply(genes, function(g) {
      x <- as.numeric(sub[g, ])
      U_b <- sum(rank(x)[(n1 + 1):N]) - n2 * (n2 + 1) / 2
      tie <- table(x); tie <- sum(tie^3 - tie)
      sigma <- sqrt(n1 * n2 / 12 * ((N + 1) - tie / (N * (N - 1))))
      z <- if (sigma > 0) (U_b - n1 * n2 / 2) / sigma else 0
      data.frame(gene = g, U_effect = U_b / (n1 * n2) - 0.5, p = 2 * pnorm(-abs(z)))
    }))
    df$comparison <- paste(cmp[1], "vs", cmp[2])
    df$p_adj <- p.adjust(df$p, "BH")
    df
  }))
}

plot_u_heatmap <- function(tab, gene_sets, limit = 0.2) {
  txt <- element_text(size = 8, colour = "black")
  plots <- lapply(names(gene_sets), function(s) {
    d <- tab %>%
      filter(gene %in% gene_sets[[s]]) %>%
      mutate(gene = factor(gene, levels = rev(gene_sets[[s]])),
             comparison = factor(comparison, levels = unique(tab$comparison)),
             fill = pmax(pmin(U_effect, limit), -limit),
             star = ifelse(p_adj < 0.05, "*", ""),
             set = s)
    ggplot(d, aes(comparison, gene, fill = fill)) +
      geom_tile(colour = "white", linewidth = 0.4) +
      geom_text(aes(label = star), size = 8 / .pt, colour = "black", vjust = 0.75) +
      scale_fill_distiller(palette = "RdBu", limits = c(-limit, limit),
                           breaks = seq(-limit, limit, length.out = 5),
                           name = expression(italic(U)*"-statistic")) +
      facet_grid(set ~ ., switch = "y") +
      coord_fixed(expand = FALSE) +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_size = 8) +
      theme(text = txt, axis.text = txt,
            axis.text.y = element_text(size = 8, colour = "black", face = "italic"),
            axis.text.x = if (s == tail(names(gene_sets), 1))
              element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5)
            else element_blank(),
            legend.text = txt, legend.title = txt,
            strip.text.y.left = element_text(size = 8, colour = "black", angle = 90, face = "bold"),
            strip.placement = "outside",
            panel.grid = element_blank())
  })
  wrap_plots(plots, ncol = 1, heights = lengths(gene_sets), guides = "collect") &
    theme(legend.position = "top", legend.key.width = unit(0.4, "cm"))
}

expr   <- GetAssayData(FAPs, assay = "RNA", layer = "data")
groups <- as.character(FAPs@meta.data[[GROUP_COL]])
keep   <- groups %in% GROUP_LEVELS

tab <- u_effect_table(expr[, keep], groups[keep], COMPARISONS, unlist(gene_sets))
p   <- plot_u_heatmap(tab, lapply(gene_sets, intersect, rownames(expr)), COLOR_LIMIT)
p

ggsave(filename = "Figures/YAP_TAZ_heatmap.svg",
       height = 35, width = 6, units = "cm")