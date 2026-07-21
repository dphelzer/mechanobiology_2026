install.packages("Rfast2")
if (!requireNamespace("remotes", quietly = T)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
BiocManager::install("SPOTlight")
BiocManager::install("BayesSpace")

library(Seurat)
library(ggplot2)
library(patchwork)
library(SPOTlight)
library(tidyverse)
library(clustree)
library(SeuratData)
library(SeuratDisk)
library(enrichR)
library(SingleCellExperiment)
library(BayesSpace)

# Guide from https://satijalab.org/seurat/articles/spatial_vignette.html

### WT data ####
#### Load in the data ####
# need to have the filtered_feature_bc_matrix.h5 file and "spatial" folder in the directory
WT = Load10X_Spatial("WT")
WT=readRDS("WT_processed.rds")

#### Data Processing ####
# plot data
p1 = VlnPlot(
  WT, 
  features = "nCount_Spatial", 
  pt.size = 0.1
  ) + 
  NoLegend()
p2 = SpatialFeaturePlot(
  WT, 
  features = "nCount_Spatial"
  ) + 
  theme(legend.position = "right")

wrap_plots(p1, p2)

# Normalize
WT = SCTransform(WT, assay = "Spatial", verbose = T)

#### Gene Expression Visualization ####
# plot any features of interest
# use pt.size.factor = to change the size of the spots
## default is 1.6
# use alpha = c(1,1) to change the transparency
## default is c(1,1)
## changing this to, for example, c(0.1,1) lowers transparency of lower expressed genes
SpatialFeaturePlot(WT, features = c("Col1a1","Lgals3","Cd68"), alpha = c(1,1), pt.size.factor = 2.8,
                   stroke = 2)
SpatialFeaturePlot(
  WT, 
  features = "Myh4", 
  pt.size.factor = 2.8,
  stroke = 2)

SpatialFeaturePlot(
  WT, 
  features = c("Yap1"), 
  pt.size.factor = 2.8,
  stroke = 2,
  alpha = c(1,1)
)


#### Dimensionality reduction, clustering, and visualization ####
WT = RunPCA(WT, assay = "SCT", verbose = T)
WT = FindNeighbors(WT, reduction = "pca", dims = 1:30)

# find optimal clusters by repeating with resolution form 0 to 1 in 0.1 increments
WT = FindClusters(WT,
                  algorithm = 4,
                  verbose = T,
                  resolution = 0.6)

# use clustree to see clusters
head(WT@meta.data,3)
clustree(WT, prefix = "SCT_snn_res.", node_colour = "sc3_stability")

# remove unwanted clusters
WT[['SCT_snn_res.1']] = NULL

WT = RunUMAP(WT, reduction = "pca", dims = 1:30)

# save the file
SaveH5Seurat(WT, "WT_processed.h5Seurat") 
WT = LoadH5Seurat("WT_processed.h5Seurat") # cannot load back in image data
saveRDS(WT, "WT_processed.rds")

# visualize with DimPlot and by overlaying with image
p1 = DimPlot(WT, 
             reduction = "umap", 
             label = T,
             repel = F, 
             pt.size = 1.5, 
             label.box = T
             )+
  NoLegend()
p2 = SpatialDimPlot(WT, 
                    label = T, 
                    label.size = 3,
                    pt.size.factor = 2.8,
                    stroke = 2)
p1+p2

# facet by cluster for easier viewing
SpatialDimPlot(
  WT, 
  cells.highlight = CellsByIdentities(object = WT, idents = c(2,4,5,6,7)), 
  facet.highlight = T, 
  ncol = 3
  )

#### Interactive Plotting ####
# both SpatialDimPlot() and SpatialFeaturePlot() have an interactive mode
# set interactive = T to open a pane in RStudio
## SpatialDimPlot interaction can hover over spots and view the cell name and current identity class
## SpatialFeaturePlot interaction can adjust the transparency of the spots, the point size, as well as the Assay
# these are a little slow, but can be cool
SpatialDimPlot(WT, interactive = T)
SpatialFeaturePlot(WT, features = "Myh4", interactive = T)

# LinkedDimPlot() links the UMAP representation to the tissue image representation and allows interaction
# For example, you can select a region in the UMAP plot and the corresponding spots in the image representation will be highlighted
# again, it can be slow, but might be cool to explore
LinkedDimPlot(WT)

#### Differential Expression ####
##### standard method by finding marker genes ####
all_markers = FindAllMarkers(WT)
all_markers = all_markers %>% relocate(gene)
write.csv(all_markers, "WT_clusters.csv")

##### Finding genes based on spatial patterning ####
# search for features exhibiting spatial patterning in the absence of pre-annotation
?FindSpatiallyVariableFeatures
WT = FindSpatiallyVariableFeatures(
  WT, 
  assay = "SCT", 
  features = VariableFeatures(WT)[1:1000],
  selection.method = "markvariogram",
  verbose = T
  )

# Visualize the expression of the top features identified by this measure
# For sample, 10 features
top_features_WT = head(SpatiallyVariableFeatures(WT, selection.method = "markvariogram"), 10)
SpatialFeaturePlot(
  WT, 
  features = top_features_WT, 
  ncol = 5,
  alpha = c(1, 1),
  pt.size.factor = 2.8
  )

#### Integration with single cell data ####
browseVignettes("SPOTlight")
mdx <- readRDS("mdx_processed.rds")
mdx_sc <- LoadH5Seurat("/Users/pranav/Desktop/mdx nuclei /processed_mdx_annotated.h5Seurat")
mdx_marker_genes <- readRDS("marker genes for deconvolution.rds")
mdx_sc_sample <- readRDS("mdx nuclei sample.rds")

hvg <- mdx_sc_sample@assays$SCT@var.features
mdx_marker_genes <- FindAllMarkers(
  mdx_sc,
  assay = "SCT",
  slot = "data",
  verbose = T,
  only.pos = T
)

mdx_marker_genes <- mdx_marker_genes %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.5) %>% 
  arrange(desc(avg_log2FC))
saveRDS(mdx_marker_genes, "marker genes for deconvolution.rds")

cluster_ID_mdx <- levels(mdx_sc)
cluster_ID_mdx = c(
  "Myonuclei", "Myonuclei", "FAPs", "Myocytes", 
  "Leukocytes", "SCs", "MTJ", "Myoblasts",
  "ECs", "Smooth Muscle", "FAPs", "NMJ",
  "Adipocytes"
  )
names(cluster_ID_mdx) <- levels(mdx_sc)
mdx_sc <- RenameIdents(mdx_sc, cluster_ID_mdx) 
saveRDS(mdx_sc, "mdx nuclei for deconvolution.rds")

idx <- split(seq(ncol(mdx_sc)), Idents(mdx_sc))
n_cells <- 100
cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells)
    n_cells <- n
  sample(i, n_cells)
})
mdx_sc_sample <- mdx_sc[, unlist(cs_keep)]
saveRDS(mdx_sc_sample, "mdx nuclei sample.rds")

mdx_decon <- SPOTlight(
  x = mdx_sc_sample,
  y = mdx,
  groups = Idents(mdx_sc_sample),
  mgs = mdx_marker_genes,
  hvg = hvg,
  weight_id = "avg_log2FC",
  group_id = "cluster",
  gene_id = "gene"
  )

saveRDS(mdx_decon, "mdx decon.rds")

### mdx data ####
#### Load in the data ####
# need to have the filtered_feature_bc_matrix.h5 file and "spatial" folder in the directory
mdx = Load10X_Spatial("mdx")
mdx=readRDS("mdx_processed.rds")

#### Data Processing ####
# plot data
p1 = VlnPlot(
  mdx, 
  features = "nCount_Spatial", 
  pt.size = 0.1
) + 
  NoLegend()
p2 = SpatialFeaturePlot(
  mdx, 
  features = "nCount_Spatial"
) + 
  theme(legend.position = "right")

wrap_plots(p1, p2)

# Normalize
mdx = SCTransform(mdx, assay = "Spatial", verbose = T) # error message suggests some cells have 0 UMI counts
# SCTransform assumes all cells have >=1 UMI
# filter out any that have 0 UMI counts with code below, then re-run SCTransform
min(mdx$nCount_Spatial) # checks if there are 0 UMI counts in any cells
sum(mdx$nCount_Spatial == 0)
mdx = subset(mdx, nCount_Spatial>0) # filter out those with 0 UMI

#### Gene Expression Visualization ####
SpatialFeaturePlot(
  mdx, 
  features = c("Cd68","Mrc1","Lgals3"), 
  pt.size.factor = 2.8,
  stroke = 2)

SpatialFeaturePlot(
  mdx, 
  features = c("Lgals3", "Pdgfra", "Pax7", "Cd80"), 
  pt.size.factor = 2.8,
  stroke = 2,
  alpha = c(1,1)
)


#### Dimensionality reduction, clustering, and visualization ####
mdx = RunPCA(mdx, assay = "SCT", verbose = T)
mdx = FindNeighbors(mdx, reduction = "pca", dims = 1:30)

# find optimal clusters by repeating with resolution form 0 to 1 in 0.1 increments
mdx = FindClusters(mdx,
                  algorithm = 4,
                  verbose = T,
                  resolution = 0.6)

# use clustree to see clusters
head(mdx@meta.data,3)
clustree(mdx, prefix = "SCT_snn_res.", node_colour = "sc3_stability")

# remove unwanted clusters
mdx[['SCT_snn_res.1']] = NULL

mdx = RunUMAP(mdx, reduction = "pca", dims = 1:30)

# save the file
SaveH5Seurat(mdx, "mdx_processed.h5Seurat")
mdx = LoadH5Seurat("mdx_processed.h5Seurat")
saveRDS(mdx, "mdx_processed.rds")

# visualize with DimPlot and by overlaying with image
p1 = DimPlot(mdx, 
             reduction = "umap", 
             label = T,
             repel = F, 
             pt.size = 1.5, 
             label.box = T
)+
  NoLegend()
p2 = SpatialDimPlot(mdx, 
                    label = T, 
                    label.size = 3,
                    pt.size.factor = 2.8,
                    stroke = 2)
p1+p2

# facet by cluster for easier viewing
SpatialDimPlot(
  mdx, 
  cells.highlight = CellsByIdentities(object = mdx, idents = c(2,4,5,6,7)), 
  facet.highlight = T, 
  ncol = 3
)

#### Interactive Plotting ####
# both SpatialDimPlot() and SpatialFeaturePlot() have an interactive mode
# set interactive = T to open a pane in RStudio
## SpatialDimPlot interaction can hover over spots and view the cell name and current identity class
## SpatialFeaturePlot interaction can adjust the transparency of the spots, the point size, as well as the Assay
# these are a little slow, but can be cool
SpatialDimPlot(mdx, interactive = T)
SpatialFeaturePlot(mdx, features = "Myh4", interactive = T)

# LinkedDimPlot() links the UMAP representation to the tissue image representation and allows interaction
# For example, you can select a region in the UMAP plot and the corresponding spots in the image representation will be highlighted
# again, it can be slow, but might be cool to explore
LinkedDimPlot(mdx)

#### Differential Expression ####
##### standard method by finding marker genes ####
all_markers = FindAllMarkers(mdx)
all_markers = all_markers %>% relocate(gene)
write.csv(all_markers, "mdx_clusters.csv")

##### Finding genes based on spatial patterning ####
# search for features exhibiting spatial patterning in the absence of pre-annotation
?FindSpatiallyVariableFeatures
mdx = FindSpatiallyVariableFeatures(
  mdx, 
  assay = "SCT", 
  features = VariableFeatures(mdx)[1:1000],
  selection.method = "markvariogram"
)

# Visualize the expression of the top features identified by this measure
# For sample, 10 features
top_features_mdx = head(SpatiallyVariableFeatures(mdx, selection.method = "markvariogram"), 10)
SpatialFeaturePlot(
  mdx, 
  features = top_features_mdx, 
  ncol = 5,
  alpha = c(1, 1),
  pt.size.factor = 2.8
)

#### Integration with single cell data ####
# since the single cell data is already processed, just load in the dataset here
# note the single cell data must be normalized and contain dim reductions pca and umap
# integrate the data
# reference is the single cell data
# query is the spatial
anchors = FindTransferAnchors(
  reference = mdx_SingeCell, 
  query = mdx, 
  normalization.method = "SCT"
)

# set refdata equal to the cluster annotations
predictions_assay = TransferData(
  anchorset = anchors, 
  refdata = mdx_SingleCell$seurat_clusters, 
  prediction.assay = T,
  weight.reduction = mdx[["pca"]], 
  dims = 1:30
)

mdx[["predictions"]] = predictions.assay

### mdxTG Data ####
#### Load in the data ####
# need to have the filtered_feature_bc_matrix.h5 file and "spatial" folder in the directory
mdx_TG = Load10X_Spatial("mdxTG")
mdxTG=readRDS("mdx_TG_processed.rds")

#### Data Processing ####
# plot data
p1 = VlnPlot(
  mdx_TG, 
  features = "nCount_Spatial", 
  pt.size = 0.1
) + 
  NoLegend()
p2 = SpatialFeaturePlot(
  mdx_TG, 
  features = "nCount_Spatial"
) + 
  theme(legend.position = "right")

wrap_plots(p1, p2)

# Normalize
mdx_TG = SCTransform(mdx_TG, assay = "Spatial", verbose = T) # error message suggests some cells have 0 UMI counts
# SCTransform assumes all cells have >=1 UMI
# filter out any that have 0 UMI counts with code below, then re-run SCTransform
min(mdx_TG$nCount_Spatial)
sum(mdx_TG$nCount_Spatial == 0) # some have 0 UMI
mdx_TG = subset(mdx_TG, nCount_Spatial>0)

#### Gene Expression Visualization ####
SpatialFeaturePlot(
  mdxTG, 
  features = c("Yap1","Ccn2","Pdgfra"), 
  pt.size.factor = 2.8,
  stroke = 2,
  alpha = c(0.5,1)
  )

SpatialFeaturePlot(
  mdxTG, 
  features = c("Sspn"), 
  pt.size.factor = 2.8,
  stroke = 2,
  alpha = c(1,1)
)


#### Dimensionality reduction, clustering, and visualization ####
mdx_TG = RunPCA(mdx_TG, assay = "SCT", verbose = T)
mdx_TG = FindNeighbors(mdx_TG, reduction = "pca", dims = 1:30)

# find optimal clusters by repeating with resolution form 0 to 1 in 0.1 increments
mdx_TG = FindClusters(mdx_TG,
                   algorithm = 4,
                   verbose = T,
                   resolution = 0.6)

# use clustree to see clusters
head(mdxTG@meta.data,3)
clustree(mdxTG, prefix = "SCT_snn_res.", node_colour = "sc3_stability")

# remove unwanted clusters
mdxTG[['SCT_snn_res.1']] = NULL

mdx_TG = RunUMAP(mdx_TG, reduction = "pca", dims = 1:30)

# save the file
SaveH5Seurat(mdx_TG, "mdx_TG_processed.h5Seurat")
saveRDS(mdx_TG,"mdx_TG_processed.rds")

# visualize with DimPlot and by overlaying with image
p1 = DimPlot(mdx_TG, 
             reduction = "umap", 
             label = T,
             repel = F, 
             pt.size = 1.5, 
             label.box = T
)+
  NoLegend()
p2 = SpatialDimPlot(mdx_TG, 
                    label = T, 
                    label.size = 3,
                    pt.size.factor = 2.8,
                    stroke = 2)
p1+p2

# facet by cluster for easier viewing
SpatialDimPlot(
  mdxTG, 
  cells.highlight = CellsByIdentities(object = mdxTG, idents = c(2,4,5,6,7)), 
  facet.highlight = T, 
  ncol = 3
)

#### Interactive Plotting ####
# both SpatialDimPlot() and SpatialFeaturePlot() have an interactive mode
# set interactive = T to open a pane in RStudio
## SpatialDimPlot interaction can hover over spots and view the cell name and current identity class
## SpatialFeaturePlot interaction can adjust the transparency of the spots, the point size, as well as the Assay
# these are a little slow, but can be cool
SpatialDimPlot(mdxTG, interactive = T)
SpatialFeaturePlot(mdxTG, features = "Myh4", interactive = T)

# LinkedDimPlot() links the UMAP representation to the tissue image representation and allows interaction
# For example, you can select a region in the UMAP plot and the corresponding spots in the image representation will be highlighted
# again, it can be slow, but might be cool to explore
LinkedDimPlot(mdxTG)

#### Differential Expression ####
##### standard method by finding marker genes ####
all_markers = FindAllMarkers(mdx_TG)
all_markers = all_markers %>% relocate(gene)
write.csv(all_markers, "mdx_TG_clusters.csv")

##### Finding genes based on spatial patterning ####
# search for features exhibiting spatial patterning in the absence of pre-annotation
?FindSpatiallyVariableFeatures
mdxTG = FindSpatiallyVariableFeatures(
  mdxTG, 
  assay = "SCT", 
  features = VariableFeatures(mdxTG)[1:1000],
  selection.method = "markvariogram"
)

# Visualize the expression of the top features identified by this measure
# For sample, 10 features
top_features_mdxTG = head(SpatiallyVariableFeatures(mdxTG, selection.method = "markvariogram"), 10)
SpatialFeaturePlot(
  mdxTG, 
  features = top_features_mdxTG, 
  ncol = 5,
  alpha = c(1, 1),
  pt.size.factor = 2.8
)

#### Integration with single cell data ####
# since the single cell data is already processed, just load in the dataset here
# note the single cell data must be normalized and contain dim reductions pca and umap
# integrate the data
# reference is the single cell data
# query is the spatial
anchors = FindTransferAnchors(
  reference = mdxTG_SingeCell, 
  query = mdxTG, 
  normalization.method = "SCT"
)

# set refdata equal to the cluster annotations
predictions_assay = TransferData(
  anchorset = anchors, 
  refdata = mdxTG_SingleCell$seurat_clusters, 
  prediction.assay = T,
  weight.reduction = mdxTG[["pca"]], 
  dims = 1:30
)

mdx[["predictions"]] = predictions.assay

### All groups merged ####
WT <- readRDS("WT_processed.rds")
mdx <- readRDS("mdx_processed.rds")
mdx_TG <- readRDS("mdx_TG_processed.rds")

merged <- merge(
  x = WT,
  y = c(mdx, mdx_TG)
)

DefaultAssay(merged) <- "Spatial"

merged <- NormalizeData(merged)
saveRDS(merged, "merged_data.rds")
merged <- readRDS("merged_data.rds")

SpatialFeaturePlot(
  merged,
  features = "Tgfb1",
  pt.size.factor = 2.8,
  stroke = 1
)

SpatialFeaturePlot(
  merged,
  features = "Lama2"
) &
  scale_fill_viridis_c(limits = c(0,3.5))

## mechano pathways ####
merged <- readRDS("merged_clusters.rds")
YAP <- readxl::read_excel("YAP target genes.xlsx",
                          sheet = 2)


merged <- AddModuleScore(
  merged,
  features = list(YAP$lower),
  assay = "Spatial",
  name = "YAP_score"
)

min_max <- function(x) {
  ( x - min(x) ) / ( max(x) - min(x) )
}

merged$YAP_scaled <- min_max(merged$YAP1)

quantile(merged$YAP_scaled, probs = seq(0.1, 1, 0.1))

p1 <- SpatialFeaturePlot(merged, features = "YAP_scaled",
                         image.alpha = 0,
                         alpha = c(1,1),
                         min.cutoff = 0.25,
                         max.cutoff = 0.7,
                         pt.size.factor = 3,
                         stroke = 0.1) &
  NoLegend()
# viridis::scale_fill_viridis() &
# scale_fill_viridis_c(limits = c(0,1))

p1 <- SpatialFeaturePlot(merged, features = "Cd44",
                         image.alpha = 0,
                         pt.size.factor = 3,
                         stroke = 0.1) &
  NoLegend() 
  # scale_fill_continuous(
  #   palette = rev(RColorBrewer::brewer.pal(name = "RdYlBu",
  #                                     n = 11)),
  #   limits = c(0,3.5))

p1_raster <- list()
for (i in seq(p1)) {
  p1_raster[[i]] <- ggrastr::rasterise(p1[[i]], 
                                       layers = "Spatial", 
                                       dpi = 800)
}

p1_raster[[1]] | p1_raster[[2]] | p1_raster[[3]]

ggsave(file = "Figures/Cd44_raster.svg",
       height = 10, width = 10,
       units = "cm")


merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged)
merged <- RunUMAP(merged, dims = 1:30)
merged <- FindNeighbors(merged)
merged <- FindClusters(
  merged,
  algorithm = 4,
  resolution = 0.8
)

saveRDS(object = merged, file = "merged_clusters.rds")

### immune pathways ####
# genes related to immune fibrosis axis
genes <- c('Lgals3', 'Fabp5', 'Spp1', 
          'Trem2', 'Thbs1', 'Gpnmb') # pro-fibrotic macrophage markers

merged <- AddModuleScore(
  merged,
  features = list(genes),
  assay = "Spatial",
  name = "fibro_score"
)

merged$fibro_scaled <- min_max(merged$fibro_score1)

quantile(merged$fibro_scaled, probs = seq(0.1, 1, 0.1))

p1 <- SpatialFeaturePlot(merged, features = "fibro_scaled",
                         image.alpha = 0,
                         min.cutoff = 0.05,
                         # max.cutoff = 0.7,
                         pt.size.factor = 3,
                         stroke = 0.1) &
  NoLegend()

p1_raster <- list()
for (i in seq(p1)) {
  p1_raster[[i]] <- ggrastr::rasterise(p1[[i]], 
                                       layers = "Spatial", 
                                       dpi = 800)
}

p1_raster[[1]] | p1_raster[[2]] | p1_raster[[3]]

ggsave(file = "Figures/pro_fibro_raster.svg",
       height = 20, width = 20,
       units = "cm")

## Area and Expression Quantification##
WT <- readRDS("WT_processed.rds")
mdx <- readRDS("mdx_processed.rds")
mdx_TG <- readRDS("mdx_TG_processed.rds")

# add in genotype column
WT_counts <- ncol(WT@assays$Spatial$counts)
mdx_counts <- ncol(mdx@assays$Spatial$counts)
mdxTG_counts <- ncol(mdx_TG@assays$Spatial$counts)

Genotype <- c(rep("WT", WT_counts), rep("mdx", mdx_counts), rep("mdx_TG", mdxTG_counts))
merged$Genotype <- Genotype
head(merged@meta.data,3)

data <- DotPlot(
  merged,
  features = c("Col1a1", "Col3a1", "Col5a1", "Yap1", "Tgfb1"),
  group.by = "Genotype"
)


df <- as.data.frame(data$data)
row.names(df) <- NULL
df <- df %>% dplyr::relocate(features.plot)
View(df)

# multiple pathways with corrleation ----
YAP_targets <- readxl::read_xlsx('YAP target genes.xlsx', sheet = 2)
gene_sets <- list(
  
  YAP = YAP_targets$lower,
  
  Inflammatory = c(
    'Lgals3', 'Fabp5', 'Spp1', 'Ctss',
    'Trem2', 'Thbs1', 'Gpnmb', 'Cd68',
    'Mmp19', 'Adgre1', 'Stat3', 'Lgals1', 'Il2rb',
    'Ifngr1', 'S100a6', 'Ptprc'
    ),
  
  Fibrotic = c(
    'Col1a1', 'Col1a2', 'Col3a1', 'Fn1', 
    'Tgfb1', 'Ccn2', 'Mmp3', 'Tnc', 'Bgn',
    'Timp1', 'Timp2', 'Pdgfra',
    'Pdgfrb', 'Postn', 'Cd44', 'Thbs4'
    )
  )

merged <- AddModuleScore(object = merged,
                         features = gene_sets,
                         name = names(gene_sets))


SpatialFeaturePlot(merged, features = "Fibrotic3",
                   image.alpha = 0,
                   alpha = c(1,1),
                   pt.size.factor = 3,
                   stroke = 0.1) &
  NoLegend()

#### Gene expression based on spots ####
merged <- readRDS("merged.rds")

## Define Col1 positive and negative spots
# try changing expression cutoff for more stringent results
Yap1_row <- which(rownames(merged) == "Yap1")
Ccn2_row <- which(rownames(merged) == "Ccn2")

# shows value of 1 is >70th percentile for Col1a1 and >50th percentile for Col1a2
quantile(merged@assays$Spatial$data[Yap1_row,])
quantile(merged@assays$Spatial$data[Ccn2_row,])
ecdf(merged@assays$Spatial$data[Yap1_row,])(1)
ecdf(merged@assays$Spatial$data[Ccn2_row,])(1)

# add Col1 ID data to the cells based on thresholds
Idents(
  merged,
  WhichCells(object = merged,
             expression = Yap1 > 1,
             slot = 'data')
) <- "YAP_pos"

Idents(
  merged,
  WhichCells(object = merged,
             expression = Yap1 <= 1,
             slot = 'data')
) <- "YAP_neg"

# DE general pos vs negative
# not super useful
genes <- FindMarkers(
  merged,
  ident.1 = 'YAP_pos',
  ident.2 = 'YAP_neg'
)
genes <- genes %>% dplyr::filter(p_val_adj < 0.05)
genes$gene <- row.names(genes)
row.names(genes) <- NULL
genes <- genes %>% relocate(gene)
View(genes)

# add meta data column to ID genotype and collagen
merged$Genotype_YAP <- paste0(
  merged$Genotype,
  "_",
  merged@active.ident
)
head(merged@meta.data,3)
merged$Genotype_YAP <- as.factor(merged$Genotype_YAP)
merged@active.ident <- merged$Genotype_YAP
saveRDS(merged, "merged.rds")

DEenrichRPlot(
  merged,
  ident.1 = "mdx_YAP_pos",
  ident.2 = "mdx_YAP_neg",
  balanced = F,
  max.genes = 1000,
  enrich.database = "GO_Cellular_Component_2023"
)

pathway_list <- DEenrichRPlot(
  merged,
  ident.1 = "mdx_TG_YAP_pos",
  ident.2 = "WT_YAP_pos",
  balanced = F,
  max.genes = 1000,
  enrich.database = "GO_Cellular_Component_2023",
  return.gene.list = T
)

pathway_list <- as.data.frame(pathway_list$pos)
collagen <- pathway_list$GO_Cellular_Component_2023.Genes[1]
focal_adhesion <- pathway_list$GO_Cellular_Component_2023.Genes[3]

### BayesSpace Analysis ####
merged <- readRDS("merged.rds")
head(merged@meta.data,3)

#### WT ####
WT <- subset(
  merged,
  Genotype == "WT"
)

WT_prep <- DietSeurat(WT)
WT_sce <- as.SingleCellExperiment(WT_prep)
colData(WT_sce) <- cbind(colData(WT_sce), WT_prep@images$slice1@coordinates)

set.seed(001)
WT_sce <- spatialPreprocess(
  WT_sce, 
  platform = "Visium",
  skip.PCA = F,
  log.normalize = F
  )

WT_sce <- qTune(WT_sce, qs=seq(2, 10), platform="Visium", d=15)
qPlot(WT_sce)
# 5 clusters looks ok

set.seed(002)
WT_sce <- spatialCluster(
  WT_sce, 
  q = 5, 
  platform = "Visium",
  d = 15,
  init.method = "mclust",
  model="t", 
  nrep = 10000
  )

head(colData(WT_sce))

# color pal for visualization 
colors_new <- c(
  "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", 
  "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", 
  "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", 
  "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", 
  "yellow4", "yellow3", "darkorange4", "brown"
)

clusterPlot(WT_sce, palette = colors_new) &
  theme(legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(face = "bold", size = 15))

saveRDS(WT_sce, "WT_BayesSpace.rds")
WT_sce <- readRDS("WT_BayesSpace.rds")

# tune jitter_scale parameter to result in chain stabilization 0.25-0.40
WT_enhance <- spatialEnhance(
  WT_sce, 
  q = 5, 
  platform = "Visium", 
  d = 15,
  model = "t", 
  jitter_scale = 1,
  nrep=100000,
  save.chain = T)
mcmcChain(WT_enhance, "Ychange")

clusterPlot(WT_sce, palette = colors_new) &
  theme(legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(face = "bold", size = 15))

saveRDS(WT_enhance, "WT_EnhanceBayesSpace.rds")

#### mdx ####
mdx <- subset(
  merged,
  Genotype == "mdx"
)

mdx_prep <- DietSeurat(mdx)
mdx_sce <- as.SingleCellExperiment(mdx_prep)
colData(mdx_sce) <- cbind(colData(mdx_sce), mdx_prep@images$slice1@coordinates)

set.seed(003)
mdx_sce <- spatialPreprocess(
  mdx_sce, 
  platform = "Visium",
  skip.PCA = F,
  log.normalize = F
)

mdx_sce <- qTune(mdx_sce, qs=seq(2, 10), platform = "Visium", d = 15)
qPlot(mdx_sce)
# 5 clusters looks ok

set.seed(004)
mdx_sce <- spatialCluster(
  mdx_sce, 
  q = 5, 
  platform = "Visium",
  d = 15,
  init.method = "mclust",
  model="t", 
  nrep = 10000
)

head(colData(mdx_sce))

clusterPlot(mdx_sce, palette = colors_new) &
  theme(legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(face = "bold", size = 15))

saveRDS(mdx_sce, "mdx_BayesSpace.rds")
mdx_sce <- readRDS("mdx_BayesSpace.rds")

# tune jitter_scale parameter to result in chain stabilization 0.25-0.40
mdx_enhance <- spatialEnhance(
  mdx_sce, 
  q = 5, 
  platform = "Visium", 
  d = 15,
  model = "t", 
  jitter_scale = 1.6,
  nrep = 100000
  )
mcmcChain(mdx_enhance, "Ychange")
# 1.6 looks fine with nrep = 1800, stabilizes to ~0.38
# proceeded to run with nrep = 100,000
# did not save chain (potentially saves time?)
# total runtime was 

clusterPlot(mdx_enhance, palette = colors_new) &
  theme(legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(face = "bold", size = 15))

saveRDS(mdx_enhance, "mdx_EnhanceBayesSpace.rds")

## enhanced features
markers <- c("Yap1", "Ccn2", "Col1a1", "Col1a2", "Pdgfra", "Cd68")
?enhanceFeatures
mdx_enhance <- enhanceFeatures(
  mdx_enhance, 
  mdx_sce,
  feature_names = markers,
  nrounds = 0
  )

?featurePlot
p1 <- featurePlot(mdx_enhance, "Yap1") &
  theme(legend.text = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 15)) &
  viridis::scale_fill_viridis()

p2 <- featurePlot(mdx_enhance, "Pdgfra") &
  theme(legend.text = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 15)) &
  viridis::scale_fill_viridis()

p3 <- featurePlot(mdx_enhance, "Cd68") &
  theme(legend.text = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 15)) &
  viridis::scale_fill_viridis()

p4 <- featurePlot(mdx_enhance, "Ccn2") &
  theme(legend.text = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 15)) &
  viridis::scale_fill_viridis()

p5 <- featurePlot(mdx_enhance, "Col1a1") &
  theme(legend.text = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 15)) &
  viridis::scale_fill_viridis()

p6 <- featurePlot(mdx_enhance, "Col1a2") &
  theme(legend.text = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 15)) &
  viridis::scale_fill_viridis()

(p1|p2|p3)/(p4|p5|p6)