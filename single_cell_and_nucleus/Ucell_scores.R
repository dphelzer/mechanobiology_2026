pak::pak("cellgeni/sceasy")
library(sceasy)
library(Seurat)
library(ggplot2)
library(patchwork)
## get correct FAPs from single cell data
FAPs <- convertFormat(
  "single_cell_FAPs.h5ad", 
  from = "anndata", 
  to = "seurat", 
  outFile = NULL)

# FAPs UCell scoring ------------------------------------------------------
# BiocManager::install("UCell")
library(UCell)

rownames(FAPs) # all capital

# YAP_targets <- readxl::read_xlsx('YAP target genes.xlsx', sheet = 2)

gene_sets <- list(
  
  ## from PMID: 29317486, PMID: 19324877, PMID: 30380420
  YAP = c(
    'ANKRD1', 'CCN1', 'CCN2', 'CDC42EP3', 'PARVA', 'NEDD9', 
    'CSRP1', 'FMN2', 'THBS1', 'DCN', 'SPARC', 'COL6A2', 'MFAP5', 
    'LTBP1', 'PXDN', 'INHBA', 'TGFB2', 'CDH13', 'DAB2', 'MYOF', 
    'ASAP1', 'CCDC80', 'DOCK5', 'ARHGEF17', 'NUAK2', 'LATS2', 
    'AMOTL2'
    ),
  
  Inflammatory = c(
    'LGALS3', 'FABP5', 'SPP1', 'CTSS',
    'TREM2', 'THBS1', 'GPNMB', 'CD68',
    'MMP19', 'ADGRE1', 'STAT3', 'LGALS1', 'IL2RB',
    'IFNGR1', 'S100A6', 'PTPRC'
  ),
  
  Fibrotic = c(
    'COL1A1', 'COL1A2', 'COL3A1', 'FN1', 
    'TGFB1', 'CCN2', 'MMP3', 'TNC', 'BGN',
    'TIMP1', 'TIMP2', 'PDGFRA',
    'PDGFRB', 'POSTN', 'CD44', 'THBS4'
  )
)

FAPs@active.ident <- FAPs$condition

FAPs <- AddModuleScore(FAPs, gene_sets)

FAPs <- AddModuleScore_UCell(
  FAPs, 
  features=gene_sets,
  name = NULL)

head(FAPs@meta.data)

## plot UCell as vln plot ----
VlnPlot(FAPs, features = "YAP", 
        group.by = 'condition',
        pt.size = 0) +
  scale_fill_manual(values = c('black', 'red', 'blue')) +
  NoLegend() +
  labs(x = NULL,
       y = 'Score',
       title = NULL) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.35, colour = 'black'),
    axis.ticks = element_line(linewidth = 0.35, colour = 'black'),
    axis.text.y = element_text(colour = 'black', size = 8),
    axis.text.x = element_text(colour = 'black', size = 8, 
                               angle = 0, hjust = 0.5),
    axis.title = element_text(colour = 'black', size = 8),
    title = element_text(colour = 'black', size = 8),
    legend.title = element_text(colour = 'black', size = 8),
    legend.text = element_text(colour = 'black', size = 8),
    legend.ticks = element_line(linewidth = 0.35, colour = 'black'),
    axis.ticks.length.y = unit(-0.1, 'cm'),
    axis.ticks.length.x = unit(0, 'cm')
  )

ggsave("Figures/YAP_score_FAPs_sc.svg",
       height = 3.5, width = 4.5, units = 'cm')

FeatureScatter(FAPs, "YAP", "Fibrotic",
               split.by = "condition")

saveRDS(FAPs, "seurat_FAPs_single_cell.rds")

# yap vs fibrotic score plots ---------------------------------------------
FAPs <- readRDS('seurat_FAPs_single_cell.rds')
