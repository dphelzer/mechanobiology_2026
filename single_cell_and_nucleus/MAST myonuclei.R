# load in anndata to sce
ad <- reticulate::import("anndata", convert = FALSE)
anndata <- ad$read_h5ad("myonuclei_var_genes.h5ad")
myonuclei <- zellkonverter::AnnData2SCE(anndata)


# run MAST model
library(MAST)

## gene detection rate
colData(myonuclei)$cngeneson <- scale(colSums(assay(myonuclei) > 0))
head(myonuclei@colData)

## conver to sca for input to MAST
myonuclei_assay <- SceToSingleCellAssay(myonuclei)

## run model
fit <- zlm(~ comparisons + cngeneson, myonuclei_assay)

# predict for expected gene expression values
MM <- model.matrix(
  ~comparisons, 
  unique(colData(myonuclei_assay)[,c("comparisons"),drop=FALSE])
  )

predicted <- predict(fit, modelmatrix = MM)
head(predicted)

predicted_mat <- as.matrix(predicted)

# take mean across cells per condition

# find correlation between genes