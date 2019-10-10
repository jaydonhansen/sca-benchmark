library(dplyr)
library(Seurat)
library(SingleCellExperiment)

setwd("~/advanced-bioinformatics/")
load("data/sincell_with_class.RData")

sce_10x_seurat <- sce_sc_10x_qc %>% as.Seurat(counts = "counts", data = "logcounts") %>%
    FindVariableFeatures() %>% ScaleData() %>% RunPCA()

p1 <- DimPlot(sce_10x_seurat, reduction = "pca", group.by = "cell_line")

p1
