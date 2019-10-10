library(dplyr)
library(Seurat)
library(SingleCellExperiment)

setwd("~/advanced-bioinformatics/")
load("data/sincell_with_class.RData")

seurat_10x_qc <- sce_sc_10x_qc %>% as.Seurat(counts = "counts", data = "logcounts") %>%
    FindVariableFeatures() %>% ScaleData() %>% RunPCA()

seurat_CELseq2_qc <- sce_sc_CELseq2_qc %>% as.Seurat(counts = "counts", data = "logcounts") %>%
    FindVariableFeatures() %>% ScaleData() %>% RunPCA()

seurat_Dropseq_qc <- sce_sc_Dropseq_qc %>% as.Seurat(counts = "counts", data = "logcounts") %>%
    FindVariableFeatures() %>% ScaleData() %>% RunPCA()

p1 <- DimPlot(seurat_10x_qc, reduction = "pca", group.by = "cell_line")
p2 <- DimPlot(seurat_CELseq2_qc, reduction = "pca", group.by = "cell_line")
p3 <- DimPlot(seurat_Dropseq_qc, reduction = "pca", group.by = "cell_line")

p1
p2
p3
