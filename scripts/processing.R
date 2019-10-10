library(dplyr)
library(Seurat)
library(SingleCellExperiment)

setwd("~/advanced-bioinformatics/")
load("data/sincell_with_class.RData")

# Pipe the 10x dataset through the regular PCA pipeline in Seurat
seurat_10x_qc <- sce_sc_10x_qc %>% as.Seurat(counts = "counts", data = "logcounts") %>%
    FindVariableFeatures() %>% ScaleData() %>% RunPCA()

rm(sce_sc_10x_qc) # Save memory

# Pipe the CELLseq2 dataset through the regular PCA pipeline in Seurat
seurat_CELseq2_qc <- sce_sc_CELseq2_qc %>% as.Seurat(counts = "counts", data = "logcounts") %>%
    FindVariableFeatures() %>% ScaleData() %>% RunPCA()

rm(sce_sc_CELseq2_qc) # Save memory

# Pipe the Dropseq dataset through the regular PCA pipeline in Seurat
seurat_Dropseq_qc <- sce_sc_Dropseq_qc %>% as.Seurat(counts = "counts", data = "logcounts") %>%
    FindVariableFeatures() %>% ScaleData() %>% RunPCA()

rm(sce_sc_Dropseq_qc) # Save memory

# Some PCA plots
p1 <- DimPlot(seurat_10x_qc, reduction = "pca", group.by = "cell_line")
p2 <- DimPlot(seurat_CELseq2_qc, reduction = "pca", group.by = "cell_line")
p3 <- DimPlot(seurat_Dropseq_qc, reduction = "pca", group.by = "cell_line")

p1
p2
p3
