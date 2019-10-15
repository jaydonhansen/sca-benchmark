library(dplyr)
library(Seurat)
library(SingleCellExperiment)
library(scater)

setwd("~/advanced-bioinformatics/")
load("data/sincell_with_class_5cl.RData")

# generate log counts (not sure that data is QC'd properly still)
sc_Celseq2_5cl_p1@assays$data$logcounts <- sc_Celseq2_5cl_p1@assays$data$counts

# coax into a seurat object
seurat_Celseq2_5cl_p1 <- sc_Celseq2_5cl_p1 %>% 
    as.Seurat(counts = "counts", data = "logcounts") %>%
    FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunTSNE

# constructing plots
pca_plot_Celseq2_5cl_p1 <- DimPlot(seurat_Celseq2_5cl_p1, reduction = "pca", group.by = "cell_line_demuxlet")
tsne_plot_Celseq2_5cl_p1 <- DimPlot(seurat_Celseq2_5cl_p1, reduction = "tsne", group.by = "cell_line_demuxlet")

pca_plot_Celseq2_5cl_p1
tsne_plot_Celseq2_5cl_p1
