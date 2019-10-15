library(dplyr)
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(phyclust)

setwd("~/advanced-bioinformatics/")
load("data/sincell_with_class.RData")

#example code
sc_Celseq2_5cl_p1@assays$data$logcounts <- sc_Celseq2_5cl_p1@assays$data$counts
# Output the plots into images
sc_Celseq2_5cl_p1 <- runUMAP(sc_Celseq2_5cl_p1, n_neighbors=15, pca=NULL)
sc_Celseq2_5cl_p1 <- runTSNE(sc_Celseq2_5cl_p1)
# plot UMAP/TSNE/PCA images
plot(plotUMAP(sc_Celseq2_5cl_p1))
plot(plotTSNE(sc_Celseq2_5cl_p1)) 
plot(plotPCA(sc_Celseq2_5cl_p1))