library(scran)
library(scater)
library(dplyr)
library(Seurat)
library(SingleCellExperiment)

setwd("~/advanced-bioinformatics/")
load("data/sincell_with_class.RData")
load("data/sincell_with_class_5cl.RData")

# Pipe the 10x datasets through the regular PCA/t-SNE pipeline in Seurat
seurat_10x_3cl <- sce_sc_10x_qc %>% 
    as.Seurat(counts = "counts", data = "logcounts") %>%
    FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% FindNeighbors() %>% FindClusters(reduction.type = "pca", dims.use = 1:5) %>%
    RunTSNE %>% RunUMAP(dims = 1:5)
colData(sce_sc_10x_qc)$clustering_res <- as.factor(seurat_10x_3cl@active.ident)

# Pipe the 10x datasets through the regular PCA/t-SNE pipeline in Seurat
seurat_10x_5cl <- sce_sc_10x_5cl_qc %>% 
    as.Seurat(counts = "counts", data = "logcounts") %>%
    FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% FindNeighbours() %>% FindClusters(reduction.type = "pca") %>% RunTSNE %>% RunUMAP(dims = 1:5)
colData(sce_sc_10x_5cl_qc)$clustering_res <- as.factor(seurat_10x_5cl@active.ident)

# Pipe the CELLseq2 datasets through the regular PCA/t-SNE pipeline in Seurat
seurat_CELseq2_3cl <- sce_sc_CELseq2_qc %>%
    as.Seurat(counts = "counts", data = "logcounts") %>%
    FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% FindNeighbors() %>% FindClusters(reduction.type = "pca") %>% RunTSNE %>% RunUMAP(dims = 1:5)
colData(sce_sc_CELseq2_qc)$clustering_res <- as.factor(seurat_CELseq2_3cl@active.ident)


seurat_CELseq2_5cl_p1 <- sc_Celseq2_5cl_p1 %>% 
    as.Seurat(counts = "counts", data = "logcounts") %>%
    FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% FindNeighbours() %>% FindClusters(reduction.type = "pca") %>% RunTSNE %>% RunUMAP(dims = 1:5)
colData(sc_Celseq2_5cl_p1)$clustering_res <- as.factor(seurat_CELseq2_5cl_p1@active.ident)

seurat_CELseq2_5cl_p2 <- sc_Celseq2_5cl_p2 %>% 
    as.Seurat(counts = "counts", data = "logcounts") %>%
    FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% FindNeighbours() %>% FindClusters(reduction.type = "pca") %>% RunTSNE %>% RunUMAP(dims = 1:5)
colData(sc_Celseq2_5cl_p2)$clustering_res <- as.factor(seurat_CELseq2_5cl_p2@active.ident)

seurat_CELseq2_5cl_p3 <- sc_Celseq2_5cl_p3 %>% 
    as.Seurat(counts = "counts", data = "logcounts") %>%
    FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% FindNeighbours() %>% FindClusters(reduction.type = "pca") %>% RunTSNE %>% RunUMAP(dims = 1:5)
colData(sc_Celseq2_5cl_p3)$clustering_res <- as.factor(seurat_CELseq2_5cl_p3@active.ident)
 
# Pipe the Dropseq datasets through the regular PCA/t-SNE pipeline in Seurat
seurat_Dropseq_3cl <- sce_sc_Dropseq_qc %>% 
    as.Seurat(counts = "counts", data = "logcounts") %>%
    FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% FindNeighbours() %>% FindClusters(reduction.type = "pca") %>% RunTSNE %>% RunUMAP(dims = 1:5)
colData(sce_sc_Dropseq_qc)$clustering_res <- as.factor(seurat_Dropseq_3cl@active.ident)

# PCA/t-SNE plots for 10x
pca_plot_10x_3cl <- DimPlot(seurat_10x_3cl, reduction = "pca", group.by = "cell_line")
pca_plot_10x_5cl <- DimPlot(seurat_10x_5cl, reduction = "pca", group.by = "cell_line")

tsne_plot_10x_3cl <- DimPlot(seurat_10x_3cl, reduction = "tsne", group.by = "cell_line")
tsne_plot_10x_5cl <- DimPlot(seurat_10x_5cl, reduction = "tsne", group.by = "cell_line")

# UMAP plots for 10x
umap_plot_10x_3cl <- DimPlot(seurat_10x_3cl, reduction = "umap", group.by = "cell_line")
umap_plot_10x_5cl <- DimPlot(seurat_10x_5cl, reduction = "umap", group.by = "cell_line")

# PCA/t-SNE plots for CELseq2
pca_plot_CELseq2_3cl <- DimPlot(seurat_CELseq2_3cl, reduction = "pca", group.by = "cell_line")
pca_plot_CELseq2_5cl_p1 <- DimPlot(seurat_CELseq2_5cl_p1, reduction = "pca", group.by = "cell_line")
pca_plot_CELseq2_5cl_p2 <- DimPlot(seurat_CELseq2_5cl_p2, reduction = "pca", group.by = "cell_line")
pca_plot_CELseq2_5cl_p3 <- DimPlot(seurat_CELseq2_5cl_p3, reduction = "pca", group.by = "cell_line")

tsne_plot_CELseq2_3cl <- DimPlot(seurat_CELseq2_3cl, reduction = "tsne", group.by = "cell_line")
tsne_plot_CELseq2_5cl_p1 <- DimPlot(seurat_CELseq2_5cl_p1, reduction = "tsne", group.by = "cell_line")
tsne_plot_CELseq2_5cl_p2 <- DimPlot(seurat_CELseq2_5cl_p2, reduction = "tsne", group.by = "cell_line")
tsne_plot_CELseq2_5cl_p3 <- DimPlot(seurat_CELseq2_5cl_p3, reduction = "tsne", group.by = "cell_line")

pca_plot_Dropseq <- DimPlot(seurat_Dropseq_3cl, reduction = "pca", group.by = "cell_line")
tsne_plot_Dropseq <- DimPlot(seurat_Dropseq_3cl, reduction = "tsne", group.by = "cell_line")


# 10x samples

png("output/seurat/seurat-pca-3cl-10x.png", height = 1000, width = 1000)
    pca_plot_10x_3cl # 3-cell PCA
dev.off()

png("output/seurat/seurat-tsne-3cl-10x.png", height = 1000, width = 1000)
    tsne_plot_10x_3cl # 3-cell t-SNE
dev.off()

png("output/seurat/seurat-umap-3cl-10x.png", height = 1000, width = 1000)
    umap_plot_10x_3cl # 3-cell UMAP
dev.off()

png("output/seurat/seurat-pca-5cl-10x.png", height = 1000, width = 1000)
    pca_plot_10x_5cl # 5-cell PCA
dev.off()

png("output/seurat/seurat-tsne-5cl-10x.png", height = 1000, width = 1000)
    tsne_plot_10x_5cl # 5-cell t-SNE
dev.off()

png("output/seurat/seurat-umap-5cl-10x.png", height = 1000, width = 1000)
    umap_plot_10x_5cl # 5-cell UMAP
dev.off()


# CELseq2 samples

png("output/seurat/seurat-pca-3cl-CELseq2.png", height = 1000, width = 1000)
    pca_plot_CELseq2_3cl # 3-cell PCA
dev.off()

png("output/seurat/seurat-tsne-3cl-CELseq2.png", height = 1000, width = 1000)
    tsne_plot_CELseq2_3cl # 3-cell t-SNE
dev.off()

png("output/seurat/seurat-pca-5cl-p1-CELseq2.png", height = 1000, width = 1000)
    pca_plot_CELseq2_5cl_p1 # 5-cell p1 PCA
dev.off()

png("output/seurat/seurat-tsne-5cl-p1-CELseq2.png", height = 1000, width = 1000)
    tsne_plot_CELseq2_5cl_p1 # 5-cell p1 t-SNE
dev.off()

png("output/seurat/seurat-pca-5cl-p2-CELseq2.png", height = 1000, width = 1000)
    pca_plot_CELseq2_5cl_p2 # 5-cell p2 PCA
dev.off()

png("output/seurat/seurat-tsne-5cl-p2-CELseq2.png", height = 1000, width = 1000)
    tsne_plot_CELseq2_5cl_p2 # 5-cell p2 t-SNE
dev.off()

png("output/seurat/seurat-pca-5cl-p3-CELseq2.png", height = 1000, width = 1000)
    pca_plot_CELseq2_5cl_p3 # 5-cell p3 PCA
dev.off()

png("output/seurat/seurat-tsne-5cl-p3-CELseq2.png", height = 1000, width = 1000)
    tsne_plot_CELseq2_5cl_p3 # 5-cell p3 t-SNE
dev.off()

# Dropseq samples

png("output/seurat/seurat-pca-3cl-Dropseq.png", height = 1000, width = 1000)
    pca_plot_Dropseq # 3-cell PCA
dev.off()

png("output/seurat/seurat-tsne-3cl-Dropseq.png", height = 1000, width = 1000)
    tsne_plot_Dropseq # 3-cell t-SNE
dev.off()
