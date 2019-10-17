#Analysis of Three Cell Line Mixtures using Monocle3

library(monocle3)
library(dplyr)
library(devtools)

setwd("~/UQ/Semester 2 2019/Adv. Bioinformatics/Project")
x = load("sincell_with_class.RData")
sincell_3cl = get(x)
rm(x)

# 10X Chromium Mixture

sce10x_counts <- sce_sc_10x_qc@assays$data$counts
sce10x_cell_metadata <- colData(sce_sc_10x_qc)
sce10x_gene_metadata <- rowData(sce_sc_10x_qc)

cds_10x <- new_cell_data_set(sce10x_counts,
                             cell_metadata = sce10x_cell_metadata,
                             gene_metadata = sce10x_gene_metadata)

cds_10x <- preprocess_cds(cds_10x, num_dim = 100)

# 10X UMAP

cds_10x <- reduce_dimension(cds_10x, reduction_method = "UMAP")

cds_10x = cluster_cells(cds_10x)
learn_graph(cds_10x)

umap_3cl_10x <- plot_cells(cds_10x, color_cells_by = "cell_line", group_label_size = 3.5)

png("monocle3-umap-3cl-10x.png", height = 1000, width = 1000)
umap_3cl_10x
dev.off()

# 10X tSNE

cds_10x_tSNE <- reduce_dimension(cds_10x, reduction_method = "tSNE")

cds_10x_tSNE = cluster_cells(cds_10x_tSNE)
learn_graph(cds_10x_tSNE)

tsne_3cl_10x <-plot_cells(cds_10x_tSNE, reduction_method = "tSNE", color_cells_by = "cell_line", group_label_size = 3.5)

png("monocle3-tsne-3cl-10x.png", height = 1000, width = 1000)
tsne_3cl_10x
dev.off()

# 10X PCA

cds_10x_PCA <- reduce_dimension(cds_10x, reduction_method = "PCA")

cds_10x_PCA = cluster_cells(cds_10x_PCA)
learn_graph(cds_10x_PCA)

pca_3cl_10x <- plot_cells(cds_10x_PCA, reduction_method = "PCA", color_cells_by = "cell_line", group_label_size = 3.5)

png("monocle3-pca-3cl-10x.png", height = 1000, width = 1000)
pca_3cl_10x
dev.off()



# CEL-seq2 Mixture

celseq2_counts <- sce_sc_CELseq2_qc@assays$data$counts
celseq2_cell_metadata <- colData(sce_sc_CELseq2_qc)
celseq2_gene_metadata <- rowData(sce_sc_CELseq2_qc)

cds_celseq2 <- new_cell_data_set(celseq2_counts,
                             cell_metadata = celseq2_cell_metadata,
                             gene_metadata = celseq2_gene_metadata)

cds_celseq2 <- preprocess_cds(cds_celseq2, num_dim = 100)

# CEL-seq2 UMAP

cds_celseq2 <- reduce_dimension(cds_celseq2, reduction_method = "UMAP")

cds_celseq2 = cluster_cells(cds_celseq2)
learn_graph(cds_celseq2)

umap_3cl_celseq2 <- plot_cells(cds_celseq2, color_cells_by = "cell_line", group_label_size = 3.5)

png("monocle3-umap-3cl-celseq2.png", height = 1000, width = 1000)
umap_3cl_celseq2
dev.off()

# CEL-seq2 tSNE

cds_celseq2_tSNE <- reduce_dimension(cds_celseq2, reduction_method = "tSNE")

cds_celseq2_tSNE = cluster_cells(cds_celseq2_tSNE)
learn_graph(cds_celseq2_tSNE)

tsne_3cl_celseq2 <- plot_cells(cds_celseq2_tSNE, reduction_method = "tSNE", color_cells_by = "cell_line", group_label_size = 3.5)

png("monocle3-tsne-3cl-celseq2.png", height = 1000, width = 1000)
tsne_3cl_celseq2
dev.off()

# CEL-seq2 PCA

cds_celseq2_PCA <- reduce_dimension(cds_celseq2, reduction_method = "PCA")

cds_celseq2_PCA = cluster_cells(cds_celseq2_PCA)
learn_graph(cds_celseq2_PCA)

pca_3cl_celseq2 <- plot_cells(cds_celseq2_PCA, reduction_method = "PCA", color_cells_by = "cell_line", group_label_size = 3.5)

png("monocle3-pca-3cl-celseq2.png", height = 1000, width = 1000)
pca_3cl_celseq2
dev.off()



# Drop-seq Mixture

dropseq_counts <- sce_sc_Dropseq_qc@assays$data$counts
dropseq_cell_metadata <- colData(sce_sc_Dropseq_qc)
dropseq_gene_metadata <- rowData(sce_sc_Dropseq_qc)

cds_dropseq <- new_cell_data_set(dropseq_counts,
                                 cell_metadata = dropseq_cell_metadata,
                                 gene_metadata = dropseq_gene_metadata)

cds_dropseq <- preprocess_cds(cds_dropseq, num_dim = 100)

# Drop-seq UMAP

cds_dropseq <- reduce_dimension(cds_dropseq, reduction_method = "UMAP")

cds_dropseq = cluster_cells(cds_dropseq)
learn_graph(cds_dropseq)

umap_3cl_dropseq <- plot_cells(cds_dropseq, color_cells_by = "cell_line", group_label_size = 3.5)

png("monocle3-umap-3cl-dropseq.png", height = 1000, width = 1000)
umap_3cl_dropseq
dev.off()

# Drop-seq tSNE

cds_dropseq_tSNE <- reduce_dimension(cds_dropseq, reduction_method = "tSNE")

cds_dropseq_tSNE = cluster_cells(cds_dropseq_tSNE)
learn_graph(cds_dropseq_tSNE)

tsne_3cl_dropseq <- plot_cells(cds_dropseq_tSNE, reduction_method = "tSNE", color_cells_by = "cell_line", group_label_size = 3.5)

png("monocle3-tsne-3cl-dropseq.png", height = 1000, width = 1000)
tsne_3cl_dropseq
dev.off()

# Drop-seq PCA

cds_dropseq_PCA <- reduce_dimension(cds_dropseq, reduction_method = "PCA")

cds_dropseq_PCA = cluster_cells(cds_dropseq_PCA)
learn_graph(cds_dropseq_PCA)

pca_3cl_dropseq <- plot_cells(cds_dropseq_PCA, reduction_method = "PCA", color_cells_by = "cell_line", group_label_size = 3.5)

png("monocle3-pca-3cl-dropseq.png", height = 1000, width = 1000)
pca_3cl_dropseq
dev.off()
