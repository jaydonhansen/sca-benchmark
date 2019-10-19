#Analysis of Three Cell Line Mixtures using Monocle3

library(monocle3)
library(dplyr)
library(devtools)

setwd("~/UQ/Semester 2 2019/Adv. Bioinformatics/Project")
x = load("sincell_with_class_5cl.RData")
sincell_5cl = get(x)
rm(x)



# 10X 5cl Chromium Mixture (sc_10x_5cl)

sce10x_5cl_counts <- sce_sc_10x_5cl_qc@assays$data$counts
sce10x_5cl_cell_metadata <- colData(sce_sc_10x_5cl_qc)
sce10x_5cl_gene_metadata <- rowData(sce_sc_10x_5cl_qc)

cds_5cl_10x <- new_cell_data_set(sce10x_5cl_counts,
                             cell_metadata = sce10x_5cl_cell_metadata,
                             gene_metadata = sce10x_5cl_gene_metadata)

# Projects data onto top principal components using 100 dimensions, log normalization,
# and scaling
cds_5cl_10x <- preprocess_cds(cds_5cl_10x, num_dim = 100, norm_method = "log", scaling = TRUE)

# 10X 5cl UMAP

# Reduces dimensionality based on chosen method. Default method is UMAP.
cds_5cl_10x <- reduce_dimension(cds_5cl_10x, reduction_method = "UMAP")

# Clusters cells in accordance with identified cell types, and partitions cells as a way to
# detect and handle outliers using a kNN pruning method
cds_5cl_10x <- cluster_cells(cds_5cl_10x, reduction_method = "UMAP")

# Organizes each group into a separate trajectory
cds_5cl_10x <- learn_graph(cds_5cl_10x)

umap_5cl_10x <- plot_cells(cds_5cl_10x, color_cells_by = "cell_line", group_label_size = 3.5,
                           show_trajectory_graph = FALSE)

png("monocle3-umap-5cl-10x.png", height = 1000, width = 1000)
umap_5cl_10x
dev.off()

# 10X 5cl tSNE

cds_5cl_10x_tSNE <- reduce_dimension(cds_5cl_10x, reduction_method = "tSNE")

cds_5cl_10x_tSNE <- cluster_cells(cds_5cl_10x_tSNE, reduction_method = "tSNE")

cds_5cl_10x_tSNE <- learn_graph(cds_5cl_10x_tSNE)

tsne_5cl_10x <- plot_cells(cds_5cl_10x_tSNE, reduction_method = "tSNE", 
                           color_cells_by = "cell_line", group_label_size = 3.5, 
                           show_trajectory_graph = FALSE)

png("monocle3-tsne-5cl-10x.png", height = 1000, width = 1000)
tsne_5cl_10x
dev.off()

# 10X 5cl PCA

cds_5cl_10x_PCA <- reduce_dimension(cds_5cl_10x, reduction_method = "PCA")

cds_5cl_10x_PCA <- cluster_cells(cds_5cl_10x_PCA, reduction_method = "PCA")

cds_5cl_10x_PCA <- learn_graph(cds_5cl_10x_PCA)

pca_5cl_10x <- plot_cells(cds_5cl_10x_PCA, reduction_method = "PCA", 
                          color_cells_by = "cell_line", group_label_size = 3.5, 
                          show_trajectory_graph = FALSE, labels_per_group = 2)

png("monocle3-pca-5cl-10x.png", height = 1000, width = 1000)
pca_5cl_10x
dev.off()



# CEL-seq2 5cl Mixture 1 (sc_CEL-seq2_5cl_p1)

celseq2_p1_counts <- sc_Celseq2_5cl_p1@assays$data$counts
celseq2_p1_cell_metadata <- colData(sc_Celseq2_5cl_p1)
celseq2_p1_gene_metadata <- rowData(sc_Celseq2_5cl_p1)

row.names(celseq2_p1_cell_metadata) = colnames(celseq2_p1_counts)

cds_celseq2_p1 <- new_cell_data_set(celseq2_p1_counts,
                                 cell_metadata = celseq2_p1_cell_metadata,
                                 gene_metadata = celseq2_p1_gene_metadata)

cds_celseq2_p1 <- preprocess_cds(cds_celseq2_p1, num_dim = 100, norm_method = "log", 
                                 scaling = TRUE)

# CEL-seq2 P1 5cl UMAP

cds_celseq2_p1 <- reduce_dimension(cds_celseq2_p1, reduction_method = "UMAP")

cds_celseq2_p1 <- cluster_cells(cds_celseq2_p1, reduction_method = "UMAP")

cds_celseq2_p1 <- learn_graph(cds_celseq2_p1)

umap_celseq2_p1 <- plot_cells(cds_celseq2_p1, color_cells_by = "cell_line_demuxlet", 
                           group_label_size = 3.5, show_trajectory_graph = FALSE)

png("monocle3-umap-celseq2-p1.png", height = 1000, width = 1000)
umap_celseq2_p1
dev.off()

# CEL-seq2 P1 5cl tSNE

cds_celseq2_p1_tSNE <- reduce_dimension(cds_celseq2_p1, reduction_method = "tSNE")

cds_celseq2_p1_tSNE <- cluster_cells(cds_celseq2_p1_tSNE, reduction_method = "tSNE")

cds_celseq2_p1_tSNE <- learn_graph(cds_celseq2_p1_tSNE)

tsne_celseq2_p1 <- plot_cells(cds_celseq2_p1_tSNE, reduction_method = "tSNE", 
                           color_cells_by = "cell_line_demuxlet", group_label_size = 3.5, 
                           show_trajectory_graph = FALSE, labels_per_group = 5)

png("monocle3-tsne-celseq2-p1.png", height = 1000, width = 1000)
tsne_celseq2_p1
dev.off()

# CEL-seq2 P1 5cl PCA

cds_celseq2_p1_PCA <- reduce_dimension(cds_celseq2_p1, reduction_method = "PCA")

cds_celseq2_p1_PCA <- cluster_cells(cds_celseq2_p1_PCA, reduction_method = "PCA")

cds_celseq2_p1_PCA <- learn_graph(cds_celseq2_p1_PCA)

pca_celseq2_p1 <- plot_cells(cds_celseq2_p1_PCA, reduction_method = "PCA", 
                          color_cells_by = "cell_line_demuxlet", group_label_size = 3.5, 
                          show_trajectory_graph = FALSE, labels_per_group = 5)

png("monocle3-pca-celseq-p1.png", height = 1000, width = 1000)
pca_celseq2_p1
dev.off()



# CEL-seq2 5cl Mixture 2 (sc_CEL-seq2_5cl_p2)

celseq2_p2_counts <- sc_Celseq2_5cl_p2@assays$data$counts
celseq2_p2_cell_metadata <- colData(sc_Celseq2_5cl_p2)
celseq2_p2_gene_metadata <- rowData(sc_Celseq2_5cl_p2)

row.names(celseq2_p2_cell_metadata) = colnames(celseq2_p2_counts)

cds_celseq2_p2 <- new_cell_data_set(celseq2_p2_counts,
                                    cell_metadata = celseq2_p2_cell_metadata,
                                    gene_metadata = celseq2_p2_gene_metadata)

cds_celseq2_p2 <- preprocess_cds(cds_celseq2_p2, num_dim = 100, norm_method = "log", 
                                 scaling = TRUE)

# CEL-seq2 P2 5cl UMAP

cds_celseq2_p2 <- reduce_dimension(cds_celseq2_p2, reduction_method = "UMAP")

cds_celseq2_p2 <- cluster_cells(cds_celseq2_p2, reduction_method = "UMAP")

cds_celseq2_p2 <- learn_graph(cds_celseq2_p2)

umap_celseq2_p2 <- plot_cells(cds_celseq2_p2, color_cells_by = "cell_line_demuxlet", 
                              group_label_size = 3.5, show_trajectory_graph = FALSE)

png("monocle3-umap-celseq2-p2.png", height = 1000, width = 1000)
umap_celseq2_p2
dev.off()

# CEL-seq2 P2 5cl tSNE

cds_celseq2_p2_tSNE <- reduce_dimension(cds_celseq2_p2, reduction_method = "tSNE")

cds_celseq2_p2_tSNE <- cluster_cells(cds_celseq2_p2_tSNE, reduction_method = "tSNE")

cds_celseq2_p2_tSNE <- learn_graph(cds_celseq2_p2_tSNE)

tsne_celseq2_p2 <- plot_cells(cds_celseq2_p2_tSNE, reduction_method = "tSNE", 
                              color_cells_by = "cell_line_demuxlet", group_label_size = 3.5, 
                              show_trajectory_graph = FALSE, labels_per_group = 5)

png("monocle3-tsne-celseq2-p2.png", height = 1000, width = 1000)
tsne_celseq2_p2
dev.off()

# CEL-seq2 P2 5cl PCA

cds_celseq2_p2_PCA <- reduce_dimension(cds_celseq2_p2, reduction_method = "PCA")

cds_celseq2_p2_PCA <- cluster_cells(cds_celseq2_p2_PCA, reduction_method = "PCA")

cds_celseq2_p2_PCA <- learn_graph(cds_celseq2_p2_PCA)

pca_celseq2_p2 <- plot_cells(cds_celseq2_p2_PCA, reduction_method = "PCA", 
                             color_cells_by = "cell_line_demuxlet", group_label_size = 3.5, 
                             show_trajectory_graph = FALSE, labels_per_group = 5)

png("monocle3-pca-celseq-p2.png", height = 1000, width = 1000)
pca_celseq2_p2
dev.off()



# CEL-seq2 5cl Mixture 3 (sc_CEL-seq2_5cl_p3)

celseq2_p3_counts <- sc_Celseq2_5cl_p3@assays$data$counts
celseq2_p3_cell_metadata <- colData(sc_Celseq2_5cl_p3)
celseq2_p3_gene_metadata <- rowData(sc_Celseq2_5cl_p3)

row.names(celseq2_p3_cell_metadata) = colnames(celseq2_p3_counts)

cds_celseq2_p3 <- new_cell_data_set(celseq2_p3_counts,
                                    cell_metadata = celseq2_p3_cell_metadata,
                                    gene_metadata = celseq2_p3_gene_metadata)

cds_celseq2_p3 <- preprocess_cds(cds_celseq2_p3, num_dim = 100, norm_method = "log", 
                                 scaling = TRUE)

# CEL-seq2 P3 5cl UMAP

cds_celseq2_p3 <- reduce_dimension(cds_celseq2_p3, reduction_method = "UMAP")

cds_celseq2_p3 <- cluster_cells(cds_celseq2_p3, reduction_method = "UMAP")

cds_celseq2_p3 <- learn_graph(cds_celseq2_p3)

umap_celseq2_p3 <- plot_cells(cds_celseq2_p3, color_cells_by = "cell_line_demuxlet", 
                              group_label_size = 3.5, show_trajectory_graph = FALSE)

png("monocle3-umap-celseq2-p3.png", height = 1000, width = 1000)
umap_celseq2_p3
dev.off()

# CEL-seq2 P3 5cl tSNE

cds_celseq2_p3_tSNE <- reduce_dimension(cds_celseq2_p3, reduction_method = "tSNE")

cds_celseq2_p3_tSNE <- cluster_cells(cds_celseq2_p3_tSNE, reduction_method = "tSNE")

cds_celseq2_p3_tSNE <- learn_graph(cds_celseq2_p3_tSNE)

tsne_celseq2_p3 <- plot_cells(cds_celseq2_p3_tSNE, reduction_method = "tSNE", 
                              color_cells_by = "cell_line_demuxlet", group_label_size = 3.5, 
                              show_trajectory_graph = FALSE, labels_per_group = 4)

png("monocle3-tsne-celseq2-p3.png", height = 1000, width = 1000)
tsne_celseq2_p3
dev.off()

# CEL-seq2 P3 5cl PCA

cds_celseq2_p3_PCA <- reduce_dimension(cds_celseq2_p3, reduction_method = "PCA")

cds_celseq2_p3_PCA <- cluster_cells(cds_celseq2_p3_PCA, reduction_method = "PCA")

cds_celseq2_p3_PCA <- learn_graph(cds_celseq2_p3_PCA)

pca_celseq2_p3 <- plot_cells(cds_celseq2_p3_PCA, reduction_method = "PCA", 
                             color_cells_by = "cell_line_demuxlet", group_label_size = 3.5, 
                             show_trajectory_graph = FALSE, labels_per_group = 5)

png("monocle3-pca-celseq-p3.png", height = 1000, width = 1000)
pca_celseq2_p3
dev.off()
