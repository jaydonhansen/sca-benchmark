library(scran)
library(scater)
library(CellBench)
set_cellbench_threads(nthreads = 1)
setwd("~/advanced-bioinformatics/")
load("~/Downloads/advanced-bioinformatics-master/data/sincell_with_class.RData")
load("~/Downloads/advanced-bioinformatics-master/data/sincell_with_class_5cl.RData")

#make a list of dataset for cellbench
datasets <- list(sc_Celseq2_5cl_p1,
                 sc_Celseq2_5cl_p2,
                 sc_Celseq2_5cl_p3,
                 sce_sc_10x_5cl_qc,
)
#set the gene filter
gene_filter = function(sce){
  keep1 = (apply(counts(sce), 1, function(x) mean(x[x>0])) > 1.1)  # average count larger than 1.1
  keep2 = (rowSums(counts(sce)>0) > 10) # expressed in more than 10 cells
  sce = sce[(keep1 & keep2), ]
  return(sce)
}
tp = system.time({
logcounts(sce) = counts(sce)
})

#set the normalization methods
library(scran)

scran_norm = function(sce){
  tp = system.time({
    sce = computeSumFactors(sce)
    sce = normalize(sce) # goes to `logcounts` by default
  })
  
  method_name = "scran"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}


# for clustering
library(dplyr)
library(Seurat)
library(SingleCellExperiment)


# generate log counts 
sc_Celseq2_5cl_p1@assays$data$logcounts <- sc_Celseq2_5cl_p1@assays$data$counts
sc_Celseq2_5cl_p2@assays$data$logcounts <- sc_Celseq2_5cl_p2@assays$data$counts
sc_Celseq2_5cl_p3@assays$data$logcounts <- sc_Celseq2_5cl_p3@assays$data$counts
sce_sc_10x_5cl_qc@assays$data$logcounts <- sce_sc_10x_5cl_qc@assays@data@counts

# Pipe the CELLseq2 datasets through the regular PCA/t-SNE pipeline in Seurat
seurat_CELseq2_5cl_p1 <- sc_Celseq2_5cl_p1 %>% 
  as.Seurat(counts = "counts", data = "logcounts") %>%
  FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunTSNE %>% RunUMAP(dims = 1:5)

seurat_CELseq2_5cl_p2 <- sc_Celseq2_5cl_p2 %>% 
  as.Seurat(counts = "counts", data = "logcounts") %>%
  FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunTSNE %>% RunUMAP(dims = 1:5)

seurat_CELseq2_5cl_p3 <- sc_Celseq2_5cl_p3 %>% 
  as.Seurat(counts = "counts", data = "logcounts") %>%
  FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunTSNE %>% RunUMAP(dims = 1:5)

rm(sc_Celseq2_5cl_p1, sc_Celseq2_5cl_p2, sc_Celseq2_5cl_p3) # Save memory

# Pipe the 10x datasets through the regular PCA/t-SNE pipeline in Seurat
seurat_10x_5cl <- sce_sc_10x_5cl_qc %>% 
  as.Seurat(counts = "counts", data = "logcounts") %>%
  FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunTSNE %>% RunUMAP(dims = 1:5)

rm(sce_sc_10x_5cl_qc) # Save memory


# PCA/t-SNE/UMAP plots for CELseq2
pca_plot_CELseq2_5cl_p1 <- DimPlot(seurat_CELseq2_5cl_p1, reduction = "pca", group.by = "cell_line_demuxlet")
pca_plot_Celseq2_5cl_p1
pca_plot_CELseq2_5cl_p2 <- DimPlot(seurat_CELseq2_5cl_p2, reduction = "pca", group.by = "cell_line_demuxlet")
pca_plot_CELseq2_5cl_p2
pca_plot_CELseq2_5cl_p3 <- DimPlot(seurat_CELseq2_5cl_p3, reduction = "pca", group.by = "cell_line_demuxlet")
pca_plot_CELseq2_5cl_p3

tsne_plot_CELseq2_5cl_p1 <- DimPlot(seurat_CELseq2_5cl_p1, reduction = "tsne", group.by = "cell_line_demuxlet")
tsne_plot_Celseq2_5cl_p1
tsne_plot_CELseq2_5cl_p2 <- DimPlot(seurat_CELseq2_5cl_p2, reduction = "tsne", group.by = "cell_line_demuxlet")
tsne_plot_CELseq2_5cl_p2
tsne_plot_CELseq2_5cl_p3 <- DimPlot(seurat_CELseq2_5cl_p3, reduction = "tsne", group.by = "cell_line_demuxlet")
tsne_plot_CELseq2_5cl_p3


umap_plot_CELseq2_5cl_p1 <- DimPlot(seurat_CELseq2_5cl_p1, reduction = "umap", group.by = "cell_line_demuxlet")
umap_plot_Celseq2_5cl_p1
umap_plot_CELseq2_5cl_p2 <- DimPlot(seurat_CELseq2_5cl_p2, reduction = "umap", group.by = "cell_line_demuxlet")
umap_plot_CELseq2_5cl_p2
umap_plot_CELseq2_5cl_p3 <- DimPlot(seurat_CELseq2_5cl_p3, reduction = "umap", group.by = "cell_line_demuxlet")
umap_plot_CELseq2_5cl_p3


# PCA/t-SNE/UMAP plots for 10x
pca_plot_10x_5cl <- DimPlot(seurat_10x_5cl, reduction = "pca", group.by = "cell_line_demuxlet")
pca_plot_10x_5cl
tsne_plot_10x_5cl <- DimPlot(seurat_10x_5cl, reduction = "tsne", group.by = "cell_line_demuxlet")
tsne_plot_10x_5cl
umap_plot_10x_5cl <- DimPlot(seurat_10x_5cl, reduction = "umap", group.by = "cell_line_demuxlet")
umap_plot_10x_5cl

#save the plots as width=1000,height=1000

