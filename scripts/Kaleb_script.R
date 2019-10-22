library(ggplot2)
library(data.table)
library(clusteval)
library(Seurat)

iris_means <- kmeans(iris[,-5], centers = 3)$cluster
iris_means

iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
iris_hclust

cluster_similarity(iris_means, iris_hclust)
########################################################################3
library(scran)
library(cluster)
library(scater)
library(mclust)

## Method for determining the Rand index of each clustering method 

ARI_matricc= function(sce){
  if(!("clustering_res" %in% colnames(colData(sce)))){
    return(NA)
  }
  if ("group" %in% colnames(colData(sce))){
    ari_val = adjustedRandIndex(sce$group, sce$clustering_res)
    print("here")
  }else{
    ari_val1 = adjustedRandIndex(sce$cell_line_demuxlet, 
                                 sce$clust_mono)
    ari_val2 = adjustedRandIndex(sce$cell_line_demuxlet, 
                                 sce$clust_race)
    ari_val3 = adjustedRandIndex(sce$cell_line_demuxlet, 
                                 sce$clust_seur)
    
    ari_val_col = c("Mono", ari_val1, 
                    
                    "Race", ari_val2, "Seur", ari_val3)
    
      }
  
  return(ari_val_col)
}

get_cluster_purity=function(sce, t){
  if(!("clustering_res" %in% colnames(colData(sce)))){
    return(NA)
  }
  if ("group" %in% colnames(colData(sce))){
    return(mean(unlist(lapply(unique(sce$group),function(x){cal_entropy(sce$clustering_res[sce$group==x])}))))
  }else{
    seurat_pur <- mean(unlist(lapply(unique(sce$cell_line_demuxlet),function(x){cal_entropy(sce$clust_seur[sce$cell_line_demuxlet==x])})))
    monocle_pur <- mean(unlist(lapply(unique(sce$cell_line_demuxlet),function(x){cal_entropy(sce$clust_mono[sce$cell_line_demuxlet==x])})))
    raceid_pur <-mean(unlist(lapply(unique(sce$cell_line_demuxlet),function(x){cal_entropy(sce$clust_race[sce$cell_line_demuxlet==x])})))
  }
  
  x <- data.frame("Method" = c("Seurat", "Monocle3", "RaceID"), 
                  "Purity" = c(seurat_pur, monocle_pur,raceid_pur),
                  "Cell Sample" = t)
  
  return (x)
}

get_cluster_accuracy=function(sce, t){
  if(!("clustering_res" %in% colnames(colData(sce)))){
    return(NA)
  }
  if ("group" %in% colnames(colData(sce))){
    return(mean(unlist(lapply(unique(sce$clustering_res),function(x){cal_entropy(sce$group[sce$clustering_res==x])}))))
  }else{
    seurat_acc <- mean(unlist(lapply(unique(sce$clust_seur),function(x){cal_entropy(sce$cell_line_demuxlet[sce$clust_seur==x])})))
    monocle_acc <- mean(unlist(lapply(unique(sce$clust_mono),function(x){cal_entropy(sce$cell_line_demuxlet[sce$clust_mono==x])})))
    raceid_acc <- mean(unlist(lapply(unique(sce$clust_race),function(x){cal_entropy(sce$cell_line_demuxlet[sce$clust_race==x])})))
    
  }
  
  x <- data.frame("Method" = c("Seurat", "Monocle3", "RaceID"), 
                  "Accuracy" = c(seurat_acc, monocle_acc,raceid_acc),
                  "Cell Sample" = t)
  
  return (x)
  
}

cluster_number = function(sce){
  if(!("clustering_res" %in% colnames(colData(sce)))){
    return(NA)
  }
  return(length(table(sce$clustering_res)))
}

cal_entropy=function(x){
  freqs <- table(x)/length(x)
  freqs = freqs[freqs>0]
  return(-sum(freqs * log(freqs)))
}




#####################################################################
# ____DROPSeq
#####################################################################
# Load the clustered variables into new vars 
#   > Note that names of the clustered seqs may have changed. 
#       Make sure to look it up
RaceID_DropSeq = sc_CEL
seruat_dropseq = seurat_Dropseq_3cl
mono_dropseq = cds_dropseq

# Save the clustering results to the original dataset 

sce_sc_Dropseq_qc$clust_mono <- as.factor(mono_dropseq@clusters$UMAP$clusters)
sce_sc_Dropseq_qc$clust_seur <- as.factor(seruat_dropseq$clustering_res)
sce_sc_Dropseq_qc$clust_race <- as.factor(RaceID_DropSeq@cpart)

# run the rand index formula 
ri_dropseq <- ARI_matricc(sce_sc_Dropseq_qc)
purity_drop <- get_cluster_purity(sce_sc_Dropseq_qc, "dropseq")
accuracy_drop <- get_cluster_accuracy(sce_sc_Dropseq_qc, "dropseq")
total_drop <- data.frame("Method" = c("Seruat", "Monocle3", "RaceID"), 
                "Purity" = c(purity_drop$Purity),
                "Accuracy" = c(accuracy_drop$Accuracy))
#####################################################################
# ____CELSeq
#####################################################################
# Load the clustered variables into new vars 
#   > Note that names of the clustered seqs may have changed. 
#       Make sure to look it up
RaceID_CELseq = sc_CEL
seurat_CELseq = seurat_CELseq2_3cl
mono_CELseq = cds_celseq2

# Save the clustering results to the original dataset 

sce_sc_CELseq2_qc$clust_mono <- as.factor(mono_CELseq@clusters$UMAP$clusters)
sce_sc_CELseq2_qc$clust_seur <- as.factor(seurat_CELseq$clustering_res)
sce_sc_CELseq2_qc$clust_race <- as.factor(RaceID_CELseq@cpart)

# run the rand index formula 

ri_celseq <- ARI_matricc(sce_sc_CELseq2_qc)
purity_celseq <- get_cluster_purity(sce_sc_CELseq2_qc, "celseq")
accuracy_celseq <- get_cluster_accuracy(sce_sc_CELseq2_qc, "celseq")
total_celseq <- data.frame("Method" = c("Seruat", "Monocle3", "RaceID"), 
                "Purity" = c(purity_celseq$Purity),
                "Accuracy" = c(accuracy_celseq$Accuracy))
#####################################################################
# ____5cl_Seq
#####################################################################
# Load the clustered variables into new vars 
#   > Note that names of the clustered seqs may have changed. 
#       Make sure to look it up
RaceID_5cl_p3 = sc_5_cell
seurat_5cl_p3 = seurat_CELseq2_5cl_p3
mono_5cl_p3 = cds_celseq2_p3
View(mono_5cl_p3@clusters$UMAP$clusters)

# Save the clustering results to the original dataset 

sc_Celseq2_5cl_p3$clust_mono <- as.factor(mono_5cl_p3@clusters$UMAP$clusters)
sc_Celseq2_5cl_p3$clust_seur <- as.factor(seurat_5cl_p3$clustering_res)
sc_Celseq2_5cl_p3$clust_race <- as.factor(RaceID_5cl_p3@cpart)

# Run the rand index formula

ri_5cl <- ARI_matricc(sc_Celseq2_5cl_p3)
purity_5cl <- get_cluster_purity(sc_Celseq2_5cl_p3, "5cl")
accuracy_5cl <- get_cluster_accuracy(sc_Celseq2_5cl_p3, "5cl")
total_5cl <- data.frame("Method" = c("Seruat", "Monocle3", "RaceID"), 
                "Purity" = c(purity_5cl$Purity),
                "Accuracy" = c(accuracy_5cl$Accuracy))
#####################################################################
# ____10x_Seq
#####################################################################
RaceID_10x_qc = sc
seurat_10_qc = seurat_10x_3cl
mono_10x = cds_10x

# Save the clustering results to the original dataset 

sce_sc_10x_qc$clust_race <- as.factor(RaceID_10x_qc@cpart)
sce_sc_10x_qc$clust_seur <- as.factor(seurat_10_qc$clustering_res)
sce_sc_10x_qc$clust_mono <- as.factor(mono_10x@clusters$UMAP$clusters)

# Run the rand index formula
ri_10x <- ARI_matricc(sce_sc_10x_qc)
purity_10x <- get_cluster_purity(sce_sc_10x_qc, "10x")
accuracy_10x <- get_cluster_accuracy(sce_sc_10x_qc, "10x")
total_10x <- data.frame("Method" = c("Seruat", "Monocle3", "RaceID"), 
                "Purity" = c(purity_10x$Purity),
                "Accuracy" = c(accuracy_10x$Accuracy))
#####################################################################
# Totals by method
#####################################################################
total_seur <- data.frame("Sample" = c("sc_10x_5cl", "sc_Drop-seq", 
                                      "sc_10x", "sc_CEL-seq2"),
                         
                         "Purity" =  c(0.22536546, 0.8267437, 
                                       1.3112102, 0.40728487),
                         
                         "Accuracy" = c(0.1259798, 0.09952274,
                                        0.005133893, 0))
  
total_mono <-data.frame("Sample" = c("sc_10x_5cl", "sc_Drop-seq", 
                                     "sc_10x", "sc_CEL-seq2"),
                        
                        "Purity" =  c(0.08502911, 0,
                                      0.227893, 0.052852),
                        
                        "Accuracy" = c(0.107844, 0.34356025, 
                                       0, 0.3744704))


total_race <- data.frame("Sample" = c("sc_10x_5cl", "sc_Drop-seq", 
                                      "sc_10x", "sc_CEL-seq2"),
                         
                         "Purity" =  c(0.63936587, 0.8267437,
                                       1.3112102, 0.53556053),
                         
                         "Accuracy" = c(0.1118639, 0.09952274,
                                        0.005133893, 0))




#####################################################################
# Graphs and etc 
#####################################################################

plt10x <- ggplot(total_10x, aes(x = Purity,
              y = Accuracy,
              colour = Method, 
              label = Method)) + geom_point(size = 5, aes(shape = Method)) + 
  ggtitle("Accuracy and Purity of Clusters in sc_10x Sample")+ 
  labs(x = "Entropy of Cluster Purity",
       y = "Entropy of Cluster Accuracy")

plt5cl <- ggplot(total_5cl, aes(x = Purity,
                        y = Accuracy,
                        colour = Method)) + geom_point(size = 5, aes(shape = Method)) + 
  ggtitle("Accuracy and Purity of Clusters in sc_10x_5cl Sample") + 
  labs(x = "Entropy of Cluster Purity",
       y = "Entropy of Cluster Accuracy")


pltdrop <- ggplot(total_drop, aes(x = Purity,
                         y = Accuracy,
                         colour = Method)) + geom_point(size = 5, aes(shape = Method)) + 
  ggtitle("Accuracy and Purity of Clusters in sc_Drop-seq Sample")+ 
  labs(x = "Entropy of Cluster Purity",
       y = "Entropy of Cluster Accuracy")


pltcelseq <- ggplot(total_celseq, aes(x = Purity,
                           y = Accuracy,
                           colour = Method)) + geom_point(size = 5, aes(shape = Method)) + 
  ggtitle("Accuracy and Purity of Clusters in sc_CEL-seq2 Sample")+ 
  labs(x = "Entropy of Cluster Purity",
       y = "Entropy of Cluster Accuracy")

plt10x
plt5cl
pltdrop
pltcelseq

pltseur <- ggplot(total_seur, aes(x = Purity, 
                                  y = Accuracy,
                                  colour = Sample)) + geom_point(size = 5, aes(shape = Sample)) + 
  ggtitle("Accuracy and Purity of Clusters Using Seurat ")+ 
  labs(x = "Entropy of Cluster Purity",
       y = "Entropy of Cluster Accuracy")

pltmono <- ggplot(total_mono, aes(x = Purity, 
                                  y = Accuracy,
                                  colour = Sample)) + geom_point(size = 5, aes(shape = Sample)) + 
  ggtitle("Accuracy and Purity of Clusters Using Monocle ")+ 
  labs(x = "Entropy of Cluster Purity",
       y = "Entropy of Cluster Accuracy")

pltrace <- ggplot(total_race, aes(x = Purity, 
                                  y = Accuracy,
                                  colour = Sample)) + geom_point(size = 5, aes(shape = Sample)) + 
  ggtitle("Accuracy and Purity of Clusters Using RaceID ")+ 
  labs(x = "Entropy of Cluster Purity",
       y = "Entropy of Cluster Accuracy")

pltseur
pltmono
pltrace

rangg <- ggplot(randindex, aes(x = Sample, 
                               y = `Rand Index`,
                               fill = Method)) + 
  geom_col(stat="identity", width=.5, position = "dodge") + 
  ggtitle('Rand Index for Method-Based Clustering of Each Sample Type')

rangg2 <- ggplot(randindex, aes(x = Method, 
                               y = `Rand Index`,
                               fill = Sample)) + 
  geom_col(stat="identity", width=.5, position = "dodge") + 
  ggtitle('Rand Index for Method-Based Clustering of Each Sample Type')
rangg
rangg2
