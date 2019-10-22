library(scran)
library(cluster)
library(scater)
library(mclust)

ARI_matric = function(sce){
    if(!("clustering_res" %in% colnames(colData(sce)))){
        return(NA)
    }
    if ("group" %in% colnames(colData(sce))){
        ari_val = adjustedRandIndex(sce$group, sce$clustering_res)
    }else{
        ari_val = adjustedRandIndex(sce$cell_line_demuxlet, sce$clustering_res)
    }
    
    return(ari_val)
}

ARI_matric(sce_sc_10x_qc)
