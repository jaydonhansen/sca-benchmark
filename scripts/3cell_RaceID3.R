library("scran")
library("DESeq2")
library("biomaRt")
library("tsne")
library("scater")
library("RaceID")
library("Seurat")
setGeneric("plottsne", function(object,final=TRUE) standardGeneric("plottsne"))

setMethod("plottsne",
          signature = "SCseq",
          definition = function(object,final){
            if ( length(object@tsne) == 0 ) stop("run comptsne before plottsne")
            if ( final & length(object@cpart) == 0 ) stop("run findoutliers before plottsne")
            if ( !final & length(object@cluster$kpart) == 0 ) stop("run clustexp before plottsne")
            part <- if ( final ) object@cpart else object@cluster$kpart
            plot(object@tsne,xlab="Dim 1",ylab="Dim 2",pch=20,cex=1.5,col="lightgrey")
            for ( i in 1:max(part) ){
              if ( sum(part == i) > 0 ) text(object@tsne[part == i,1],object@tsne[part == i,2],i,col=object@fcol[i],cex=.75,font=4)
            }
          }
)



topn=1000
#sce_sc_10x_qc
logcounts(sce_sc_10x_qc) = log2(logcounts(sce_sc_10x_qc)+
                                  1-min(min(logcounts(sce_sc_10x_qc)),0))
var.fit <- trendVar(sce_sc_10x_qc, method="loess", use.spikes=FALSE)
var.out <- decomposeVar(sce_sc_10x_qc, var.fit)
hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:topn], ]
rowData(sce_sc_10x_qc)$hi_var = FALSE
rowData(sce_sc_10x_qc)$hi_var[rownames(rowData(sce_sc_10x_qc)) %in% rownames(hvg.out)] = TRUE

var.genes = rownames(sce_sc_10x_qc)[rowData(sce_sc_10x_qc)$hi_var]
tmp = logcounts(sce_sc_10x_qc)[var.genes,]
sc <- SCseq(as.data.frame(as.matrix(tmp)))
sc <- filterdata(sc, mintotal=1, minexpr = 1, minnumber = 1,
                 LBatch = NULL, knn = 10, CGenes = NULL, FGenes = NULL, ccor = 0.4,
                 bmode = "RaceID")
sc@ndata = sc@expdata
sc@genes = rownames(sc@ndata)
sc@counts = rep(1,ncol(sce_sc_10x_qc))
names(sc@counts) = colnames(sc@ndata)
sc@cluster$features = sc@genes
sc <- compdist(sc, metric="pearson", FSelect = FALSE, knn = NULL)
sc <- clustexp(sc, sat = TRUE, samp = NULL, cln = NULL, clustnr = 30,
               bootnr = 50, rseed = 17000, FUNcluster = "kmedoids")
sc <- findoutliers(sc, probthr = 0.001, outminc = 5, outlg = 2,
                   outdistquant = 0.95)
colData(sce_sc_10x_qc)$clustering_res = as.factor(sc@cpart)
sc <- comptsne(sc, dimRed = FALSE,initial_cmd = TRUE, perplexity = 30,rseed = 15555)
plottsne(sc,final=FALSE)


#sce_sc_CELseq2_qc
logcounts(sce_sc_CELseq2_qc) = log2(logcounts(sce_sc_CELseq2_qc)+
                                  1-min(min(logcounts(sce_sc_CELseq2_qc)),0))
var.fit <- trendVar(sce_sc_CELseq2_qc, method="loess", use.spikes=FALSE)
var.out <- decomposeVar(sce_sc_CELseq2_qc, var.fit)
hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:topn], ]
rowData(sce_sc_CELseq2_qc)$hi_var = FALSE
rowData(sce_sc_CELseq2_qc)$hi_var[rownames(rowData(sce_sc_CELseq2_qc)) %in% rownames(hvg.out)] = TRUE

var.genes = rownames(sce_sc_CELseq2_qc)[rowData(sce_sc_CELseq2_qc)$hi_var]
tmp = logcounts(sce_sc_CELseq2_qc)[var.genes,]
sc_CEL <- SCseq(as.data.frame(as.matrix(tmp)))
sc_CEL <- filterdata(sc_CEL, mintotal=1, minexpr = 1, minnumber = 1,
                 LBatch = NULL, knn = 10, CGenes = NULL, FGenes = NULL, ccor = 0.4,
                 bmode = "RaceID")
sc_CEL@ndata = sc_CEL@expdata
sc_CEL@genes = rownames(sc_CEL@ndata)
sc_CEL@counts = rep(1,ncol(sce_sc_CELseq2_qc))
names(sc_CEL@counts) = colnames(sc_CEL@ndata)
sc_CEL@cluster$features = sc_CEL@genes
sc_CEL <- compdist(sc_CEL, metric="pearson", FSelect = FALSE, knn = NULL)
sc_CEL <- clustexp(sc_CEL, sat = TRUE, samp = NULL, cln = NULL, clustnr = 30,
               bootnr = 50, rseed = 17000, FUNcluster = "kmedoids")
sc_CEL <- findoutliers(sc_CEL, probthr = 0.001, outminc = 5, outlg = 2,
                   outdistquant = 0.95)
colData(sce_sc_CELseq2_qc)$clustering_res = as.factor(sc_CEL@cpart)
sc_CEL <- comptsne(sc_CEL, dimRed = FALSE,initial_cmd = TRUE, perplexity = 30,rseed = 15555)
plottsne(sc_CEL,final=FALSE)

#sce_sc_Dropseq_qc
logcounts(sce_sc_Dropseq_qc) = log2(logcounts(sce_sc_Dropseq_qc)+
                                      1-min(min(logcounts(sce_sc_Dropseq_qc)),0))
var.fit <- trendVar(sce_sc_Dropseq_qc, method="loess", use.spikes=FALSE)
var.out <- decomposeVar(sce_sc_Dropseq_qc, var.fit)
hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:topn], ]
rowData(sce_sc_Dropseq_qc)$hi_var = FALSE
rowData(sce_sc_Dropseq_qc)$hi_var[rownames(rowData(sce_sc_Dropseq_qc)) %in% rownames(hvg.out)] = TRUE

var.genes = rownames(sce_sc_Dropseq_qc)[rowData(sce_sc_Dropseq_qc)$hi_var]
tmp = logcounts(sce_sc_Dropseq_qc)[var.genes,]
sc_CEL <- SCseq(as.data.frame(as.matrix(tmp)))
sc_CEL <- filterdata(sc_CEL, mintotal=1, minexpr = 1, minnumber = 1,
                     LBatch = NULL, knn = 10, CGenes = NULL, FGenes = NULL, ccor = 0.4,
                     bmode = "RaceID")
sc_CEL@ndata = sc_CEL@expdata
sc_CEL@genes = rownames(sc_CEL@ndata)
sc_CEL@counts = rep(1,ncol(sce_sc_Dropseq_qc))
names(sc_CEL@counts) = colnames(sc_CEL@ndata)
sc_CEL@cluster$features = sc_CEL@genes
sc_CEL <- compdist(sc_CEL, metric="pearson", FSelect = FALSE, knn = NULL)
sc_CEL <- clustexp(sc_CEL, sat = TRUE, samp = NULL, cln = NULL, clustnr = 30,
                   bootnr = 50, rseed = 17000, FUNcluster = "kmedoids")
sc_CEL <- findoutliers(sc_CEL, probthr = 0.001, outminc = 5, outlg = 2,
                       outdistquant = 0.95)
colData(sce_sc_Dropseq_qc)$clustering_res = as.factor(sc_CEL@cpart)
sc_CEL <- comptsne(sc_CEL, dimRed = FALSE,initial_cmd = TRUE, perplexity = 30,rseed = 15555)
plottsne(sc_CEL,final=FALSE)




  