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
sc_Celseq2_5cl_p3@assays$data$logcounts<-sc_Celseq2_5cl_p3@assays$data$counts
logcounts(sc_Celseq2_5cl_p3) = log2(logcounts(sc_Celseq2_5cl_p3)+1-min(min(logcounts(sc_Celseq2_5cl_p3)),0))
var.fit <- trendVar(sc_Celseq2_5cl_p3, method="loess", use.spikes=FALSE)
var.out <- decomposeVar(sc_Celseq2_5cl_p3, var.fit)
hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:topn], ]
rowData(sc_Celseq2_5cl_p3)$hi_var = FALSE
rowData(sc_Celseq2_5cl_p3)$hi_var[rownames(rowData(sc_Celseq2_5cl_p3)) %in% rownames(hvg.out)] = TRUE

var.genes = rownames(sc_Celseq2_5cl_p3)[rowData(sc_Celseq2_5cl_p3)$hi_var]
tmp = logcounts(sc_Celseq2_5cl_p3)[var.genes,]
sc <- SCseq(as.data.frame(as.matrix(tmp)))
sc <- filterdata(sc, mintotal=1, minexpr = 1, minnumber = 1,
                 LBatch = NULL, knn = 10, CGenes = NULL, FGenes = NULL, ccor = 0.4,
                 bmode = "RaceID")
sc@ndata = sc@expdata
sc@genes = rownames(sc@ndata)
sc@counts = rep(1,ncol(sc_Celseq2_5cl_p3))
names(sc@counts) = colnames(sc@ndata)
sc@cluster$features = sc@genes
sc <- compdist(sc, metric="pearson", FSelect = FALSE, knn = NULL)
sc <- clustexp(sc, sat = TRUE, samp = NULL, cln = NULL, clustnr = 30,
               bootnr = 50, rseed = 17000, FUNcluster = "kmedoids")
sc <- findoutliers(sc, probthr = 0.001, outminc = 5, outlg = 2,
                   outdistquant = 0.95)
colData(sc_Celseq2_5cl_p3)$clustering_res = as.factor(sc@cpart)
sc <- comptsne(sc, dimRed = FALSE,initial_cmd = TRUE, perplexity = 30,rseed = 15555)
plottsne(sc,final=FALSE)
