library(dplyr)
library(monocle3)
library(SingleCellExperiment)
library(reticulate)
import("louvain")

setwd("~/advanced-bioinformatics/")
load("data/sincell_with_class.RData")
load("data/sincell_with_class_5cl.RData")

