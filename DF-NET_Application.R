# RFNET Application (linekd Omics)

source("~/GitLab/rf-feature-selection-on-graphs/DF-NET_main.R")
source("~/GitLab/rf-feature-selection-on-graphs/DF-NET_generate_Graph_MM.R")
source("~/GitLab/rf-feature-selection-on-graphs/DF-NET_generate_Graph_SM.R")
source("~/GitLab/rf-feature-selection-on-graphs/DF-NET_accuracy.R")
source("~/GitLab/rf-feature-selection-on-graphs/DF-NET_modules.R")
source("~/GitLab/rf-feature-selection-on-graphs/DFNET_generate_graph_Omics.R")
source("~/GitLab/rf-feature-selection-on-graphs/DF-NET_Edge_Importance.R")


library(ranger)
library(igraph)
library(pROC)


PPI      <- read.table("~/LinkedOmics/KIRC/KIDNEY_PPI.txt")
mRNA     <- read.table("~/LinkedOmics/KIRC/KIDNEY_mRNA_FEATURES.txt")
Methy    <- read.table("~/LinkedOmics/KIRC/KIDNEY_Methy_FEATURES.txt")
TARGET   <- read.table("~/LinkedOmics/KIRC/KIDNEY_SURVIVAL.txt")

#@FIXME -- UGLY
# Replace NANs with mean
na.ids <- which(apply(Methy,2,function(x){any(is.na(x))}))

for (xx in na.ids){

	ids <- which(is.na(Methy[,xx]))
	Methy[ids,xx] <- mean(Methy[,xx], na.rm=TRUE)

}
#-----------------------------

# Perform DFNET

DFNET_graph  <- DFNET_generate_graph_Omics(PPI, list(mRNA, Methy), TARGET, 0.95)

DFNET_object <- DFNET(DFNET_graph, ntrees=100, niter=10, init.mtry=15)

DFNET_acc    <- DFNET_accuracy(DFNET_graph, DFNET_object)

DFNET_Eimp   <- DFNET_Edge_Importance(DFNET_graph, DFNET_object)

DFNET_mod    <- DFNET_modules(DFNET_graph, DFNET_object, DFNET_Eimp)

