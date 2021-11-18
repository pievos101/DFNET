# RFNET Application (linekd Omics)

library(ranger)
library(igraph)
library(pROC)

## PPI
PPI      <- read.table("~/LinkedOmics/BRCA/BRCA_PPI.txt")

## Features
mRNA     <- read.table("~/LinkedOmics/BRCA/BRCA_mRNA_FEATURES.txt")
Methy    <- read.table("~/LinkedOmics/BRCA/BRCA_Methy_FEATURES.txt")
Mut      <- read.table("~/LinkedOmics/BRCA/BRCA_Mut_FEATURES.txt")

# Outcome class
TARGET   <- read.table("~/LinkedOmics/BRCA/BRCA_SURVIVAL.txt")

#@FIXME -- UGLY
# Replace NANs with mean
na.ids <- which(apply(Methy,2,function(x){any(is.na(x))}))

for (xx in na.ids){

	ids <- which(is.na(Methy[,xx]))
	Methy[ids,xx] <- mean(Methy[,xx], na.rm=TRUE)

}
#-----------------------------

# Perform DFNET

DFNET_graph  <- DFNET_generate_graph_Omics(PPI, list(mRNA, Methy, Mut), TARGET, 0.95)

DFNET_object <- DFNET(DFNET_graph, ntrees=200, niter=300, init.mtry=20)

DFNET_acc    <- DFNET_accuracy(DFNET_graph, DFNET_object)

DFNET_Eimp   <- DFNET_Edge_Importance(DFNET_graph, DFNET_object)

DFNET_mod    <- DFNET_modules(DFNET_graph, DFNET_object, DFNET_Eimp)

# Get the nodes of the top-1 module
Nodes        <- as.numeric(strsplit(DFNET_mod[1,1]," ")[[1]])

DFNET_Fimp   <- DFNET_calc_feature_importance(Nodes, DFNET_object, DFNET_graph)


#Generate some plots

## GGPLOT
library(ggplot2)
library(reshape)

RES1 <- cbind(colnames(DFNET_Fimp),DFNET_Fimp[1,])
RES2 <- cbind(colnames(DFNET_Fimp),DFNET_Fimp[2,])
RES1 <- cbind(RES1,"mRNA")
RES2 <- cbind(RES2,"Methylation")

RES <- rbind(RES1,RES2)
rownames(RES) <- NULL
colnames(RES) <- c("Gene","IMP","Type")
RES     <- as.data.frame(RES)
RES$IMP <- as.numeric(RES$IMP)

p <- ggplot(RES, aes(fill=Type, y=IMP, x=Gene)) + 
    geom_bar(position="dodge", stat="identity") +
    ylab("Feature Importance") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot(p)    
