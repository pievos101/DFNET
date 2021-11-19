# RFNET Application (linekd Omics)

library(ranger)
library(igraph)
library(pROC)

## PPI
#PPI      <- read.table("~/LinkedOmics/BRCA/BRCA_PPI.txt")
PPI      <- read.table("~/LinkedOmics/KIRC/KIDNEY_PPI.txt")

## Features
#mRNA     <- read.table("~/LinkedOmics/BRCA/BRCA_mRNA_FEATURES.txt")
mRNA     <- read.table("~/LinkedOmics/KIRC/KIDNEY_mRNA_FEATURES.txt")

#Methy    <- read.table("~/LinkedOmics/BRCA/BRCA_Methy_FEATURES.txt")
Methy    <- read.table("~/LinkedOmics/KIRC/KIDNEY_Methy_FEATURES.txt")

#Mut      <- read.table("~/LinkedOmics/BRCA/BRCA_Mut_FEATURES.txt")

# Outcome class
#TARGET   <- read.table("~/LinkedOmics/BRCA/BRCA_SURVIVAL.txt")
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

# Make data balanced -------------------------------------------- #
TT        <- table(unlist(TARGET))
id        <- which.min(TT)
down_samp <- TT[id]
class     <- as.numeric(names(TT[id]))
down_ids  <- sample(which(unlist(TARGET)!=class), down_samp)

ids1 <- down_ids
ids2 <- which(unlist(TARGET)==class)

TARGET2 <- unlist(TARGET)[c(ids1,ids2)]

DFNET_graph$Feature_Matrix[[1]] <- DFNET_graph$Feature_Matrix[[1]][c(ids1,ids2),]
DFNET_graph$Feature_Matrix[[2]] <- DFNET_graph$Feature_Matrix[[2]][c(ids1,ids2),]


# Create TRAIN set ----------------------------------- #
DFNET_graph_train <- DFNET_graph
## 80% of the sample size
smp_size  <- floor(0.80 * nrow(DFNET_graph$Feature_Matrix[[1]]))
train_ids <- sample(seq_len(nrow(DFNET_graph$Feature_Matrix[[1]])), size = smp_size)
DFNET_graph_train$Feature_Matrix[[1]] <- DFNET_graph$Feature_Matrix[[1]][train_ids,]
DFNET_graph_train$Feature_Matrix[[2]] <- DFNET_graph$Feature_Matrix[[2]][train_ids,]

table(TARGET2[train_ids])

# Create TEST set ------------------------------------ #
DFNET_graph_test  <- DFNET_graph
test_ids <- (1:nrow(DFNET_graph$Feature_Matrix[[1]]))[-train_ids]
DFNET_graph_test$Feature_Matrix[[1]] <- DFNET_graph$Feature_Matrix[[1]][test_ids,]
DFNET_graph_test$Feature_Matrix[[2]] <- DFNET_graph$Feature_Matrix[[2]][test_ids,]

table(TARGET2[test_ids])

DFNET_object <- DFNET(DFNET_graph_train, ntrees=100, niter=300, init.mtry=15)

DFNET_pred   <- DFNET_predict(DFNET_object, DFNET_graph_test)

DFNET_perf   <- DFNET_performance(DFNET_pred, TARGET2[test_ids])

DFNET_perf$byClass

# Sensitivity is more relevant because hit means death outcome
Sensitivity  <- DFNET_perf$byClass["Sensitivity"]
Specificity  <- DFNET_perf$byClass["Specificity"]

# Feature Selection - Importance Measures

DFNET_Eimp   <- DFNET_Edge_Importance(DFNET_graph_train, DFNET_object)

DFNET_mod    <- DFNET_modules(DFNET_graph_train, DFNET_object, DFNET_Eimp)

# Get the nodes of the top-1 module
Nodes        <- as.numeric(strsplit(DFNET_mod[1,1]," ")[[1]])

DFNET_Fimp   <- DFNET_calc_feature_importance(Nodes, DFNET_object, DFNET_graph_train)

# Test this module on the TEST set @TODO

# Generate some plots

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
