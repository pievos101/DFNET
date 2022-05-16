# DFNET Application (linekd Omics)
library(ranger)
library(igraph)
library(pROC)
library(DFNET)

## PPI
# PPI      <- read.table("~/LinkedOmics/BRCA/BRCA_PPI.txt")
# PPI      <- read.table("~/LinkedOmics/KIRC/KIDNEY_PPI.txt")
PPI <- read.table("~/LinkedOmics/KIRC-RANDOM/KIDNEY_RANDOM_PPI.txt")

## Features
# mRNA     <- read.table("~/LinkedOmics/BRCA/BRCA_mRNA_FEATURES.txt")
# mRNA     <- read.table("~/LinkedOmics/KIRC/KIDNEY_mRNA_FEATURES.txt")
mRNA <- read.table("~/LinkedOmics/KIRC-RANDOM/KIDNEY_RANDOM_mRNA_FEATURES.txt")

# Methy    <- read.table("~/LinkedOmics/BRCA/BRCA_Methy_FEATURES.txt")
# Methy    <- read.table("~/LinkedOmics/KIRC/KIDNEY_Methy_FEATURES.txt")
Methy <- read.table("~/LinkedOmics/KIRC-RANDOM/KIDNEY_RANDOM_Methy_FEATURES.txt")

# Mut      <- read.table("~/LinkedOmics/BRCA/BRCA_Mut_FEATURES.txt")

# Outcome class
# TARGET   <- read.table("~/LinkedOmics/BRCA/BRCA_SURVIVAL.txt")
# TARGET   <- read.table("~/LinkedOmics/KIRC/KIDNEY_SURVIVAL.txt")

TARGET <- read.table("~/LinkedOmics/KIRC-RANDOM/KIDNEY_RANDOM_TARGET.txt")
# @FIXME -- UGLY
# Replace NANs with mean
na.ids <- which(apply(Methy, 2, function(x) {
    any(is.na(x))
}))

for (xx in na.ids) {
    ids <- which(is.na(Methy[, xx]))
    Methy[ids, xx] <- mean(Methy[, xx], na.rm = TRUE)
}
#-----------------------------

# Read Data
DFNET_graph <- DFNET_generate_graph_Omics(PPI, list(mRNA, Methy), TARGET, 0.99)

# Make data balanced -------------------------------------------- #
TT <- table(unlist(TARGET))
id <- which.min(TT)
down_samp <- TT[id]
class <- as.numeric(names(TT[id]))
down_ids <- sample(which(unlist(TARGET) != class), down_samp)

ids1 <- down_ids
ids2 <- which(unlist(TARGET) == class)

TARGET2 <- unlist(TARGET)[c(ids1, ids2)]

for (xx in 1:length(DFNET_graph$feature.matrix)) {
    DFNET_graph$feature.matrix[[xx]] <- DFNET_graph$feature.matrix[[xx]][c(ids1, ids2), ]
}
# ---------------------------------------- #

# Pre-Processing
DFNET_graph <- DFNET_preprocess(DFNET_graph)


# Create TRAIN set ----------------------------------- #
DFNET_graph_train <- DFNET_graph
## 80% of the sample size
smp_size <- floor(0.80 * nrow(DFNET_graph$feature.matrix[[1]]))
train_ids <- sample(seq_len(nrow(DFNET_graph$feature.matrix[[1]])), size = smp_size)
for (xx in 1:length(DFNET_graph_train$feature.matrix)) {
    DFNET_graph_train$feature.matrix[[xx]] <- DFNET_graph$feature.matrix[[xx]][train_ids, ]
}
table(TARGET2[train_ids])

# Create TEST set ------------------------------------ #
DFNET_graph_test <- DFNET_graph
test_ids <- (1:nrow(DFNET_graph$feature.matrix[[1]]))[-train_ids]
for (xx in 1:length(DFNET_graph_test$feature.matrix)) {
    DFNET_graph_test$feature.matrix[[xx]] <- DFNET_graph$feature.matrix[[xx]][test_ids, ]
}
table(TARGET2[test_ids])

## DFNET GREEDY FOREST
# default values are:
# ntrees = 500
# niter = 20
# init.mtry = 20

# Perform DFNET
DFNET_object <- DFNET(DFNET_graph_train, ntrees = 50, niter = 10, init.mtry = 20)

print(length(DFNET_object$DFNET_trees))

# PREDICTION
DFNET_pred <- DFNET_predict(DFNET_object, DFNET_graph_test)

# PERFORMANCE
DFNET_perf <- DFNET_performance(DFNET_pred, as.factor(DFNET_graph_test$feature.matrix[[1]][, "target"]))

DFNET_perf

DFNET_perf$byClass

############################################
# Feature Selection - Importance Measures
############################################

DFNET_Eimp <- DFNET_Edge_Importance(DFNET_graph_train, DFNET_object)

DFNET_mod <- DFNET_modules(DFNET_graph_train, DFNET_object, DFNET_Eimp)

# Get the nodes of the top-1 module
Nodes <- as.numeric(strsplit(DFNET_mod[1, 1], " ")[[1]])

DFNET_Fimp <- DFNET_calc_feature_importance(Nodes, DFNET_object, DFNET_graph_train)

# Test this module on the TEST set @TODO
TREEID <- as.numeric(rownames(DFNET_mod)[1])

DFNET_pred_best <- DFNET_predict(DFNET_object, DFNET_graph_test, tree.ID = TREEID)
DFNET_perf_best <- DFNET_performance(DFNET_pred_best, as.factor(DFNET_graph_test$feature.matrix[[1]][, "target"]))

DFNET_perf_best

# TREE SHAP
require(treeshap)
forest_shap <- DFNET_explain(DFNET_object, DFNET_graph_test)
sv <- forest_shap$shaps


# Generate some plots

## GGPLOT
library(ggplot2)
library(reshape)

RES1 <- cbind(colnames(DFNET_Fimp), DFNET_Fimp[1, ])
RES2 <- cbind(colnames(DFNET_Fimp), DFNET_Fimp[2, ])
RES1 <- cbind(RES1, "mRNA")
RES2 <- cbind(RES2, "Methylation")

RES <- rbind(RES1, RES2)
rownames(RES) <- NULL
colnames(RES) <- c("Gene", "IMP", "Type")
RES <- as.data.frame(RES)
RES$IMP <- as.numeric(RES$IMP)

p <- ggplot(RES, aes(fill = Type, y = IMP, x = Gene)) +
    geom_bar(position = "dodge", stat = "identity") +
    ylab("Feature Importance") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

plot(p)
