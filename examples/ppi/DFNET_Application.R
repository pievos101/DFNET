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

target <- read.table("~/LinkedOmics/KIRC-RANDOM/KIDNEY_RANDOM_TARGET.txt")
target <- as.numeric(target)

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

## new code

graph <- graph_from_edgelist(as.matrix(PPI[, 1:2]), directed = FALSE)

features <- list(
    mRNA = as.matrix(mRNA),
    Methy = as.matrix(Methy)
)

dfnet_graph <- launder(graph, features, threshold = NaN)

# train the GDF
dfnet_forest <- train(,
    dfnet_graph$graph,
    dfnet_graph$features, target,
    importance="impurity",
    ntrees = 100, niter = 10,
    initial.walk.depth = 10
)

# get the accuracy of the selected modules
last_gen <- tail(dfnet_forest, 1)
trees <- last_gen$trees
tree_imp <- attr(last_gen, "last.performance")

# edge importance
e_imp <- edge_importance(dfnet_graph$graph, trees, tree_imp)

# feature importance
feat_imp <- feature_importance(last_gen, dfnet_graph$features)

# module importance
mod_imp <- module_importance(
    dfnet_graph$graph,
    last_gen$modules,
    e_imp,
    tree_imp
)

# get the best module
ids <- which.max(as.numeric(mod_imp[,"total"]))
best_DT <- last_gen$trees[[ids]]

# Prediction @TODO
forest = concatenate(last_gen$trees)

d1 = dfnet_graph$features[,,1]
colnames(d1) = paste(colnames(d1),"$","mRNA", sep="")
d2 = dfnet_graph$features[,,2]
colnames(d2) = paste(colnames(d2),"$","Methy", sep="")
predict(forest, cbind(d1,d2))









#### OLD CODE PREDICTION & SHAP values

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
