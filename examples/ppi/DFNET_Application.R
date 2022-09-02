# DFNET Application (linekd Omics)
library(ranger)
library(igraph)
library(pROC)
library(DFNET)

## PPI
PPI <- read.table("~/LinkedOmics/KIRC-RANDOM/KIDNEY_RANDOM_PPI.txt")

## Features
mRNA <- read.table("~/LinkedOmics/KIRC-RANDOM/KIDNEY_RANDOM_mRNA_FEATURES.txt")

Methy <- read.table("~/LinkedOmics/KIRC-RANDOM/KIDNEY_RANDOM_Methy_FEATURES.txt")

# outcome class
target <- read.table("~/LinkedOmics/KIRC-RANDOM/KIDNEY_RANDOM_TARGET.txt")
target <- as.numeric(target)

# Replace NANs with mean
na.ids <- which(apply(Methy, 2, function(x) {
    any(is.na(x))
}))

for (xx in na.ids) {
    ids <- which(is.na(Methy[, xx]))
    Methy[ids, xx] <- mean(Methy[, xx], na.rm = TRUE)
}

# reduce data dimension for test purposes
mRNA  = mRNA[,1:100]
Methy = Methy[,1:100]

#-----------------------------

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
trees    <- last_gen$trees
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

# convert2ranger
# forest = convert2ranger(last_gen$trees)

# Prepare test data
d1 = dfnet_graph$features[,,1]
colnames(d1) = paste(colnames(d1),"$","mRNA", sep="")
d2 = dfnet_graph$features[,,2]
colnames(d2) = paste(colnames(d2),"$","Methy", sep="")
DATA = as.data.frame(cbind(d1,d2))

# Predict using the best DT
pred_best = predict(best_DT, DATA)$predictions

# predict using all detected modules
pred_all = predict(last_gen, DATA)$predictions

pred_best
pred_all

# Check the performance of the predictions
ModelMetrics::auc(pred_best, target)
ModelMetrics::auc(pred_all, target)

# TREE SHAP explanationa
require(treeshap)
forest_shap <- DFNET_explain(DFNET_object, DFNET_graph_test)
sv <- forest_shap$shaps



