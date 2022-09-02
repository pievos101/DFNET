# DFNET Application (data from linkedOmics)
library(ranger)
library(igraph)
library(pROC)
library(DFNET)

## ppi
PPI <- read.table("~/LinkedOmics/KIRC-RANDOM/KIDNEY_RANDOM_PPI.txt")

## features
mRNA <- read.table("~/LinkedOmics/KIRC-RANDOM/KIDNEY_RANDOM_mRNA_FEATURES.txt")

Methy <- read.table("~/LinkedOmics/KIRC-RANDOM/KIDNEY_RANDOM_Methy_FEATURES.txt")

# outcome class
target <- read.table("~/LinkedOmics/KIRC-RANDOM/KIDNEY_RANDOM_TARGET.txt")
target <- as.numeric(target)

# reduce data dimension for test purposes
mRNA  = mRNA[,1:100]
Methy = Methy[,1:100]

# train test split (80-20)
train_ids = sample(1:length(target), length(target)*0.80, replace=FALSE)
test_ids  = setdiff(1:length(target), train_ids) 
 
mRNA_train  = mRNA[train_ids,]
Methy_train = Methy[train_ids,]

mRNA_test  = mRNA[test_ids,]
Methy_test = Methy[test_ids,]
#-----------------------------

graph <- graph_from_edgelist(as.matrix(PPI[, 1:2]), directed = FALSE)

features <- list(
    mRNA  = as.matrix(mRNA_train),
    Methy = as.matrix(Methy_train)
)

dfnet_graph <- launder(graph, features, threshold = NaN)

# train the GDF
dfnet_forest <- train(,
    dfnet_graph$graph,
    dfnet_graph$features, target[train_ids],
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
colnames(mRNA_test)  = paste(colnames(mRNA_test),"$","mRNA", sep="")
colnames(Methy_test) = paste(colnames(Methy_test),"$","Methy", sep="")
DATA = as.data.frame(cbind(mRNA_test, Methy_test))

# Predict using the best DT
pred_best = predict(best_DT, DATA)$predictions

# predict using all detected modules
pred_all = predict(last_gen, DATA)$predictions

pred_best
pred_all

# Check the performance of the predictions
ModelMetrics::auc(pred_best, target[test_ids])
ModelMetrics::auc(pred_all, target[test_ids])

# TREE SHAP explanationa
require(treeshap)
#forest_shap <- DFNET_explain(DFNET_object, DFNET_graph_test)
#ssv <- forest_shap$shaps



