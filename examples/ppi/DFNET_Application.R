# DFNET Application (data from linkedOmics)
# http://www.linkedomics.org/login.php#dataSource 

library(ranger)
library(igraph)
library(pROC)
library(DFNET)

## ppi
load(here::here("examples", "ppi", "ppi.rda"))

## features
load(here::here("examples", "ppi", "mRNA.rda"))
load(here::here("examples", "ppi", "Methy.rda"))

# outcome class
load(here::here("examples", "ppi", "target.rda"))

# train test split (80-20)
train_ids = sample(1:length(target), length(target)*0.80, replace=FALSE)
test_ids  = setdiff(1:length(target), train_ids) 
 
mRNA_train  = mRNA[train_ids,]
Methy_train = Methy[train_ids,]

mRNA_test  = mRNA[test_ids,]
Methy_test = Methy[test_ids,]
#-----------------------------

# get the induced graph structure
graph <- graph_from_edgelist(as.matrix(ppi), directed = FALSE)

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

# get the selected modules
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
test_features <- list(
    mRNA  = as.matrix(mRNA_test),
    Methy = as.matrix(Methy_test)
)
DATA_test <- DFNET:::flatten2ranger(common_features(test_features))

# Predict using the best DT
pred_best = predict(best_DT, DATA_test)$predictions

# predict using all detected modules
pred_all = predict(last_gen, DATA_test)$predictions

pred_best
pred_all

# Check the performance of the predictions
ModelMetrics::auc(pred_best, target[test_ids])
ModelMetrics::auc(pred_all, target[test_ids])

# TREE SHAP explanations
require(treeshap)
source(here::here("extensions", "treeshap.R"))

# unify dfnet forest to treeSHAP object
forest_unified = dfnet.unify(last_gen$trees, DATA_test)

# generate shapley values
forest_shap = treeshap(forest_unified, DATA_test)
