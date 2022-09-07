![DFNETLogo](https://github.com/pievos101/DFNET/blob/cran/DFNET_logo.png)

# Graph-guided random forest for gene set selection

## Paper 

@misc{https://doi.org/10.48550/arxiv.2108.11674,
  doi = {10.48550/ARXIV.2108.11674},
  
  url = {https://arxiv.org/abs/2108.11674},
  
  author = {Pfeifer, Bastian and Baniecki, Hubert and Saranti, Anna and Biecek, Przemyslaw and Holzinger, Andreas},
  
  keywords = {Artificial Intelligence (cs.AI), Machine Learning (cs.LG), Molecular Networks (q-bio.MN), FOS: Computer and information sciences, FOS: Computer and information sciences, FOS: Biological sciences, FOS: Biological sciences},
  
  title = {Graph-guided random forest for gene set selection},
  
  publisher = {arXiv},
  
  year = {2021},
  
  copyright = {Creative Commons Attribution Non Commercial Share Alike 4.0 International}
}

### https://arxiv.org/abs/2108.11674

## Installation
The DFNET R-package can be installed using devtools.

```{r}
install.packages("devtools")
library(devtools)
install_github("pievos101/DFNET")
```

## Usage

Lets load some data (PPI Network, Gene expression, and DNA Methylation)

```{r}
## ppi network
data(ppi)

## features
data(mRNA)
data(Methy)

# outcome class
data(target)
```

Splitting the dat into train and test set

```{r}
# train test split (80-20)
train_ids = sample(1:length(target), length(target)*0.80, replace=FALSE)
test_ids  = setdiff(1:length(target), train_ids) 
 
mRNA_train  = mRNA[train_ids,]
Methy_train = Methy[train_ids,]

mRNA_test  = mRNA[test_ids,]
Methy_test = Methy[test_ids,]
```

Now, we need to define the dfnet object.

```{r}
# get the induced graph structure
graph <- graph_from_edgelist(as.matrix(ppi), directed = FALSE)

features <- list(
    mRNA  = as.matrix(mRNA_train),
    Methy = as.matrix(Methy_train)
)

dfnet_graph <- launder(graph, features, threshold = NaN)
```

We are now ready to train the model using the Greedy Decision Forest.

```{r}
# train the GDF
dfnet_forest <- train(,
    dfnet_graph$graph,
    dfnet_graph$features, target[train_ids],
    importance="impurity",
    ntrees = 100, niter = 10,
    initial.walk.depth = 10
)
```

The selected modules (gene sets) are within the last ntrees=100.
We store these trees using the function tail().

```{r}
# get the selected modules
last_gen <- tail(dfnet_forest, 1)
trees    <- last_gen$trees
tree_imp <- attr(last_gen, "last.performance")
```

From these trees we compute edge, feature and module importance


```{r}
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
```

The best performing gene set (module) can be obtained as follows


```{r}
# get the best module
ids <- which.max(as.numeric(mod_imp[,"total"]))
best_DT <- last_gen$trees[[ids]]
```

Now, lets check the performance of that module on the independent test data set. We compare the results with the performance of all trees selected.

```{r}
# Prepare test data
colnames(mRNA_test)  = paste(colnames(mRNA_test),"$","mRNA", sep="")
colnames(Methy_test) = paste(colnames(Methy_test),"$","Methy", sep="")
DATA_test = as.data.frame(cbind(mRNA_test, Methy_test))

# Predict using the best DT
pred_best = predict(best_DT, DATA_test)$predictions

# predict using all detected modules
pred_all = predict(last_gen, DATA_test)$predictions

pred_best
pred_all

# Check the performance of the predictions
ModelMetrics::auc(pred_best, target[test_ids])
ModelMetrics::auc(pred_all, target[test_ids])
```

For computing tree-based SHAP values the 'treeshap' package is required (see https://github.com/ModelOriented/treeshap).

```{r}
# TREE SHAP explanationa
require(treeshap)

# unify dfnet forest to treeSHAP object
forest_unified = dfnet_unify(last_gen$trees, DATA_test)

# generate shapley values
forest_shap = treeshap(forest_unified, DATA_test)
```