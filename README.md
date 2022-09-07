![DFNETLogo](https://github.com/pievos101/DFNET/blob/cran/DFNET_logo.png)

# Graph-guided random forest for feature subset selection

## Paper

### https://arxiv.org/abs/2108.11674

## Installation
The DFNET R-package can be installed using devtools.

```{r}
install.packages("devtools")
devtools::install_github("pievos101/DFNET")
```

## Usage

See our examples using [synthetic data sets](https://github.com/pievos101/DFNET/blob/cran/examples/barabasi/simulation.R)
or [real world cancer data](https://github.com/pievos101/DFNET/blob/cran/examples/ppi/DFNET_Application.R).

Generally speaking, DFNET follows a four step process:

1. Preparing the input data (graph and features)
2. Training the forest.
3. Finding useful decision trees.
4. Using these trees for evaluation.

### Preparing input data
DFNET expects an `igraph::igraph` and a 2D or 3D feature array, as well as a
target vector with the same number of rows as the array.
The vertex names of the graph should be the same as the column names of the array.
When in doubt, use `launder` or related functions to prepare the input data.

### Training the forest
Once you have your graph and features, you can train your forest like so:

```{r}
forest <- train(,
    graph, features, target,
    ...
)
```

If you have a pre-trained forest, you can use that for training as well:

```{r}
forest <- train(forest,
    graph, features, target,
    ...
)
```



### Finding useful trees
Since DFNET performs greedy optimization, the last generation of trees
is the best according to the provided test metric.  DFNET provides overrides for
the standard R methods `head` and `tail`, which return generation.

```{r}
# get the selected modules
last_gen <- tail(forest, 1)
tree_imp <- attr(last_gen, "last.performance")
```

Note, that performance metrics for earlier generations are not kept.
Several importance scores can be derived from these metrics.

```{r}
e_imp <- edge_importance(graph, last_gen$trees, tree_imp)

f_imp <- feature_importance(last_gen, features)

m_imp <- module_importance(
    graph,
    last_gen$modules,
    e_imp,
    tree_imp
)
```

The module importance is particularly useful for feature selection, as it
combines the importance of edges within a module with the overall accuracy
of the decision tree.  You can use it to order decision trees or simply
extract the best one.

```{r}
best <- which.max(as.numeric(m_imp[, "total"]))
best.tree <- last_gen$trees[[best]]
```

```{r}
by_importance <- order(m_imp[, "total"], decreasing = TRUE)
last_gen$trees[by_importance]
```

### Using these trees for evaluation

DFNET provides an override for the `predict` method, that functions much like ranger's.

```{r}
# Predict using the best DT
pred_best = predict(best.tree, test_data)$predictions

# predict using all detected modules
pred_all = predict(last_gen, test_data)$predictions
```

You can use [ModelMetrics](https://cran.r-project.org/web/packages/ModelMetrics/index.html) to
evaluate the accuracy, precision, recall, or other performance metrics.

```
ModelMetrics::auc(pred_best, test_target)
ModelMetrics::auc(pred_all, test_target)
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

Finally, we provide [an extension](https://github.com/pievos101/DFNET/blob/cran/extensions/treeshap.R) to
compute tree-based SHAP values via [treeshap](https://github.com/ModelOriented/treeshap).

```{r}
forest_unified = dfnet.unify(last_gen$trees, test_data)
forest_shap = treeshap(forest_unified, test_data)
```

## BibTeX Citation

```
@misc{https://doi.org/10.48550/arxiv.2108.11674,
    doi = {10.48550/ARXIV.2108.11674},
    url = {https://arxiv.org/abs/2108.11674},
    author = {Pfeifer, Bastian and Baniecki, Hubert and Saranti, Anna and Biecek, Przemyslaw and Holzinger, Andreas},
    title = {Graph-guided random forest for gene set selection},
    publisher = {arXiv},
    year = {2021},
}
```
