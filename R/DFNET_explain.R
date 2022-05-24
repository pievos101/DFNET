#' Explain the output of a DFNET classifiers via SHAP values.
#'
#' @param classifier the classifier
#' @param graph the graph on which @code{classifier} is tested
#' @param n.last.trees if given, use only the N last trees for calculating SHAP
#' @param tree.ID if given, calculate SHAP for the single tree given by this ID
#' @return the SHAP values for the decision forest in \code{classifier} with
#' restrictions applied as per \code{n.last.trees} and \code{tree.ID}.
#' \code{tree.ID} takes precedence over \code{n.last.trees}.  If neither is
#' given, the whole forest is used.
#' @examples
#' \dontrun{
#' classifier <- DFNET(training_graph)
#' prediction <- DFNET_predict(classifier, test_graph)
#' forest_shap <- DFNET_explain(classifier, test_graph)
#' # post-process ...
#' }
DFNET_explain <- function(classifier, graph, n.last.trees = NaN, tree.ID = NaN) {
    require(treeshap)
    # retrieve full data
    dataset <- do.call(cbind, graph[[2]])
    # leave only one target column
    n.mod <- length(graph[[2]])
    dataset <- dataset[, -head(which(colnames(dataset) == "target"), n.mod - 1)]

    if (!is.na(tree.ID)) {
        DECISION_TREES_ALL <- classifier$DFNET_trees[tree.ID]
    } else if (!is.na(n.last.trees)) {
        DECISION_TREES_ALL <- tail(classifier$DFNET_trees, n.last.trees)
    } else {
        DECISION_TREES_ALL <- classifier$DFNET_trees
    }

    # Unify
    forest <- unify_forest(DECISION_TREES_ALL, dataset)

    # calculate treeshap
    forest_shap <- treeshap(forest, dataset[, -dim(dataset)[2]])
    # sv <- forest_shap$shaps
    return(forest_shap)
}

#' Unify trees into a single forest model.
#'
#' @param trees A list of trees
#' @param data The data set on which trees were trained
#' @return the unified forest model as a \code{data.frame}.
unify_forest <- function(trees, data) {
    unified_model <- data.frame()
    i <- 0
    support <- 0
    for (tree in trees) {
        cat(i + 1, " of ", length(trees), " trees unified \n")
        unified_tree <- ranger.unify(tree, data)
        unified_tree$model$Tree <- i
        # XXX: Why do we add support equally to the categorizations?
        unified_tree$model$Yes <- unified_tree$model$Yes + support
        unified_tree$model$No <- unified_tree$model$No + support
        unified_tree$model$Missing <- unified_tree$model$Missing + support
        unified_model <- rbind(unified_model, unified_tree$model)
        i <- i + 1
        support <- support + nrow(unified_tree$model)
    }
    unified_model$Prediction <- unified_model$Prediction / length(trees)
    unified_forest <- unified_tree
    unified_forest$model <- unified_model
    unified_forest$data <- data
    unified_forest
}
