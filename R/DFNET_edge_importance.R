#' Calculate the importance of particular edges in the graph.
#'
#' @param graph The graph to analyse
#' @param features The features associated with each node in the graph
#' @param trees The decision trees returned by DFNET iterations
#' @param mc.cores how many cores to use in parallel
#' @return the importance of each edge in \code{graph} w.r.t. \code{trees}.
edge_importance <- function(graph, features, trees, mc.cores = 1) {
    require(parallel) # parallel is a base package, it should always exist

    mmt <- multi_modal_target(features)
    edges <- as_edgelist(graph, names = TRUE)

    tree_vars <- lapply(
        trees, function(x) {
            as.numeric(substring(names(x$variable.importance), 4))
        }
    )
    tree_aucs <- auc_per_tree(trees, mmt$target)

    edge_imp <- numeric(dim(edges)[1])
    edges_list <- lapply(apply(edges, 1, list), unlist)

    for(xx in 1:length(tree_vars)) {
        if (xx %% 100 == 0) cat (xx, " of ", length(tree_vars), " trees\n")

        pred <- function(x) {all(is.element(x, tree_vars[[xx]]))}
        res <- unlist(mclapply(edges_list, pred, mc.cores = mc.cores))

        edge_imp[res] <- edge_imp[res] + tree_aucs[xx]
    }

    return(relat(edge_imp))
}

#' Calculate the importance of particular edges in the graph.
#'
#' @param DFNET_graph a list, whose first element is a graph and whose second
#' element is a matrix of features
#' @param DFNET_object an object as returned by \code{DFNET(DFNET_graph)}
#' @param parallel TRUE to compute importances in parallel
#' @return the importance of each edge in the graph of \code{DFNET_graph}
#' w.r.t. the decision trees in \code{DFNET_object}
DFNET_edge_importance <- function(DFNET_graph, DFNET_object, parallel = FALSE) {
    graph <- DFNET_graph[[1]]
    features <- DFNET_graph[[2]]
    trees <- DFNET_object$DFNET_trees
    if (parallel)
        mc.cores <- detectCores() - 2
    else
        mc.cores <- 1
    return(edge_importance(graph, features, trees, mc.cores = mc.cores))
}

DFNET_Edge_Importance <- DFNET_edge_importance
