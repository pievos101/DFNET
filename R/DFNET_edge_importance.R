## DFNET_edge_importance - Calculate the importance of edges.

## Copyright © 2021, 2022 Bastian Pfeifer <bastianxpfeifer@gmail.com>

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

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

    tree_auc <- function(tree) area_under_curve(tree$predictions, mmt$target)

    tree_vars <- lapply(
        trees, function(x) {
            as.numeric(substring(names(x$variable.importance), 4))
        }
    )
    tree_aucs <- sapply(trees, tree_auc)

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
