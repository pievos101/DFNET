## metrics.R - Importance metrics.

## Copyright © 2021, 2022 Bastian Pfeifer <bastianxpfeifer@gmail.com>
## Copyright © 2022 Liliana Prikler <liliana.prikler@ist.tugraz.at>

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
#' @param trees The decision trees returned by DFNET iterations
#' @param tree_importances Numeric importance scores assigned to each tree
#' @param mc.cores how many cores to use in parallel
#' @param sep The seperator used to flatten the features from which
#' \code{trees} were generated.
#' @export
#' @family metrics
#' @return the importance of each edge in \code{graph} w.r.t. \code{trees}.
edge_importance <- function(graph, trees, tree_importances,
                            mc.cores = 1, sep = "$") {
    tree_vars <- lapply(
        trees, function(x) {
            sapply(
                strsplit(names(x$variable.importance), sep, fixed = TRUE),
                function(var) var[1]
            )
        }
    )

    edges <- igraph::as_edgelist(graph, names = TRUE)
    edge_imp <- numeric(dim(edges)[1])
    edges_list <- lapply(apply(edges, 1, list), unlist)

    for (xx in seq_along(tree_vars)) {
        if (xx %% 100 == 0) message(xx, " of ", length(tree_vars), " trees")

        pred <- function(x) {
            all(is.element(x, tree_vars[[xx]]))
        }
        res <- unlist(parallel::mclapply(edges_list, pred, mc.cores = mc.cores))

        edge_imp[res] <- edge_imp[res] + tree_importances[xx]
    }

    return(relat(edge_imp))
}

#' Importance of individual features
#'
#' @description
#' Calculates the average importance of features over multiple decision
#' trees.
#'
#' @param forest The trained \code{DFNET.forest}.
#' @param features matrix or 3D array. The features on which it was trained.
#' @param sep string.
#' @export
#' @family metrics
#' @return A matrix of importance scores with one row per column in
#' \code{features} and one column per matrix (in a 3D array).
feature_importance <- function(forest, features, sep = "$") {
    orig.dim <- dim(features)
    importance.dimnames <- list(
        dimnames(features)[[2]],
        if (length(dim(features)) == 2) NULL else dimnames(features)[[3]]
    )

    features <- flatten2ranger(features, sep = sep)
    feature.names <- dimnames(features)[[2]]
    trees <- forest$trees

    importance <- rep_len(0, length(feature.names))
    ## count <- rep_len(0, length(trees))
    for (tree in trees) {
        varimp <- tree$variable.importance
        varnam <- names(varimp)

        ids <- match(feature.names, varnam)
        sel <- !is.na(ids)
        ## count[sel] <- count[sel] + 1
        importance[sel] <- importance[sel] + varimp[ids[sel]]
    }
    importance <- importance / length(trees)

    return(matrix(
        importance,
        nrow = orig.dim[[2]],
        dimnames = importance.dimnames
    ))
}

#' Calculate the importance of particular modules in the graph.
#'
#' @param graph the graph used for training
#' @param modules the modules obtained by training
#' @param edge_importances the importance scores assigned to each edge in the
#' graph, e.g. by \code{edge_importance}.
#' @param tree_importances the importance scores assigned to each decision tree
#' in the forest
#' @param mc.cores how many cores to use in parallel
#' @export
#' @family metrics
#' @return the accumulated edge importance and total importance of each module.
module_importance <- function(graph, modules, edge_importances,
                              tree_importances, mc.cores = 1) {
    module_imp <- cbind(0, tree_importances)
    edges <- igraph::as_edgelist(graph, names = FALSE)

    for (e in 1:nrow(edges)) {
        edge <- edges[e, 1:2]
        pred <- function(module) {
            all(is.element(edge, module))
        }
        res <- unlist(parallel::mclapply(modules, pred, mc.cores = mc.cores))
        # XXX: What about edge weight?
        module_imp[res, 1] <- module_imp[res, 1] + edge_importances[e]
    }

    # FIXME: paper actually claims normalization by edge count, not module size
    module_imp[, 1] <- module_imp[, 1] / sapply(modules, length)
    module_imp[, 2] <- module_imp[, 1] + module_imp[, 2]
    colnames(module_imp) <- c("edge", "total")

    return(module_imp)
}

#' Calculate the importance of each unique module.
#'
#' @param graph the graph used for training
#' @param modules the modules obtained by training
#' @param edge_importances the importance scores assigned to each edge in the
#' graph, e.g. by \code{edge_importance}.
#' @param tree_importances the importance scores assigned to each decision tree
#' in the forest
#' @param mc.cores how many cores to use in parallel
#' @param collapse how accumulate multiple values per module
#' @export
#' @family metrics
#' @return the collapsed tree and edge importances of each unique module, as
#' well as their sum.
unique_module_importance <- function(graph, modules, edge_importances,
                                     tree_importances, mc.cores = 1,
                                     collapse = mean) {
    module_imp <- module_importance(
        graph, modules,
        edge_importances, tree_importances,
        mc.cores = mc.cores
    )

    modules <- lapply(modules, function(m) unique(sort(m)))
    module_name <- function(mod) paste(mod, collapse = " ")
    module_names <- sapply(modules, module_name)
    by_module_name <- order(module_names)

    module_imp <- module_imp[by_module_name, ]
    tree_importances <- tree_importances[by_module_name]
    modules <- modules[by_module_name]

    unique_modules <- unique(modules)
    unique_module_imp <- matrix(0, nrow = length(unique_modules), ncol = 3)
    colnames(unique_module_imp) <- c("tree", "edge", "total")
    rownames(unique_module_imp) <- sapply(unique_modules, module_name)

    for (module in unique_modules) {
        select <- sapply(modules, function(m) identical(m, module))

        unique_module_imp[module_name(module), 1:2] <- c(
            collapse(tree_importances[select]),
            collapse(module_imp[select, 1])
        )
    }
    unique_module_imp[, 3] <- unique_module_imp[, 1] + unique_module_imp[, 2]

    return(list(table = unique_module_imp, modules = unique_modules))
}
