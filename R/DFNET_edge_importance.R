## DFNET_edge_importance - Calculate the importance of edges.

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
