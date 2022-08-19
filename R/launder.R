## launder.R - Medieval methods for data cleaning.

## Copyright © 2021, 2022 Bastian Pfeifer <bastianxpfeifer@gmail.com>
## Copyright © 2022 Liliana Marie Prikler <liliana.prikler@ist.tugraz.at>

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

#' Common features in raw data
#'
#' Extracts common features from heterogenous data and simplifies them into a
#' three-dimensional array, thus making them homogenous.
#'
#' @param features list of matrices.
#' @export
#' @return A three-dimensional array containing the common features
#' (as per column names).
#' @family data laundering techniques
common_features <- function(features) {
    feature.names <- lapply(features, colnames)
    common.features.names <- Reduce(intersect, feature.names)

    common.features <- features
    for (xx in seq_along(common.features)) {
        common.features[[xx]] <- features[[xx]][, common.features.names]
    }

    return(simplify2array(common.features))
}

#' Retain meaningful features
#'
#' Trims \code{features} so that only columns whose names match the vertex names
#' in \code{graph} are retained.
#'
#' @param features matrix or 3D array.
#' @param graph an igraph.
#' @export
#' @return A matrix or 3D array (as features \code{features}), reduced and
#' reshaped such that only features belonging to the vertices in \code{graph}
#' appear in the order given by \code{graph}.
#' @family data laundering techniques
graphed_features <- function(features, graph) {
    if (length(dim(features)) == 2) {
        return(features[, igraph::vertex_attr(graph, "name")])
    } else {
        return(features[, igraph::vertex_attr(graph, "name"), ])
    }
}

#' Subgraph of a graph
#'
#' Creates a subgraph of a graph, containing only the specified vertices and
#' the edges between them.
#'
#' @param graph The original graph.
#' @param names vector of strings. The names of the vertices which will form
#' the subgraph.
#' @return The induced subraph.
#' @seealso \link[igraph:induced.subgraph]{igraph::induced_subgraph()}, which
#' this procedure uses internally.
#' @family data laundering techniques
#' @keywords internal
induced.subgraph.by_name <- function(graph, names) {
    igraph::induced.subgraph(
        graph,
        stats::na.omit(match(names, igraph::vertex_attr(graph, "name")))
    )
}

#' Remove lightweight edges
#'
#' Cuts off edges below a threshold.
#'
#' @param graph The original graph
#' @param threshold numeric. The minimal weight of an edge to keep.
#' @param threshold.quantile numeric. The minimal relative weight of an edge
#' to keep with respect to all other edge weights.
#' @importFrom igraph E
#' @export
#' @return A subgraph of \code{graph}, where edges have been removed according
#' to the given threshold. If no threshold is given \code{graph} is returned
#' unchanged.
#' @family data laundering techniques
cut_off <- function(graph, threshold = NaN, threshold.quantile = NaN) {
    if (is.na(threshold) && is.na(threshold.quantile)) {
        return(graph)
    }

    weights <- igraph::edge_attr(graph, "weights")

    if (is.na(threshold)) {
        threshold <- stats::quantile(weights, probs = threshold.quantile)
    }

    igraph::delete.edges(graph, E(graph)[which(weights < threshold)])
}

#' Data laundering
#'
#' Launders input graph and features.
#'
#' @param graph The original graph.
#' @param features list of matrices.
#' @param threshold see \link{cut_off}.
#' @param threshold.quantile see \link{cut_off}.
#' @export
#' @return A list \code{graph, features} with the laundered graph and features.
#' @family data laundering techniques
launder <- function(graph, features, threshold = NaN, threshold.quantile = NaN) {
    features <- common_features(features)
    graph <- induced.subgraph.by_name(graph, dimnames(features)[[2]])
    features <- graphed_features(features, graph)
    graph <- cut_off(
        graph,
        threshold = threshold, threshold.quantile = threshold.quantile
    )

    return(list(graph = graph, features = features))
}

#' Relativize data
#'
#' Maps the elements of a numeric vector to the interval [0:1].
#'
#' @param x numeric vector.
#' @param na.rm logical. Should NaNs be removed from \code{x} when calculating
#' minimal and maximal values?
#' @param default numeric. The default value(s) to use when all values in
#' \code{x} are equal. Should lie between 0 and 1.
#' @param default.na numeric. As default, but applied if minimal or maximal
#' value are NaN.  Not applicable if \code{na.rm} is true.
#' @return The relativized vector \eqn{\alpha}, such that
#' \eqn{x = (1-\alpha)\min(x) + \alpha\max(x)} holds w.r.t. floating point
#' inaccuracies.
#' @export
#' @examples
#' relat(1:3)
#' relat(c(NaN, 2, 3, 4))
#' relat(c(NaN, 2, 3, 4), na.rm = TRUE)
relat <- function(x, na.rm = FALSE, default = 1, default.na = NaN) {
    x0 <- min(x, na.rm = na.rm)
    x1 <- max(x, na.rm = na.rm)
    if (is.na(x0) || is.na(x1)) {
        return(rep_len(default.na, length(x)))
    } else if (x0 == x1) {
        return(rep_len(default, length(x)))
    } else {
        return((x - x0) / (x1 - x0))
    }
}
