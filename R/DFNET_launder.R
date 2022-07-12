## DFNET_launder.R - Medieval methods for data cleaning.

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

#' Extract common features from raw multi-omic data.
#'
#' @param features a list of feature matrices with equal rows but possibly
#' unequal columns
#' @return the common features (as per column names), simplified to an array.
common_features <- function(features) {
    feature.names <- lapply(features, colnames)
    common.features.names <- Reduce(intersect, feature.names)

    common.features <- features
    for (xx in 1:length(common.features)) {
        common.features[[xx]] <- features[[xx]][, common.features.names]
    }

    return(simplify2array(common.features))
}

#' Reduce features to the nodes of a graph.
#'
#' @param features an array of features
#' @param graph an igraph
#' @return \code{features} reduced to columns whose names appear in
#' \code{graph}
graphed_features <- function(features, graph) {
    if (length(dim(features)) == 2) {
        return(features[, igraph::vertex_attr(graph, "names")])
    } else {
        return(features[, igraph::vertex_attr(graph, "names"), ])
    }
}

#' Like \code{induced.subgraph} from the igraph package, but using names rather
#' than vertex IDs.
induced.subgraph.by_name <- function(graph, names) {
    igraph::induced.subgraph(
        graph,
        stats::na.omit(match(names, igraph::vertex_attr(graph, "names")))
    )
}

#' Cut off edges below a certain threshold.
#'
#' @param graph an igraph
#' @param threshold the minimal weight of an edge to be kept
#' @param threshold.quantile same as threshold, but derived as a quantile
#' of the edge weights
#' @return \code{graph} with edges whose weight is lower than \code{threshold}
#' removed.  If no threshold is given, return graph as-is.
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

#' Launder the input graph and features.
#'
#' @param graph the graph to launder
#' @param features the features to launder, as a list of matrices
#' @export
launder <- function(graph, features, threshold = NaN, threshold.quantile = NaN) {
    features <- common_features(features)
    graph <- induced.subgraph.by_name(graph, dimnames(features)[[2]])
    features <- graphed_features(features, graph)
    graph <- cut_off(graph)

    return(list(graph = graph, features = features))
}

#' Relativize a vector, so that all of its elements are between 0 and 1.
#'
#' @param x the vector to relativize
#' @param na.rm whether to remove NaNs from \code{x} when calculating
#' minimal and maximal values
#' @param default the default value(s) to use when all values in \code{x}
#' are equal.  Should lie between 0 and 1.
#' @param default.na as default, but applied if minimal or maximal value
#' are NaN.  Not applicable if \code{na.rm} is true.
#' @return the relativized vector alpha, such that
#' \eqn{x = (1-\alpha)\min(x) + \alpha\max(x)} holds w.r.t. floating point
#' inaccuracies.
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
