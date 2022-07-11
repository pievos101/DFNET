## treeshap.R - DFNET extensions for compatibility with treeshap.

## Copyright © 2021 Bastian Pfeifer <bastianxpfeifer@gmail.com>

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

#' Unify decision trees into a single forest model.
#'
#' @param trees A list of trees
#' @param data The data set on which trees were trained
#' @return the unified forest model as a \code{data.frame}.
#' @examples \dontrun{
#'   forest <- DFNET_init(graph, features, target)
#'   forest <- DFNET_iterate(forest, graph, features, target)
#'   unified_model <- dfnet.unify(forest$trees, features)
#' }
dfnet.unify <- function(trees, data) {
    unified_model <- data.frame()
    i <- 0
    support <- 0
    for (tree in trees) {
        message(i + 1, " of ", length(trees), " trees unified")
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
    return(unified_forest)
}
