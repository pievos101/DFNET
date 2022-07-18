## DFNET_calc_feature_importance - Calculate the importance of features.

## Copyright © 2021 Bastian Pfeifer <bastianxpfeifer@gmail.com>
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

#' Importance of individual features
#'
#' @description
#' Calculates the average importance of features over multiple decision
#' trees.
#'
#' @param forest The trained \code{DFNET.forest}.
#' @param features matrix or 3D array. The features on which it was trained.
#' @param sep string.
#' @return A matrix of importance scores with one row per column in
#' \code{features} and one column per matrix (in a 3D array).
feature_importance <- function(forest, features, sep = "$") {
    orig.dim <- dim(features)
    importance.dimnames <- list(
        dimnames(features)[[2]],
        if (length(dim(features) == 2)) NULL else dimnames(features)[[3]]
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
        importance, nrow = orig.dim[[2]],
        dimnames = importance.dimnames
    ))
}
