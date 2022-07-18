## methods.R - Operations for matured forests.

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

#' Return the first trees in the forest
#'
#' @description
#' Returns the first (= oldest) trees in the forest.  The number of trees to
#' keep is determined by multiplying the number of generations with the
#' generation size of the forest (the \code{ntrees} parameter)
#' in \link{train}.
#'
#' @details
#' This function should not be used to retrain the forest from an earlier
#' branch.  It does not return \code{last.performance} and may return a wrong
#' value for \code{walk.depth}.
#' @param x The original \code{DFNET.forest}.
#' @param n.generations integer. The number of generations to keep.
#' @param ... arguments to be passed to or from other methods.
#' @return The first \code{n.generations} generations worth of modules
#' and trees.
#' @method head DFNET.forest
#' @importFrom utils head
#' @export
head.DFNET.forest <- function(x, n.generations = 6L, ...) {
    n.trees <- attr(x, "generation_size")
    old.modules.weights <- head(x$modules.weights, n.generations * n.trees)

    structure(
        list(
            trees = head(x$trees, n.generations * n.trees),
            modules = head(x$modules, n.generations * n.trees),
            modules.weights = old.modules.weights
        ),
        class = "DFNET.forest",
        generation_size = n.trees,
        walk.depth = sapply(tail(old.modules.weights, n.trees), sum)
    )
}

#' Return the last trees in the forest
#'
#' @description
#' Returns the last (= newest) trees in the forest.  The number of trees to
#' keep is determined by multiplying the number of generations with the
#' generation size of the forest (the \code{ntrees} parameter)
#' in \link{train}.
#'
#' @details
#' This function can be used to shrink the forest while training (since only
#' the last generation will be used anyway).
#' @param x The original \code{DFNET.forest}.
#' @param n.generations integer. The number of generations to keep.
#' @param ... arguments to be passed to or from other methods.
#' @return The last \code{n.generations} generations worth of modules
#' and trees.
#' @method tail DFNET.forest
#' @importFrom utils tail
#' @export
tail.DFNET.forest <- function(x, n.generations = 6L, ...) {
    n.trees <- attr(x, "generation_size")
    return(structure(
        list(
            trees = tail(x$trees, n.generations * n.trees),
            modules = tail(x$modules, n.generations * n.trees),
            modules.weights = tail(x$modules.weights, n.generations * n.trees)
        ),
        class = "DFNET.forest",
        generation_size = n.trees,
        walk.depth = attr(x, "walk.depth"),
        last.performance = attr(x, "last.performance")
    ))
}

#' Model predictions
#'
#' Uses a \code{DFNET.forest} to run predictions on data.
#'
#' @param object The \code{DFNET.forest} to use for prediction.
#' @param data matrix or 3D array. The data to run predictions on.
#' @param ... arguments to be passed to or from other methods.
#' @return A vector of predicted classes.
#' @examples
#' \dontrun{
#' smp_size <- floor(0.80 * dim(features)[1])
#' train_ids <- sample(dim(features)[1], size = smp_size)
#' test_ids <- (1:dim(features)[1])[-train_ids]
#' forest <- DFNET_init(graph, features[train_ids, ])
#' # train ...
#' prediction <- predict(forest, features[test_ids, ])
#' }
#' @method predict DFNET.forest
#' @importFrom stats predict
#' @export
predict.DFNET.forest <- function(object, data, ...) {
    pred <- matrix(NaN, length(object$trees), dim(data)[1])
    data <- flatten2ranger(data)

    for (count in seq_along(object$trees)) {
        pred[count, ] <- predict(object$trees[[count]], data)$predictions
    }

    val <- apply(pred, 2, function(x) {
        return(as.numeric(names(sort(table(x), decreasing = TRUE))[1]))
    })
    names(val) <- dimnames(data)[[1]]
    return(val)
}
