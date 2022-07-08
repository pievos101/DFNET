## DFNET_utils - Various utilities.

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

#' Flatten data, a potentially 3-dimensional array, to a matrix to use
#' in ranger calls.
#' @param data the feature matrix or array
#' @param cols an optional column selection
#' @return the flattened data matrix
flatten2ranger <- function(data, cols) {
    if (length(dim(data)) == 2) {
        return(data[, cols])
    } else if (dim(data)[3] == 1) {
        return(as.matrix(data[, cols, 1]))
    } else {
        .data <- data[, cols, ]
        d <- dim(.data)
        # XXX: repeats colnames per mode, strips mode name
        dimnames <- list(
            dimnames(.data)[[1]],
            rep.int(dimnames(.data)[[2]], d[3])
        )
        return(matrix(.data, d[1], d[2] * d[3], dimnames = dimnames))
    }
}

#' Higher-order function to test predictors against data and target.
#'
#' @param data The data to use for testing
#' @param target The data to use for testing
#' @param metric The metric to use, a function taking (target, actual).
#' @param do.predict If TRUE (the default), run \code{predict} using the
#' predictor, otherwise access its \code{$predictions} directly.
#' @return a function which takes a predictor and evaluates \code{metric}
#' with \code{target} and the predictions of \code{predictor} on \code{data}.
tester <- function(data, target, metric = ModelMetrics::auc,
                   do.predict = TRUE) {
    if (do.predict) {
        function(predictor) {
            metric(target, predict(predictor, flatten2ranger(data))$predictions)
        }
    } else {
        function(predictor) {
            metric(target, predictor$predictions)
        }
    }
}

#' Extract the multi-modal target vector from \code{features}
#'
#' @param features potentially multi-modal features as a matrix
#' @return A list of shape \code{(target, is.multi_modal)}, where \code{target}
#' is the target vector and \code{is.multi_modal} is true if the features are
#' multi-modal.
multi_modal_target <- function(features) {
    is.multi_modal <- !is.data.frame(features)
    if (is.multi_modal) {
        target <- features[[1]][, "target"]
    } else {
        target <- features[, "target"]
    }
    return(list(target = target, is.multi_modal = is.multi_modal))
}
