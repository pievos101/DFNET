## utils.R - Various utilities.

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

#' Flatten multi-modal data for ranger calls.
#'
#' Produces a flat matrix from \code{data}, a potentially 3-dimensional array.
#' The column names of the new matrix are produced by pasting the dimnames
#' together, using \code{sep} as separator.
#'
#' @param data matrix or 3D array. The data to flatten.
#' @param cols vector. The (optional) columns to use.
#' @param sep string. A separator to use when building the new column names.
#' @return A matrix containing \code{data} reduced to \code{cols}, with
#' the third dimension inserted as extra columns.
flatten2ranger <- function(data, cols, sep = "$") {
    if (length(dim(data)) == 2) {
        return(data[, cols])
    } else if (dim(data)[3] == 1) {
        return(as.matrix(data[, cols, 1]))
    } else {
        .data <- data[, cols, ]
        d <- dim(.data)
        dimname_pasta <- matrix("", nrow = d[2] * d[3], ncol=3)
        colnames(dimname_pasta) <- c("cols", "modes", "sep")

        dimname_pasta[, 1:2] <- as.matrix(expand.grid(
            dimnames(.data)[[2]],
            if (is.null(dimnames(.data)[[3]])) {
                paste("mode", 1:d[3], sep = "")
            } else {
                dimnames(.data)[[3]]
            }
        ))
        dimname_pasta[, 3] <- sep

        dimnames <- list(
            dimnames(.data)[[1]],
            sapply(
                seq(1, d[2] * d[3]),
                function (n) do.call(paste, as.list(dimname_pasta[n, ]))
            )
        )
        return(matrix(.data, d[1], d[2] * d[3], dimnames = dimnames))
    }
}

#' Higher-order model testing
#'
#' Returns a procedure that can be used as \code{performance} metric in
#' \link{train}.
#'
#' @param data matrix or 3D array. The data to use for testing.
#' @param target numeric vector. The target values to predict.
#' @param metric binary function. A metric to compare actual values
#' with predictions.
#' @param do.predict logical. Should \code{predict} be called (TRUE) or
#' \code{$predictions} be accessed directly (FALSE)?
#' @return A unary function which, takes a ranger-like predictor and evaluates
#' \code{metric} with \code{target} and the predictions of \code{predictor}
#' on \code{data}.
#' @importFrom stats predict
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
