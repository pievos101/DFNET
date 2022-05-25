## DFNET_utils - Various utilities.

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

#' Calculate the confusion matrix, accuracy, and other preformance metrics of
#' a DFNET classifier.
#'
#' @param prediction the predicted values
#' @param reference the actual values
#' @return the confusion matrix between \code{prediction} and \code{reference},
#' along with some metrics.
#' @examples
#' \dontrun{
#' classifier <- DFNET(training_graph)
#' prediction <- DFNET_predict(classifier, test_graph)
#' target <- ground_truth(test_graph)
#' perf <- DFNET_performance(prediction, target)
#' }
DFNET_performance <- function(prediction, reference) {
    require(caret)
    require(e1071)

    res <- confusionMatrix(as.factor(pred),
        as.factor(reference),
        mode = "prec_recall",
        positive = "1"
    )

    return(res)
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

#' Calculate the area under curve w.r.t. \code{target} for the predictions made
#' by various decision trees.
#'
#' @param predictions the predictions
#' @param target the target vector optimized for
#' @return a vector of AUC values
area_under_curve <- function(predictions, target) {
    return(
        pROC::auc(
            target, predictions,
            na.rm = TRUE, levels = c(0, 1), direction = "<"
        )[1]
    )
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
