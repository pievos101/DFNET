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
