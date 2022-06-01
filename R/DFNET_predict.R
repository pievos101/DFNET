## DFNET_predict - Use greedy random forest as ML classifier.

## Copyright © 2021, 2022 Bastian Pfeifer <bastianxpfeifer@gmail.com>

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

#' Use DFNET as a classifier on graphs it wasn't trained on.
#'
#' @param classifier the classifier
#' @param graph the graph on which @code{classifier} is tested
#' @param n.last.trees if given, use only the N last trees for classification
#' @param tree.ID if given, use the tree given by this ID for classification
#' @examples
#' \dontrun{
#' training_graph <- graph
#' test_graph <- graph
#' smp_size <- floor(0.80 * nrow(graph$feature.matrix[[1]]))
#' train_ids <- sample(seq_len(nrow(graph$feature.matrix[[1]])),
#'     size = smp_size
#' )
#' test_ids <- (1:nrow(graph$feature.matrix[[1]]))[-train_ids]
#' for (xx in 1:length(training_graph$feature.matrix)) {
#'     training_graph$feature.matrix[[xx]] <-
#'         graph$feature.matrix[[xx]][train_ids, ]
#' }
#' for (xx in 1:length(testing_graph$feature.matrix)) {
#'     training_graph$feature.matrix[[xx]] <-
#'         graph$feature.matrix[[xx]][test_ids, ]
#' }
#' classifier <- DFNET(training_graph)
#' prediction <- DFNET_predict(classifier, test_graph)
#' }
DFNET_predict <- function(classifier, graph, n.last.trees = NaN, tree.ID = NaN) {
    # retrieve full data
    dataset <- do.call(cbind, graph[[2]])
    # leave only one target column
    n.mod <- length(graph[[2]])
    dataset <- dataset[, -head(which(colnames(dataset) == "target"), n.mod - 1)]

    pred <- data.frame()

    if (!is.na(tree.ID)) {
        DECISION_TREES_ALL <- classifier$DFNET_trees[tree.ID]
    } else if (!is.na(n.last.trees)) {
        DECISION_TREES_ALL <- tail(classifier$DFNET_trees, n.last.trees)
    } else {
        DECISION_TREES_ALL <- classifier$DFNET_trees
    }

    # Predict
    pred <- matrix(NaN, length(DECISION_TREES_ALL), dim(dataset)[1])
    count <- 1
    for (tree in DECISION_TREES_ALL) {
        # pred <- rbind(pred, predict(tree, dataset)$predictions)
        if (count %% 100 == 0) {
            message(count, " of ", length(DECISION_TREES_ALL), " trees")
        }
        pred[count, ] <- predict(tree, dataset)$predictions
        count <- count + 1
    }

    val <- apply(pred, 2, function(x) {
        return(as.numeric(names(sort(table(x), decreasing = TRUE))[1]))
    })

    names(val) <- rownames(dataset)
    return(val)
}
