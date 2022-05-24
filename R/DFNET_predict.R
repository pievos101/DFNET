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
#' smp_size <- floor(0.80 * nrow(graph$Feature_Matrix[[1]]))
#' train_ids <- sample(seq_len(nrow(graph$Feature_Matrix[[1]])),
#'     size = smp_size
#' )
#' test_ids <- (1:nrow(graph$Feature_Matrix[[1]]))[-train_ids]
#' for (xx in 1:length(training_graph$Feature_Matrix)) {
#'     training_graph$Feature_Matrix[[xx]] <-
#'         graph$Feature_Matrix[[xx]][train_ids, ]
#' }
#' for (xx in 1:length(testing_graph$Feature_Matrix)) {
#'     training_graph$Feature_Matrix[[xx]] <-
#'         graph$Feature_Matrix[[xx]][test_ids, ]
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
            cat(count, " of ", length(DECISION_TREES_ALL), " trees \n")
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
