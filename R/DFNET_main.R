## DFNET_main - Network module detection via greedy random forest.

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

#' Decision tree learning from modules
#'
#' \code{learn_decisions} uses \code{ranger} to perform feature selection with
#' respect to \code{raw_modules}.
#'
#' @param raw_modules list of numeric vectors. The raw modules.
#' @param features numeric matrix or 3D array. The features to train on.
#' @param target numeric vector. The target to train towards.
#' @return A list of shape (\code{trees}, \code{modules},
#' \code{modules.weights}), where \code{modules} are the sorted
#' \code{raw_modules} with individual weights \code{modules.weights}, and
#' \code{trees} contains one ranger decision tree per module.
learn_decisions <- function(raw_modules, features, target) {
    modules_rle <- lapply(raw_modules, function(m) rle(sort(m)))

    decision_trees <- lapply(modules_rle, function(m) {
        unique_nodes <- m$values
        unique_nodes_weights <- m$lengths
        weights <- unique_nodes_weights

        mm_data <- flatten2ranger(features, unique_nodes)
        weights <- rep_len(weights, dim(mm_data)[2])

        # Perform feature selection
        ranger::ranger(
            x = mm_data,
            y = target,
            split.select.weights = weights / sum(weights),
            verbose = FALSE,
            classification = TRUE,
            importance = "impurity",
            num.trees = 1,
            mtry = dim(mm_data)[2] - 1,
            replace = TRUE
        )
    })

    return(
        list(
            trees = decision_trees,
            modules = lapply(modules_rle, function(mrle) mrle$values),
            modules.weights = lapply(modules_rle, function(mrle) mrle$lengths)
        )
    )
}

#' Model initialization
#'
#' Initializes the decision forest network.
#'
#' @param graph The graph to train the network on.
#' @param features matrix or 3D array. The features to train on.
#' @param target numeric vector. the target to train towards.
#' @param ntrees integer. The number of trees to generate per iteration.
#' @param walk.depth integer. The number of nodes to select per module.
#' @param performance unary function. Called with a decision tree as argument to
#' estimate that tree's performance.
#' @importFrom igraph V
#' @return An initialized \code{DFNET.forest}.
init <- function(graph, features, target,
                 ntrees = 100, walk.depth = NaN,
                 performance = NULL) {
    nodes <- V(graph)
    n.nodes <- length(nodes)

    if (is.na(walk.depth)) {
        walk.depth <- ceiling(sqrt(n.nodes))
    }
    walk.depth <- rep_len(walk.depth, ntrees)
    if (is.null(performance)) {
        performance <- function(tree) {
            ModelMetrics::auc(target, tree$predictions)
        }
    }

    count <- 1
    selected_nodes <- list()
    repeat {
        sampled.nodes <- sample(nodes, (ntrees + 1 - count), replace = TRUE)
        for (sn in sampled.nodes) {
            depth <- walk.depth[[count]]
            selected_nodes[[count]] <- as.numeric(igraph::random_walk(graph, sn, depth))
            # Pick only walks of maximal length
            if (length(selected_nodes[[count]]) >= walk.depth[[count]]) {
                count <- count + 1
            }
        }
        if (count > ntrees) {
            break
        }
    }

    seed <- learn_decisions(selected_nodes, features, target)
    last.perf <- sapply(seed$trees, performance)

    return(
        structure(
            seed,
            class = "DFNET.forest",
            generation_size = ntrees,
            walk.depth = walk.depth,
            last.performance = last.perf
        )
    )
}

#' Model training
#'
#' @description
#' Trains a decision forest on \code{feature} and \code{target}.
#'
#' @details
#' This function generates \code{ntrees} modules and decision trees per iteration
#' and greedily selects those which improve the \code{performance} metric.
#' The trees are trained on \code{features} and \code{target}.
#' \code{performance} can use its own validation set, or default to the
#' \code{features} and \code{target} above (the default), in which case ranger
#' handles the data split.
#'
#' In each iteration, this function tries to shrink modules which have
#' previously been improved.  \code{initial.walk.depth} thus gives the maximal
#' module size, whereas \code{min.walk.depth} specifies the smallest walk depth.
#'
#' Model training can be resumed from an already trained forest, in which case
#' the attributes of that forest are used in lieu of \code{ntrees} and
#' \code{initial.walk.depth}.  When resuming this training, it might make sense
#' to also specify the \code{offset} parameter for somewhat improved logging.
#'
#' The returned \code{DFNET.forest} is a list of shape (\code{trees},
#' \code{modules}, \code{modules.weights}), where \code{trees} are the decision
#' trees created for detected \code{modules}, and \code{modules.weights} gives
#' the weights used for each node.
#'
#' As "private" attributes used for iteration, \code{generation_size} is set to
#' \code{ntrees}, \code{walk.depth} captures the walk depth for the next
#' iteration, and \code{last.performance} to a vector of length \code{ntrees},
#' containing the result of \code{performance} of each tree w.r.t. \code{target}.
#'
#' @param forest a \code{DFNET.forest} or \code{null}.
#' @param graph The graph to train the forest on.
#' @param features matrix or 3D array. The features to train on.
#' @param target numeric vector. the target to train towards.
#' @param niter integer. The number of iterations to run.
#' @param offset integer. An offset added to the iteration count for logging
#' purposes.
#' @param min.walk.depth The integer minimal number of nodes to visit per tree
#' per iteration.
#' @param ntrees integer. The number of trees to generate per iteration.
#' @param initial.walk.depth integer. The number of nodes to visit per tree
#' during initialization.
#' @param performance unary function. Called with a decision tree as argument to
#' estimate that tree's performance.
#' @return The trained \code{DFNET.forest}.
#' @importFrom utils tail
#' @export
#' @examples
#' \dontrun{
#' forest <- NULL
#' offset <- 0
#' while (keep_iterating(forest, target)) { # insert your own iteration criteria
#'     forest <- train(
#'         forest,
#'         graph,
#'         features,
#'         niter = 10,
#'         offset = offset
#'         # ...
#'     )
#'     offset <- offset + 10
#' }
#' }
#'
train <- function(forest, graph, features, target,
                  niter = 200, offset = 0, min.walk.depth = 2,
                  ntrees = 100, initial.walk.depth = NaN,
                  performance = NULL) {
    stopifnot(niter >= 0, offset >= 0)
    if (missing(forest) || is.null(forest)) {
        forest <- init(
            graph, features, target,
            ntrees = ntrees, walk.depth = initial.walk.depth,
            performance = performance
        )
    }

    if (niter == 0) {
        return(forest)
    }

    if (is.null(performance)) {
        performance <- function(tree) {
            ModelMetrics::auc(target, tree$predictions)
        }
    }

    ntrees <- attr(forest, "generation_size")
    last.walk.depth <- attr(forest, "walk.depth")

    all.trees <- forest$trees
    all.modules <- forest$modules
    all.modules.weights <- forest$modules.weights
    last.trees <- tail(all.trees, ntrees)
    last.modules <- tail(all.modules, ntrees)
    last.modules.weights <- tail(all.modules.weights, ntrees)
    last.perf <- sapply(last.trees, performance)

    iter.min <- offset + 1
    iter.max <- offset + niter
    for (iter in iter.min:iter.max) {
        message(iter, " of ", iter.max, " greedy steps")

        ids_keep <- sample(1:length(last.perf), prob = last.perf, replace = TRUE)
        kept_modules <- last.modules[ids_keep]
        walk.depth <- last.walk.depth[ids_keep]

        start_nodes <- sapply(kept_modules, function(module) {
            sample(module, 1)
        })

        modules <- lapply(1:length(start_nodes), function(sn) {
            as.numeric(
                igraph::random_walk(graph, start_nodes[sn], walk.depth[sn])
            )
        })

        next_gen <- learn_decisions(modules, features, target)
        perf <- sapply(next_gen$trees, performance)

        good_enough <- perf >= last.perf

        # Update inner state
        last.perf <- ifelse(good_enough, perf, last.perf)
        last.modules <- ifelse(good_enough, next_gen$modules, last.modules)
        last.modules.weights <-
            ifelse(good_enough, next_gen$modules.weights, last.modules.weights)
        last.trees <- ifelse(good_enough, next_gen$trees, last.trees)
        last.walk.depth <- ifelse(good_enough, walk.depth - 1, last.walk.depth)
        last.walk.depth[last.walk.depth < min.walk.depth] <- min.walk.depth

        all.modules <- c(all.modules, last.modules)
        all.modules.weights <- c(all.modules.weights, last.modules.weights)
        all.trees <- c(all.trees, last.trees)
    }


    return(
        structure(
            list(
                trees = all.trees,
                modules = all.modules,
                modules.weights = all.modules.weights
            ),
            class = "DFNET.forest",
            generation_size = ntrees,
            walk.depth = last.walk.depth,
            last.performance = last.perf
        )
    )
}

#' Return the first trees in the forest
#'
#' @description
#' Returns the first (= oldest) trees in the forest.
#'
#' @details
#' This function should not be used to retrain the forest from an earlier
#' branch.  It does not return \code{last.performance} and may return a wrong
#' value for \code{walk.depth}.
#' @param forest The original \code{DFNET.forest}.
#' @param n.generations integer. The number of generations to keep.
#' @return The first \code{n.generations} generations worth of modules
#' and trees.
#' @method head DFNET.forest
#' @importFrom utils head
#' @export
head.DFNET.forest <- function(forest, n.generations = 6L) {
    n.trees <- attr(forest, "generation_size")
    old.modules.weights <- head(
        forest$modules.weights,
        n.generations * n.trees
    )

    structure(
        list(
            trees = head(forest$trees, n.generations * n.trees),
            modules = head(forest$modules, n.generations * n.trees),
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
#' Returns the last (= newest) trees in the forest.
#'
#' @details
#' This function can be used to shrink the forest while training (since only
#' the last generation will be used anyway).
#' @param forest The original \code{DFNET.forest}.
#' @param n.generations integer. The number of generations to keep.
#' @return The last \code{n.generations} generations worth of modules
#' and trees.
#' @method tail DFNET.forest
#' @importFrom utils tail
#' @export
tail.DFNET.forest <- function(forest, n.generations = 6L) {
    n.trees <- attr(forest, "generation_size")
    return(structure(
        list(
            trees = tail(forest$trees, n.generations * n.trees),
            modules = tail(forest$modules, n.generations * n.trees),
            modules.weights = tail(
                forest$modules.weights,
                n.generations * n.trees
            )
        ),
        class = "DFNET.forest",
        generation_size = n.trees,
        walk.depth = attr(forest, "walk.depth"),
        last.performance = attr(forest, "last.performance")
    ))
}

#' Model predictions
#'
#' Uses a \code{DFNET.forest} to run predictions on data.
#'
#' @param forest The \code{DFNET.forest} to use for prediction.
#' @param data matrix or 3D array. The data to run predictions on.
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
predict.DFNET.forest <- function(forest, data) {
    pred <- matrix(NaN, length(forest$trees), dim(data)[1])
    data <- flatten2ranger(data)

    for (count in 1:length(forest$trees)) {
        pred[count, ] <- predict(forest$trees[[count]], data)$predictions
    }

    val <- apply(pred, 2, function(x) {
        return(as.numeric(names(sort(table(x), decreasing = TRUE))[1]))
    })
    names(val) <- dimnames(data)[[1]]
    return(val)
}

#' Construct a new decision forest and run some iterations on it.
#'
#' @param DFNET_graph a list, whose first element is a graph and whose second
#' element is a matrix of features
#' @param ntrees how many trees should be generated initially and per iteration
#' @param niter how many iterations to run
#' @param init.mtry how many nodes to select per initial tree
#' @return a list of shape \code{(DFNET_graph, DFNET_trees, DFNET_MODULES,
#' DFNET_MODULES_AUC)}, where \code{DFNET_graph} is as in the input,
#' \code{DFNET_trees} are the generated trees, \code{DFNET_MODULES} are the
#' modules from which the trees were generated and \code{DFNET_MODULES_AUC} is
#' the area under curve of the trees.
#' @examples
#' \dontrun{
#' DFNET_graph <- DFNET_generate_graph_omics(graph, features, target)
#' DFNET_object <- DFNET(DFNET_graph)
#' # postprocess ...
#' }
DFNET <- function(DFNET_graph, ntrees = 100, niter = 200, init.mtry = NaN) {
    if (is.data.frame(DFNET_graph[[2]])) {
        features <- as.matrix(DFNET_graph[[2]])
        target <- features[, which(dimnames(features)[[2]] == "target")]
        features <- features[, -which(dimnames(features)[[2]] == "target")]
    } else {
        features <- simplify2array(lapply(DFNET_graph[[2]], as.matrix))
        target <- features[, which(dimnames(features)[[2]] == "target"), 1]
        features <- features[, -which(dimnames(features)[[2]] == "target"), ]
    }

    forest <- train(
        NULL, DFNET_graph[[1]], features, target,
        niter = niter, ntrees = ntrees, initial.walk.depth = init.mtry
    )
    return(list(
        DFNET_graph = DFNET_graph, DFNET_trees = forest$trees,
        DFNET_MODULES = forest$modules,
        DFNET_MODULES_AUC = attr(forest, "last.performance")
    ))
}
