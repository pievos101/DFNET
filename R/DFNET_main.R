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

#' Use random forest algorithm to perform feature selection with respect to
#' given modules.
#'
#' @param modules a list of modules (node vectors)
#' @param features the (multi-modal) features to train on
#' @param target the target vector
#' @return a decision forest with one tree per module in modules
DFNET_make_forest <- function(modules, features, target) {
    modules_rle <- lapply(modules, function(m) rle(sort(m)))

    decision_trees <- lapply(modules_rle, function(m) {
        unique_nodes <- m$values
        unique_nodes_weights <- m$lengths
        weights <- unique_nodes_weights

        if (length(dim(features)) == 2) {
            mm_data <- as.matrix(features[, unique_nodes])
        } else if (dim(features)[3] == 1) {
            # unimodal data in 3D array, would otherwise collapse
            mm_data <- as.matrix(features[, unique_nodes, 1])
        } else {
            mm_data <- features[, unique_nodes, ]
            d <- dim(mm_data)
            # XXX: repeats colnames per mode, strips mode name
            dimnames <- list(
                dimnames(mm_data)[[1]],
                rep.int(dimnames(mm_data)[[2]], d[3])
            )
            mm_data <- matrix(mm_data, d[1], d[2] * d[3], dimnames = dimnames)
            weights <- rep.int(weights, d[3])
        }

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

#' Initialize the decision forest network.
#'
#' @param graph the graph to train the network on
#' @param features the features to use for training
#' @param ntrees how many trees should be generated in the initial step
#' and each following iteration
#' @param walk.depth how many nodes should be selected per tree.
#' If not a number, \code{ceiling(sqrt(length(V(graph))))} will be used instead.
#' @param performance a function to call with the predictions and target to
#' estimate the performance of a decision tree.
#' @return an initialized \code{DFNET::forest}.
#'
#' A \code{DFNET::forest} is a list of shape \code{trees, modules},
#' where \code{trees} are the decision trees created for detected \code{modules},
#' and \code{performance} is their performance as per the input metric.
#'
#' As "private" attributes used for iteration, \code{generation_size} is set to
#' \code{ntrees}, \code{walk.depth} to \code{walk.depth}, and
#' \code{last.performance} to a vector of length \code{ntrees}, containing the
#' result of \code{performance} of each tree w.r.t. \code{target}.
DFNET_init <- function(graph, features, target,
                       ntrees = 100, walk.depth = NaN,
                       performance = DFNET:::area_under_curve) {
    nodes <- V(graph)
    n.nodes <- length(nodes)

    if (is.na(walk.depth)) {
        walk.depth <- ceiling(sqrt(n.nodes))
    }
    walk.depth <- rep_len(walk.depth, ntrees)

    count <- 1
    selected_nodes <- list()
    repeat {
        sampled.nodes <- sample(nodes, (ntrees + 1 - count), replace = TRUE)
        for (sn in sampled.nodes) {
            depth <- walk.depth[[count]]
            selected_nodes[[count]] <- as.numeric(random_walk(graph, sn, depth))
            # Pick only walks of maximal length
            if (length(selected_nodes[[count]]) >= walk.depth[[count]]) {
                count <- count + 1
            }
        }
        if (count > ntrees) {
            break
        }
    }

    seed <- DFNET_make_forest(selected_nodes, features, target)
    last.perf <- sapply(seed$trees, function(tree) {
        performance(tree$predictions, target)
    })

    return(
        structure(
            seed,
            class = "DFNET::forest",
            generation_size = ntrees,
            walk.depth = walk.depth,
            last.performance = last.perf
        )
    )
}

#' Perform iterations on the decision forest network.
#'
#' @param forest a state as generated by \code{DFNET_init} or
#' \code{DFNET_iterate}.
#' @param graph the graph to train the forest on
#' @param features the features to train the forest on
#' @param target the target vector
#' @param niter the number of iterations to run.
#' @param offset an offset for the iteration count (used for logging only)
#' @param min.walk_depth the minimal random walk depth in each iteration
#' @param performance a function to call with the predictions and target to
#' estimate the performance of a decision tree.
#' @param keep.generations the number of generations to keep after iterating,
#' defaults to all generations
#' @return the updated \code{DFNET::forest}.
#' @examples
#' \dontrun{
#' forest <- DFNET_init(graph, features, ...)
#' offset <- 0
#' while (keep_iterating(forest, target)) { # insert your own iteration criteria
#'     forest <- DFNET_iterate(
#'         forest,
#'         graph,
#'         features,
#'         niter = 10,
#'         offset = offset
#'     )
#'     offset <- offset + 10
#' }
#' }
#'
DFNET_iterate <- function(forest, graph, features, target,
                          niter = 200, offset = 0, min.walk_depth = 2,
                          performance = DFNET:::area_under_curve,
                          keep.generations = NA) {
    if (offset < 0 || niter < 0) {
        stop("both offset and niter must be positive")
    } else if (niter == 0) {
        return(forest)
    }

    if (!is.na(keep.generations) && keep.generations < 1) {
        stop("need to keep at least one generation")
    }

    ntrees <- attr(forest, "generation_size")
    last.walk.depth <- attr(forest, "walk.depth")

    all.trees <- forest$trees
    all.modules <- forest$modules
    all.modules.weights <- forest$modules.weights
    last.trees <- tail(all.trees, ntrees)
    last.modules <- tail(all.modules, ntrees)
    last.modules.weights <- tail(all.modules.weights, ntrees)

    tree_performance <- function(tree) performance(tree$predictions, target)
    last.perf <- sapply(last.trees, tree_performance)

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
            as.numeric(random_walk(graph, start_nodes[sn], walk.depth[sn]))
        })

        next_gen <- DFNET_make_forest(modules, features, target)
        perf <- sapply(next_gen$trees, tree_performance)

        good_enough <- perf >= last.perf

        # Update inner state
        last.perf <- ifelse(good_enough, perf, last.perf)
        last.modules <- ifelse(good_enough, next_gen$modules, last.modules)
        last.modules.weights <-
            ifelse(good_enough, next_gen$modules.weights, last.modules.weights)
        last.trees <- ifelse(good_enough, next_gen$trees, last.trees)
        last.walk.depth <- ifelse(good_enough, walk.depth - 1, last.walk.depth)
        last.walk.depth[last.walk.depth < min.walk_depth] <- min.walk_depth

        if (is.na(keep.generations) || (keep.generations > 1)) {
            all.modules <- c(all.modules, last.modules)
            all.modules.weights <- c(all.modules.weights, last.modules.weights)
            all.trees <- c(all.trees, last.trees)
        } else {
            all.modules <- last.modules
            all.modules.weights <- last.modules.weights
            all.trees <- last.trees
        }
    }

    if (!is.na(keep.generations)) {
        nkeep <- keep.generations * ntrees

        all.modules <- tail(all.modules, nkeep)
        all.modules.weights <- tail(all.modules.weights, nkeep)
        all.trees <- tail(all.trees, nkeep)
    }

    return(
        structure(
            list(
                trees = all.trees,
                modules = all.modules,
                modules.weights = all.modules.weights
            ),
            class = "DFNET::forest",
            generation_size = ntrees,
            walk.depth = last.walk.depth,
            last.performance = last.perf
        )
    )
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

    forest <- DFNET_init(
        DFNET_graph[[1]], features, target,
        ntrees = ntrees, walk.depth = init.mtry
    )
    forest <- DFNET_iterate(forest, DFNET_graph[[1]], features, target, niter = niter)
    return(list(
        DFNET_graph = DFNET_graph, DFNET_trees = forest$trees,
        DFNET_MODULES = forest$modules,
        DFNET_MODULES_AUC = attr(forest, "last.performance")
    ))
}
