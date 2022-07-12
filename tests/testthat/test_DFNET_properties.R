## test_DFNET_properties - Property-based testing on DFNET.

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

library(hedgehog)

# Training takes time, so let's only run few tests unless the user wants more
if (is.null(getOption("hedgehog.tests"))) {
    options(hedgehog.tests = 20)
}

## Utility functions
sample.graph <- function(n.nodes, power, n.samples, n.features) {
    g <- igraph::sample_pa(n.nodes, power = power, directed = FALSE)
    features <- NULL
    if (n.features == 0) {
        features <- matrix(
            sample(0:1, n.samples * n.nodes, replace = TRUE),
            n.samples, n.nodes
        )
        colnames(features) <- c(paste("N_", 1:n.nodes, sep = ""))
    } else {
        features <- array(
            sample(0:1, n.samples * n.nodes * n.features, replace = TRUE),
            c(n.samples, n.nodes, n.features),
            dimnames = list(NULL, c(paste("N_", 1:n.nodes, sep = ""), NULL))
        )
    }
    list(graph = g, features = features)
}

performance <- function(forest) attr(forest, "last.performance")
module_weights <- function(forest) as.numeric(lapply(forest$modules.weights, sum))

## Generators
gen.graph <-
    gen.bind(
        function(params) {
            gen.impure(function(n) {
                sample.graph(
                    params$n.nodes,
                    params$power,
                    params$n.samples,
                    params$n.features
                )
            })
        },
        gen.structure(
            list(
                gen.element(seq(10, 30)),
                gen.element(seq(0.5, 3, by = 0.1)),
                gen.element(seq(100, 300, by = 10)),
                gen.element(0:3)
            ),
            names = c("n.nodes", "power", "n.samples", "n.features")
        )
    )

gen.graph_and_target <-
    gen.and_then(
        gen.graph,
        function(gf) {
            gen.impure(function(n) {
                target <- sample(0:1, dim(gf$features)[[1]], replace = TRUE)
                list(
                    graph = gf$graph,
                    features = gf$features,
                    target = target
                )
            })
        }
    )

gen.test_metric <- gen.choice(
    ModelMetrics::auc,
    ModelMetrics::precision,
    ModelMetrics::recall,
    ModelMetrics::f1Score
)

# Test parameters, such as number of trees or iterations
gen.test_niter <-
    gen.element(getOption("test_DFNET_properties.niter", 10:30))
gen.test_small_niter <-
    gen.element(getOption("test_DFNET_properties.small_niter", 5:10))
gen.test_c_small_niter <-
    do.call(
        gen.c,
        c(
            list(generator = gen.test_small_niter),
            getOption(
                "test_DFNET_properties.n_small_niter",
                list(from = 2, to = 4)
            )
        )
    )

gen.test_ntrees <-
    gen.element(
        getOption(
            "test_DFNET_properties.ntrees",
            c(2, 4, 8, 16, 32, 64)
        )
    )

test_flaky <- getOption("test_DFNET_properties.flaky", FALSE)

## Actual tests
test_that("No iteration means no iteration", {
    forall(
        list(
            gf = gen.graph_and_target,
            ntrees = gen.test_ntrees
        ),
        function(gf, niter, ntrees) {
            graph <- gf$graph
            features <- gf$features
            target <- gf$target

            forest <- train(, graph, features, target, 0, ntrees = ntrees)

            expect_equal(head(forest)$modules, tail(forest)$modules)
        }
    )
})

test_that("DFNET improves performance", {
    forall(
        list(
            gf = gen.graph_and_target,
            niter = gen.test_niter,
            ntrees = gen.test_ntrees,
            metric = gen.test_metric,
            do.predict = gen.choice(TRUE, FALSE)
        ),
        function(gf, niter, ntrees, metric, do.predict) {
            graph <- gf$graph
            features <- gf$features
            target <- gf$target

            p <- tester(features, target, metric, do.predict)

            state0 <- train(
                NULL, graph, features, target, 0,
                ntrees = ntrees, performance = p
            )
            state1 <- train(
                state0, graph, features, target, niter,
                performance = p
            )

            expect_true(
                all(performance(state0) <= performance(state1)),
                label = "no worse after iteration"
            )
            # Theoretically, this may fail if DFNET_init creates an optimal
            # configuration (unlikely).  In addition, some metrics appear
            # harder to optimize than others.
            if (test_flaky) {
                expect_true(
                    any(performance(state0) < performance(state1)),
                    label = "better after iteration"
                )
            }
        }
    )
})

test_that("DFNET performance uses target metric", {
    forall(
        list(
            gf = gen.graph_and_target,
            niter = gen.test_niter,
            ntrees = gen.test_ntrees,
            metric = gen.test_metric,
            do.predict = gen.choice(TRUE, FALSE)
        ),
        function(gf, niter, ntrees, metric, do.predict) {
            graph <- gf$graph
            features <- gf$features
            target <- gf$target

            p <- tester(features, target, metric, do.predict)

            state0 <- train(
                NULL, graph, features, target, 0,
                ntrees = ntrees, performance = p
            )
            state1 <- train(
                state0, graph, features, target, niter,
                performance = p
            )
            state1 <- tail(state1, 1)

            expect_identical(performance(state0), sapply(state0$trees, p))
            expect_identical(performance(state1), sapply(state1$trees, p))
        }
    )
})

test_that("DFNET adds up", {
    forall(
        list(
            gf = gen.graph_and_target,
            niter = gen.test_c_small_niter,
            ntrees = gen.test_ntrees,
            metric = gen.test_metric,
            do.predict = gen.choice(FALSE) # XXX: predict() invokes randomness
        ),
        function(gf, niter, ntrees, metric, do.predict) {
            graph <- gf$graph
            features <- gf$features
            target <- gf$target

            p <- tester(features, target, metric, do.predict)

            state0 <- train(
                NULL, graph, features, target, 0,
                ntrees = ntrees,
                performance = p
            )
            state1 <- state0

            saved.seed <- .Random.seed
            for (iter in niter) {
                state0 <- train(
                    state0, graph, features, target, iter,
                    performance = p
                )
            }
            assign(".Random.seed", saved.seed, envir = globalenv())
            state1 <- train(
                state1, graph, features, target, sum(niter),
                performance = p
            )
            expect_equal(state0$modules, state1$modules)
            expect_equal(performance(state0), performance(state1))
        }
    )
})

test_that("DFNET returns simplified modules", {
    forall(
        list(
            gf = gen.graph_and_target,
            niter = gen.choice(0, gen.test_niter),
            ntrees = gen.test_ntrees
        ),
        function(gf, niter, ntrees) {
            graph <- gf$graph
            features <- gf$features
            target <- gf$target

            forest <- train(, graph, features, target, niter, ntrees = ntrees)

            expect_identical(
                forest$modules,
                lapply(forest$modules, function(m) unique(sort(m)))
            )
        }
    )
})

test_that("DFNET shrinks module weights", {
    forall(
        list(
            gf = gen.graph_and_target,
            niter = gen.test_niter,
            ntrees = gen.test_ntrees
        ),
        function(gf, niter, ntrees) {
            graph <- gf$graph
            features <- gf$features
            target <- gf$target

            forest <- train(, graph, features, target, niter, ntrees = ntrees)

            weights0 <- module_weights(head(forest, 1))
            weights1 <- module_weights(tail(forest, 1))

            expect_true(
                all(weights0 >= weights1),
                label = "no worse after iteration"
            )
            # As with the other variant, this test is very flaky
            if (test_flaky) {
                expect_true(
                    any(weights0 > weights1),
                    label = "better after iteration"
                )
            }
        }
    )
})

test_that("DFNET predicts data", {
    forall(
        list(
            gf = gen.graph_and_target,
            niter = gen.test_niter,
            ntrees = gen.test_ntrees
        ),
        function(gf, niter, ntrees) {
            graph <- gf$graph
            features <- gf$features
            target <- gf$target

            forder <- sample(dim(features)[1])
            train_ids <- head(forder, floor(length(forder) * 0.8))
            test_ids <- tail(forder, -floor(length(forder) * 0.8))

            if (length(dim(features)) == 2) {
                train_features <- features[train_ids, ]
                test_features <- features[test_ids, ]
            } else {
                train_features <- features[train_ids, , ]
                test_features <- features[test_ids, , ]
            }

            forest <- train(
                NULL, graph,
                train_features, target[train_ids],
                niter,
                ntrees = ntrees
            )

            prediction <- predict(forest, test_features)
            expect_equal(length(prediction), length(target[test_ids]))
            given.labels <- unique(prediction)
            allowed.labels <- sort(unique(c(NaN, target)))
            expect_equal(
                sort(union(given.labels, allowed.labels)),
                allowed.labels
            )
        }
    )
})
