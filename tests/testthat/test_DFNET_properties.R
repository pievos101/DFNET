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

# DFNET_iter takes time, so let's only run few tests unless the user wants more
if (is.null(getOption("hedgehog.tests"))) {
    options(hedgehog.tests = 20)
}

## Utility functions
sample.graph <- function(n.nodes, power, n.samples) {
    g <- sample_pa(n.nodes, power = power, directed = FALSE)
    features <- matrix(
        sample(0:1, n.samples * n.nodes, replace = TRUE),
        n.samples, n.nodes
    )
    colnames(features) <- c(paste("N_", 1:n.nodes, sep = ""))
    list(graph = g, features = features)
}

performance <- function(forest) attr(forest, "last.performance")
module_size <- function(forest) as.numeric(lapply(forest$modules, length))

## Generators
gen.graph <-
    gen.bind(
        function(params) {
            gen.impure(function(n) {
                sample.graph(
                    params$n.nodes,
                    params$power,
                    params$n.samples
                )
            })
        },
        gen.structure(
            list(
                gen.element(seq(10, 30)),
                gen.element(seq(0.5, 3, by = 0.1)),
                gen.element(seq(100, 300, by = 10))
            ),
            names = c("n.nodes", "power", "n.samples")
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

            state0 <- DFNET_init(graph, features, target, ntrees = ntrees)
            state1 <- DFNET_iterate(state0, graph, features, target, 0)

            expect_equal(state0$modules, state1$modules)
        }
    )
})

test_that("DFNET improves performance", {
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

            state0 <- DFNET_init(graph, features, target, ntrees = ntrees)
            state1 <- DFNET_iterate(state0, graph, features, target, niter)

            expect_true(
                all(performance(state0) <= performance(state1)),
                label = "no worse after iteration"
            )
            # Theoretically, this may fail if DFNET_init creates an optimal
            # configuration (unlikely)
            expect_true(
                any(performance(state0) < performance(state1)),
                label = "better after iteration"
            )
        }
    )
})

test_that("DFNET adds up", {
    forall(
        list(
            gf = gen.graph_and_target,
            niter = gen.test_c_small_niter,
            ntrees = gen.test_ntrees
        ),
        function(gf, niter, ntrees) {
            graph <- gf$graph
            features <- gf$features
            target <- gf$target

            state0 <- DFNET_init(graph, features, target, ntrees = ntrees)
            state1 <- state0

            saved.seed <- .Random.seed
            for (iter in niter) {
                state0 <- DFNET_iterate(state0, graph, features, target, iter)
            }
            assign(".Random.seed", saved.seed, envir = globalenv())
            state1 <- DFNET_iterate(state1, graph, features, target, sum(niter))
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

            forest <- DFNET_init(graph, features, target, ntrees = ntrees)
            forest <- DFNET_iterate(forest, graph, features, target, niter)

            expect_identical(
                forest$modules,
                lapply(forest$modules, function(m) unique(sort(m)))
            )
        }
    )
})

test_that("DFNET may prune trees", {
    forall(
        list(
            gf = gen.graph_and_target,
            niter = gen.test_niter,
            ntrees = gen.test_ntrees,
            keep = gen.test_niter
        ),
        function(gf, niter, ntrees, keep) {
            graph <- gf$graph
            features <- gf$features
            target <- gf$target

            forest <- DFNET_init(graph, features, target, ntrees = ntrees)
            forest <- DFNET_iterate(
                forest, graph, features, target, niter,
                keep.generations = keep
            )

            # XXX: We only validate lengths, not selections
            expect_lte(length(forest$modules), ntrees * keep)
            expect_lte(length(forest$trees), ntrees * keep)
            expect_equal(length(forest$modules), length(forest$trees))
        }
    )
})
