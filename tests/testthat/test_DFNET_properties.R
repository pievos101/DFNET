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
                    features = as.data.frame(cbind(gf$features, target))
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

            state0 <- DFNET_init(graph, features, ntrees = ntrees)
            state1 <- DFNET_iterate(state0, graph, features, 0)

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

            state0 <- DFNET_init(graph, features, ntrees = ntrees)
            state1 <- DFNET_iterate(state0, graph, features, niter)

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

            state0 <- DFNET_init(graph, features, ntrees = ntrees)
            state1 <- state0

            saved.seed <- .Random.seed
            for (iter in niter) {
                state0 <- DFNET_iterate(state0, graph, features, iter)
            }
            assign(".Random.seed", saved.seed, envir = globalenv())
            state1 <- DFNET_iterate(state1, graph, features, sum(niter))
            expect_equal(state0$modules, state1$modules)
            expect_equal(performance(state0), performance(state1))
        }
    )
})

test_that("DFNET shrinks modules", {
    forall(
        list(
            gf = gen.graph_and_target,
            niter = gen.test_niter,
            ntrees = gen.test_ntrees
        ),
        function(gf, niter, ntrees) {
            graph <- gf$graph
            features <- gf$features

            state0 <- DFNET_init(graph, features, ntrees = ntrees)
            state1 <- DFNET_iterate(
                state0, graph, features, niter,
                keep.generations = 1
            )

            expect_true(
                all(module_size(state1) <= module_size(state0)),
                label = "no worse after iteration"
            )
            # Shaky, might fail if DFNET_init creates a near-optimal
            # solution or there are too few iterations
            expect_true(
                any(module_size(state1) < module_size(state0)),
                label = "better after iteration"
            )
        }
    )
})

test_that("DFNET shrinks modules (manual selection)", {
    forall(
        list(
            gf = gen.graph_and_target,
            niter = gen.test_niter,
            ntrees = gen.test_ntrees
        ),
        function(gf, niter, ntrees) {
            graph <- gf$graph
            features <- gf$features

            forest <- DFNET_init(graph, features, ntrees = ntrees)
            forest <- DFNET_iterate(forest, graph, features, niter)

            first <- head(forest$modules, ntrees)
            last <- tail(forest$modules, ntrees)

            expect_true(
                all(sapply(last, length) <= sapply(first, length)),
                label = "no worse after iteration"
            )
            # Shaky, might fail if DFNET_init creates a near-optimal
            # solution or there are too few iterations
            expect_true(
                any(sapply(last, length) < sapply(first, length)),
                label = "better after iteration"
            )
        }
    )
})