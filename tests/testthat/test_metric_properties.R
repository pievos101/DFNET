## test_metric_properties - Property-based testing on DFNET metrics.

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

if (is.null(getOption("hedgehog.tests"))) {
    options(hedgehog.tests = 20)
}

sample.graph <- function(n.nodes, power, n.modules) {
    g <- igraph::sample_pa(n.nodes, power = power, directed = FALSE)
    igraph::vertex_attr(g, "names") <- paste("N_", 1:n.nodes)

    modules.init <- sample(V(g), n.modules, replace = TRUE)
    modules.size <- sample(
        ceiling(sqrt(length(V(g)))), n.modules, replace = TRUE
    )
    modules.size[modules.size == 1] = 2

    modules <- apply(cbind(modules.init, modules.size), 1, function(m) {
        unique(sort(as.numeric(igraph::random_walk(g, m[1], m[2]))))
    })

    list(graph = g, modules = modules)
}

sample.fake_trees <- function(graph, modules) {
    trees <- lapply(modules, function(m) {
        varimp <- sample(seq(0, 1, by = 0.01), size = length(m), replace = TRUE)
        names(varimp) <- igraph::vertex_attr(graph, "names")[m]
        list(variable.importance = varimp)
    })

    timp <- sample(seq(0, 1, by = 0.01), size = length(trees), replace = TRUE)

    list(
        graph = graph, modules = modules, trees = trees, tree_importances = timp
    )
}

sample.features <- function(graph, modules, n.samples, n.features) {
    if (n.features == 0) {
        features <- matrix(
            sample(0:1, n.samples * length(V(graph)), replace = TRUE),
            n.samples, length(V(graph))
        )
        colnames(features) <- igraph::vertex_attr(graph, "names")
    } else {
        features <- array(
            sample(0:1, n.samples * length(V(graph)) * n.features, replace = TRUE),
            c(n.samples, length(V(graph)), n.features),
            dimnames = list(NULL, igraph::vertex_attr(graph, "names"), NULL)
        )
    }
    target <- sample(0:1, n.samples, replace = TRUE)

    return(list(
        graph = graph, modules = modules, features = features, target = target
    ))
}

gen.fake_trees <-
    gen.bind(
        function(params) {
            gen.impure(function(n) {
                gm <- sample.graph(
                    params$n.nodes, params$power, params$n.modules
                )
                sample.fake_trees(gm$graph, gm$modules)
            })
        },
        gen.structure(
            list(
                gen.element(seq(10, 30)),
                gen.element(seq(0.5, 3, by = 0.1)),
                gen.element(c(2, 4, 8, 16, 32, 64))
            ),
            names = c("n.nodes", "power", "n.modules")
        )
    )

gen.fake_features <-
    gen.bind(
        function(params) {
            gen.impure(function(n) {
                gm <- sample.graph(
                    params$n.nodes, params$power, params$n.modules
                )
                sample.features(
                    gm$graph, gm$modules,
                    params$n.samples, params$n.features
                )
            })
        },
        gen.structure(
            list(
                gen.element(seq(10, 30)),
                gen.element(seq(0.5, 3, by = 0.1)),
                gen.element(c(2, 4, 8, 16, 32, 64)),
                gen.element(seq(100, 300, by = 10)),
                gen.element(0:3)
            ),
            names = c(
                "n.nodes", "power", "n.modules", "n.samples", "n.features"
            )
        )
    )

test_that("feature_importance has the right shape", {
    forall(gen.fake_features, function(gt) {
        for (m in gt$modules) {
            if (length(m) < 2) discard()
        }

        d <- DFNET:::learn_decisions(gt$modules, gt$features, gt$target)
        fimp <- feature_importance(d, gt$features)
        expect_equal(dim(fimp), c(dim(gt$features), 1)[2:3])
    })
})

test_that("edge importance is relative", {
    forall(gen.fake_trees, function(gt) {
        graph <- gt$graph
        trees <- gt$trees
        tree_importances <- gt$tree_importances

        eimp <- edge_importance(graph, trees, tree_importances)
        expect_true(all(eimp >= 0))
        expect_true(all(eimp <= 1))
    })
})

test_that("module_importance sums up", {
    forall(gen.fake_trees, function(gt) {
        graph <- gt$graph
        trees <- gt$trees
        tree_importances <- gt$tree_importances

        eimp <- edge_importance(graph, trees, tree_importances)
        mimp <- module_importance(graph, gt$modules, eimp, tree_importances)

        expect_equal(mimp[, 2], tree_importances + mimp[, 1])
    })
})
