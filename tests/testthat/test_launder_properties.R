## test_DFNET_properties - Property-based testing on DFNET laundering.

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
library(stringr) # str_c

sample.graph <- function(names, power) {
    names <- unique(names) # let's be cautious
    g <- sample_pa(length(names), power = power, directed = FALSE)
    vertex_attr(g, "names") <- names
    return(g)
}

sample.features <- function(columns, n.samples) {
    features <- matrix(
        sample(0:1, n.samples * length(columns), replace = TRUE),
        n.samples, length(columns)
    )
    colnames(features) <- columns
    return(features)
}

gen.s <- gen.map(
    function (x) do.call(str_c, x),
    gen.list(gen.element(letters))
)

test_that("common_features are common", {
    forall(
        list(
            common = gen.c(gen.s, from=2),
            uncommon = gen.list(gen.c(gen.s), from=2),
            n.samples = gen.element(seq(100, 300, by = 10))
        ),
        function(common, uncommon, n.samples) {
            common <- unique(common)
            # matrices don't have column names if there's only one column
            if(length(common) < 2) discard()
            matrices <- lapply(uncommon, function(x) {
                sample.features(unique(c(common, x)), n.samples)
            })
            common.features <- common_features(matrices)
            expect_equal(
                # We might actually have more common features, if uncommon
                # is shaped badly.  We only test the common ones we know of,
                # because it makes no sense to discard these tests.
                intersect(dimnames(common.features)[[2]], common),
                common
            )
            for (xx in 1:length(matrices)) {
                expect_equal(
                    # Use known common again as subscript due to potentially
                    # unknown commons existing.  See the comment above.
                    common.features[,common,xx],
                    matrices[[xx]][,common]
                )
            }
        }
    )
})

test_that("graphed_features are graphed", {
    forall(
        list(
            common = gen.c(gen.s, from=2),
            ungraphed = gen.c(gen.s),
            n.samples = gen.element(seq(100, 300, by = 10)),
            power = gen.element(seq(0.5, 3, by = 0.1))

        ),
        function(common, ungraphed, n.samples, power) {
            graph <- sample.graph(common, power)
            if(length(V(graph)) < 2) discard()
            features <- sample.features(
                sample(unique(c(common, ungraphed))),
                n.samples
            )
            graphed.features <- graphed_features(features, graph)

            expect_equal(
                sort(colnames(graphed.features)),
                sort(vertex_attr(graph, "names"))
            )
        }
    )
})

test_that("launder washes data", {
    forall(
        list(
            common = gen.c(gen.s, from=2),
            graph_dummies = gen.c(gen.s),
            uncommon = gen.list(gen.c(gen.s), from=2),
            n.samples = gen.element(seq(100, 300, by = 10)),
            power = gen.element(seq(0.5, 3, by = 0.1))

        ),
        function(common, graph_dummies, uncommon, n.samples, power) {
            common <- unique(common)
            # matrices don't have column names if there's only one column
            if(length(common) < 2) discard()

            graph <- sample.graph(c(common, graph_dummies), power)
            matrices <- lapply(uncommon, function(x) {
                sample.features(unique(c(common, x)), n.samples)
            })
            laundered = launder(graph, matrices)

            expect_equal(
                dimnames(laundered$features)[[2]],
                vertex_attr(laundered$graph, "names"),
                "laundered graph and features match"
            )
            expect_equal(
                intersect(dimnames(laundered$features)[[2]], common),
                common,
                "known commons are found"
            )
        }
    )
})
