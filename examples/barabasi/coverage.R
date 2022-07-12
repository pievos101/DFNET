library(DFNET)
source("graph.R")

input <- NULL
while (is.null(input)) {
    tryCatch(
        {
            input <- simple_input_graph(n.nodes = 50)
        },
        error = function(e) {
            message("failed to create a graph, retrying...")
        }
    )
}

n.forests <- 50
niter <- c(2, 3, 5, 10, 30, 50)
cum_niter <- c(0, cumsum(niter))

## message("allocating forests, this may take some time...")
forests <- lapply(1:n.forests, function(x) NULL)

coverage <- matrix(NaN, n.forests, length(niter))
iou <- matrix(NaN, n.forests, length(niter))
auc <- matrix(NaN, n.forests, length(niter))
n_candidates <- matrix(NaN, n.forests, length(niter))

intersection_over_union <- function(as, bs) {
    length(intersect(as, bs)) / length(union(as, bs))
}

for (xx in 1:length(niter)) {
    message(xx, " of ", length(niter), " done! ----------------")

    for (yy in 1:n.forests) {
        forests[[yy]] <- train(
            forests[[yy]],
            input$graph, input$features, input$target,
            niter = niter[xx], offset = cum_niter[xx]
        )
        # keep only the last generation
        # XXX: Why do we have to qualify DFNET.forest here???
        forests[[yy]] <- tail(forests[[yy]], 1)

        tree_imp <- attr(forests[[yy]], "last.performance")
        edge_imp <- edge_importance(input$graph, forests[[yy]]$trees, tree_imp)
        umi <- unique_module_importance(
            input$graph, forests[[yy]]$modules, edge_imp, tree_imp
        )
        by_importance <- order(umi$table[, "total"], decreasing = TRUE)
        top_modules <- umi$modules[by_importance]

        sm <- unique(sort(input$module))
        bm <- unique(sort(top_modules[[1]]))

        coverage[yy, xx] <- as.numeric(identical(sm, bm))
        iou[yy, xx] <- intersection_over_union(sm, bm)
        auc[yy, xx] <- mean(tree_imp)
        n_candidates[yy, xx] <- length(top_modules)
    }
    print(coverage)
    print(n_candidates)
}

# Save the workspace
save.image("results.RData")
save(input, file = "input.RData")

# Do some plots
xnames <- cum_niter[2:length(cum_niter)]
par(mfrow = c(2, 2))
barplot(
    apply(coverage, 2, sum) / n.forests,
    ylab = "Identical modules", xlab = "Iterations",
    names.arg = xnames, col = "cadetblue", ylim = c(0, 1), las = 2
)

boxplot(
    iou,
    ylab = "Intersection over union", xlab = "Iterations",
    names = xnames, col = "cadetblue", ylim = c(0, 1), las = 2
)

boxplot(
    n_candidates,
    ylab = "Number of unique modules detected", xlab = "Iterations",
    names = xnames, col = "cadetblue", las = 2
)

boxplot(
    auc,
    ylab = "Mean area under curve", xlab = "Iterations",
    names = xnames, col = "cadetblue", ylim = c(0, 1), las = 2
)
