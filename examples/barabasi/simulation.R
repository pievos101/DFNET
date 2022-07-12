library(DFNET)

source("graph.R")

input <- NULL
while (is.null(input)) {
    tryCatch(
        {
            input <- simple_input_graph(n.nodes = 20)
        },
        error = function(e) {
            message("failed to create a graph, retrying...")
        }
    )
}

graph <- input$graph
features <- input$features
target <- input$target

forest <- train(, graph, features, target, ntrees = 100, niter = 10)

# XXX: public metrics
tree_imp <- sapply(
    forest$trees,
    function(tree) ModelMetrics::auc(target, tree$predictions)
)
edge_imp <- edge_importance(graph, forest$trees, tree_imp)
print(edge_imp)
umi <- unique_module_importance(graph, forest$modules, edge_imp, tree_imp)

by_importance <- order(umi$table[, "total"], decreasing = TRUE)
detected_module <- unlist(umi$modules[by_importance][1])

print(head(umi$table[by_importance, ]))
print(unique(sort(input$module)))

edgelist <- as_edgelist(graph, names = FALSE)
color <- rep("cadetblue", dim(edgelist)[1])

for (e in 1:dim(edgelist)[1]) {
    edge <- edgelist[e, 1:2]
    ids <- match(edge, detected_module)

    if (!any(is.na(ids))) {
        color[e] <- "coral2"
    }
}

plot(
    graph,
    edge.width = exp(edge_imp + 1), edge.color = color,
    vertex.shape = "none",
    # use numeric labels again
    vertex.label = 1:length(V(graph)),
    vertex.label.color = "black"
)
