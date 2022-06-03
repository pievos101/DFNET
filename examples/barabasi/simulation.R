library(DFNET)

source("graph.R")

input <- NULL
while(is.null(input)) {
    tryCatch({
        input <- simple_input_graph(n.nodes=20)
    }, error = function(e) {
        message("failed to create a graph, retrying...")
    })
}

graph <- input$graph
features <- input$features
target <- input$target

forest <- DFNET_init(graph, features, target, ntrees = 100)
forest <- DFNET_iterate(forest, graph, features, target, niter = 10)

# XXX: public metrics
tree_imp <- sapply(
    forest$trees,
    function(tree) DFNET:::area_under_curve(tree$predictions, target)
)
edge_imp <- edge_importance(graph, forest$trees, tree_imp)
print(edge_imp)
module_imp <- module_importance(graph, forest$modules, edge_imp, tree_imp)

modules <- sapply(forest$modules, function(module) {
    paste(unique(sort(module)), collapse = " ")
})

data <- data.frame(modules, tree_imp, module_imp)
by_importance <- order(data[,"total"], decreasing = TRUE)

detected_module <- unlist(forest$modules[by_importance][1])

print(head(data[by_importance,]))
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
    edge.width=exp(edge_imp + 1), edge.color = color,
    vertex.shape = "none",
    # use numeric labels again
    vertex.label = 1:length(V(graph)),
    vertex.label.color = "black"
)