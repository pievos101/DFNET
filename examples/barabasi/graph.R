library(igraph)

select_module <- function(graph, size = 4) {
    module <- rep.int(0, size)
    module[1] <- sample(V(graph), 1, prob = sapply(neighborhood(graph), length))
    if (size > 1) {
        N <- as.numeric(neighbors(graph, module[1]))
        for (i in 2:size) {
            module[i] <- N[sample(length(N), 1)]
            N <- setdiff(c(N, neighbors(graph, module[i])), module[1:i])
        }
    }
    return(module)
}

generate_features <- function(graph, n.samples = 1000, n.modes = 1, range = c(0, 1)) {
    n.nodes <- length(V(graph))
    data <- sample(range, n.nodes * n.samples * n.modes, replace = TRUE)
    features <- array(
        data,
        dim = c(n.samples, n.nodes, n.modes),
        dimnames = list(NULL, vertex_attr(graph, "name"), NULL)
    )
    return(features)
}

select_features <- function(features, module) {
    d <- dim(features)
    subset <- matrix(0, nrow = d[1], length(module))
    for (i in 1:length(module)) {
        subset[, i] <- features[, module[i], sample(d[3], size = 1)]
    }
    return(subset)
}

simple_input_graph <- function(n.nodes = 30, power = 1.2, n.samples = 1000,
                               n.modes = 1) {
    g <- sample_pa(n.nodes, power = power, directed = FALSE)
    vertex_attr(g, "name") <- paste("N_", 1:length(V(g)), sep = "")
    module <- select_module(g, 4)
    features <- generate_features(g, n.samples, n.modes)
    t0 <- select_features(features, module)
    target <- as.numeric(xor(t0[, 1] & t0[, 2], t0[, 3] & t0[, 4]))
    return(
        list(
            graph = g,
            features = features,
            target = target,
            module = module
        )
    )
}
