# Generate Barabasi Graph and Festures
library(igraph)
DFNET_generate_graph_SM <- function(N.Nodes = 30, power = 1.2, n.samples = 1000) {
    g <- barabasi.game(N.Nodes,
        power = power, m = NULL, out.dist = NULL,
        out.seq = NULL, out.pref = FALSE, zero.appeal = 1, directed = FALSE,
        algorithm = "psumtree", start.graph = NULL
    )


    EDGELIST <- as_edgelist(g, names = TRUE)

    # Select an XOR edge
    id <- sample(1:dim(EDGELIST)[1], 1)

    node1 <- EDGELIST[id, 1]
    node2 <- EDGELIST[id, 2]

    ## test
    id <- which(EDGELIST[, 1] == node1)
    if (length(id) < 3) {
        print("retry ... done!")
        return(NULL)
    } else {
        node1 <- EDGELIST[id[1], 1]
        node2 <- EDGELIST[id[1], 2]
        node3 <- EDGELIST[id[2], 2]
        node4 <- EDGELIST[id[3], 2]
    }

    # Generate input data
    MAT <- matrix(sample(c(0, 1), n.samples * N.Nodes, replace = TRUE), n.samples, N.Nodes)

    # Define the response
    # target <- as.numeric(xor(MAT[,node1], MAT[,node2]))
    target <- as.numeric(xor(MAT[, node1] & MAT[, node2], MAT[, node3] & MAT[, node4]))
    # target <- as.numeric(xor(xor(MAT[,node1], MAT[,node2]), MAT[,node3]))

    IN <- cbind(MAT, target)
    IN <- as.data.frame(IN)
    colnames(IN) <- c(paste("N_", 1:N.Nodes, sep = ""), "target")


    return(list(graph = g, Feature_Matrix = IN, Selected_Module = c(node1, node2, node3, node4)))
}
