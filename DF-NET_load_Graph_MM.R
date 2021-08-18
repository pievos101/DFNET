# Generate Barabasi Graph and Festures
library(igraph)
DFNET_load_graph_MM <- function(N.Nodes=30, power=1.2, n.samples=1000){

#g <- barabasi.game(N.Nodes, power = power, m = NULL, out.dist = NULL, 
#     	out.seq = NULL, out.pref = FALSE, zero.appeal = 1, directed = FALSE,
#  	 	algorithm ="psumtree", start.graph = NULL)

load("Graph.RData")
g       <- DFNET_graph[[1]]
N.Nodes <- length(V(g)) 

n.samples=1000

node1 <- DFNET_graph$Selected_Module[1]
node2 <- DFNET_graph$Selected_Module[2]
node3 <- DFNET_graph$Selected_Module[3]
node4 <- DFNET_graph$Selected_Module[4]


# Generate input data 
MAT_a    <- matrix(sample(c(0,1), n.samples*N.Nodes, replace=TRUE), n.samples, N.Nodes)
MAT_b    <- matrix(sample(c(0,1), n.samples*N.Nodes, replace=TRUE), n.samples, N.Nodes)


# Define the response
#target <- as.numeric(xor(MAT[,node1], MAT[,node2]))
target <- as.numeric(xor(MAT_a[,node1] & MAT_b[,node2], MAT_a[,node3] & MAT_b[,node4]))
#target <- as.numeric(xor(xor(MAT[,node1], MAT[,node2]), MAT[,node3]))

IN_a <- cbind(MAT_a, target)
IN_b <- cbind(MAT_b, target)

IN_a   <- as.data.frame(IN_a)
colnames(IN_a) <- c(paste("AN_", 1:N.Nodes, sep=""),"target")
IN_b   <- as.data.frame(IN_b)
colnames(IN_b) <- c(paste("BN_", 1:N.Nodes, sep=""),"target")

IN <- list(IN_a, IN_b)

return(list(graph=g, Feature_Matrix=IN, Selected_Module=c(node1,node2,node3,node4)))

}