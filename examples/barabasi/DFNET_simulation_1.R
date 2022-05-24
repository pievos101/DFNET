#
# Random Forests on Graphs
#
library(DFNET)
source("DFNET_generate_graph_sm.R")
source("DFNET_generate_graph_mm.R")
# source("DFNET_accuracy.R")

N.Nodes <- 20

DFNET_generate_graph_1 <- DFNET_generate_graph_SM # DFNET_generate_graph_MM
DFNET_graph <- NULL
while (length(DFNET_graph) == 0) {
    DFNET_graph <- DFNET_generate_graph_1(N.Nodes = N.Nodes)
}

DFNET_object <- DFNET(DFNET_graph, ntrees = 100, niter = 10)

# DFNET_acc      <- DFNET_accuracy(DFNET_graph, DFNET_object)

DFNET_Eimp <- DFNET_Edge_Importance(DFNET_graph, DFNET_object)

DFNET_mod <- DFNET_modules(DFNET_graph, DFNET_object, DFNET_Eimp)

DFNET_mod
DFNET_graph$Selected_Module

# PLOTS #############################################################

EDGELIST <- as_edgelist(DFNET_graph$graph, names = TRUE)
COLOR <- rep("cadetblue", dim(EDGELIST)[1])

for (xx in 1:dim(EDGELIST)[1]) {
    edge <- EDGELIST[xx, ]
    ids <- match(edge, as.numeric(strsplit(DFNET_mod[1, 1], " ")[[1]]))

    # print(ids)

    if (!any(is.na(ids))) {
        COLOR[xx] <- "coral2"
    }
}

plot(DFNET_graph$graph,
    edge.width = exp(DFNET_Eimp + 1), vertex.shape = "none",
    vertex.label.color = "black", edge.color = COLOR
)
