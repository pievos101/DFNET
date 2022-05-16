#' Calculate the feature importances of particular vertices.
#'
#' @param vertices the vertices
#' @param detector the module detector
#' @param graph the graph on which \code{detector} was trained
#' @return a matrix with one row per gene in \code{graph} and one column
#' per vertex in \code{vertices}.
#' @examples
#' \dontrun{
#' graph <- get_some_graph()
#' detector <- DFNET(graph)
#' eimp <- DFNET_edge_importance(graph, detector)
#' candidates <- DFNET_modules(graph, detector, eimp)
#' top_module_vertices <- as.numeric(strsplit(candidates[1, 1], " ")[[1]])
#' fimp <- DFNET_calc_feature_importance(top_module_vertices, detector, graph)
#' }
DFNET_calc_feature_importance <- function(vertices, detector, graph) {
    VERT <- vector("list", length(graph$feature.matrix))
    for (xx in 1:length(graph$feature.matrix)) {
        VERT[[xx]] <- paste(LETTERS[xx], "N_", vertices, sep = "")
    }

    SUMIMP <- vector("list", length(graph$feature.matrix))
    COUNT <- vector("list", length(graph$feature.matrix))

    for (xx in 1:length(SUMIMP)) {
        SUMIMP[[xx]] <- numeric(length(vertices))
        COUNT[[xx]] <- numeric(length(vertices))
    }


    for (xx in 1:length(detector$DFNET_trees)) {
        VARIMP <- detector$DFNET_trees[[xx]]$variable.importance
        NN <- names(VARIMP)

        for (yy in 1:length(SUMIMP)) {
            ids <- match(VERT[[yy]], NN)
            if (!all(is.na(ids))) {
                COUNT[[yy]][!is.na(ids)] <- COUNT[[yy]][!is.na(ids)] + 1
                SUMIMP[[yy]][!is.na(ids)] <- SUMIMP[[yy]][!is.na(ids)] + VARIMP[ids[!is.na(ids)]]
            }
        }

        # ids    <- match(vert1, NN)

        # if(!all(is.na(ids))){
        # COUNT1[!is.na(ids)]  <- COUNT1[!is.na(ids)] + 1
        # SUMIMP1[!is.na(ids)] <- SUMIMP1[!is.na(ids)] + VARIMP[ids[!is.na(ids)]]
        # }

        # ids    <- match(vert2, NN)

        # if(!all(is.na(ids))){
        # COUNT2[!is.na(ids)]  <- COUNT2[!is.na(ids)] + 1
        # SUMIMP2[!is.na(ids)] <- SUMIMP2[!is.na(ids)] + VARIMP[ids[!is.na(ids)]]
        # }
    }

    # Normalize
    FEATURE_imp <- vector("list", length(SUMIMP))
    for (xx in 1:length(SUMIMP)) {
        FEATURE_imp[[xx]] <- SUMIMP[[xx]] / length(detector$DFNET_trees) # COUNT[[xx]]
        names(FEATURE_imp[[xx]]) <- graph$feature.names[vertices]
    }

    # FEATURE1_imp <- SUMIMP1/COUNT1
    # FEATURE2_imp <- SUMIMP2/COUNT2

    # names(FEATURE1_imp) <- graph$feature.names[vertices]
    # names(FEATURE2_imp) <- graph$feature.names[vertices]


    RES <- Reduce(rbind, FEATURE_imp)
    rownames(RES) <- paste("omic", 1:dim(RES)[1], sep = "")

    return(RES)
}

## GGPLOT
# library(ggplot2)
# library(reshape)

# RES1 <- cbind(names(FEATURE1_imp),FEATURE1_imp)
# RES2 <- cbind(names(FEATURE2_imp),FEATURE2_imp)
# RES1 <- cbind(RES1,"mRNA")
# RES2 <- cbind(RES2,"Methylation")

# RES <- rbind(RES1,RES2)
# rownames(RES) <- NULL
# colnames(RES) <- c("Gene","IMP","Type")
# RES <- as.data.frame(RES)
# RES$IMP <- as.numeric(RES$IMP)

# p <- ggplot(RES, aes(fill=Type, y=IMP, x=Gene)) +
#    geom_bar(position="dodge", stat="identity") +
#    ylab("Impurity Importance") +
#    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# plot(p)

# }
