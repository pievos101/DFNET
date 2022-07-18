## DFNET_calc_feature_importance - Calculate the importance of features.

## Copyright © 2021 Bastian Pfeifer <bastianxpfeifer@gmail.com>
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

#' Importance of individual features
#'
#' @description
#' Calculates the average importance of features over multiple decision
#' trees.
#'
#' @param forest The trained \code{DFNET.forest}.
#' @param features matrix or 3D array. The features on which it was trained.
#' @param sep string.
#' @return A matrix of importance scores with one row per column in
#' \code{features} and one column per matrix (in a 3D array).
feature_importance <- function(forest, features, sep = "$") {
    orig.dim <- dim(features)
    importance.dimnames <- list(
        dimnames(features)[[2]],
        if (length(dim(features) == 2)) NULL else dimnames(features)[[3]]
    )

    features <- flatten2ranger(features, sep = sep)
    feature.names <- dimnames(features)[[2]]
    trees <- forest$trees

    importance <- rep_len(0, length(feature.names))
    ## count <- rep_len(0, length(trees))
    for (tree in trees) {
        varimp <- tree$variable.importance
        varnam <- names(varimp)

        ids <- match(feature.names, varnam)
        sel <- !is.na(ids)
        ## count[sel] <- count[sel] + 1
        importance[sel] <- importance[sel] + varimp[ids[sel]]
    }
    importance <- importance / length(trees)

    return(matrix(
        importance, nrow = orig.dim[[2]],
        dimnames = importance.dimnames
    ))
}

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
