## DFNET_modules - Sort modules by importance.

## Copyright © 2021 Bastian Pfeifer <bastianxpfeifer@gmail.com>

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

module_importance <- function(graph, modules, edge_importances,
                              tree_importances, mc.cores = 1) {
    module_imp <- cbind(0, tree_importances)
    edges <- igraph::as_edgelist(graph, names = FALSE)

    for (e in 1:nrow(edges)) {
        edge <- edges[e, 1:2]
        pred <- function(module) {
            all(is.element(edge, module))
        }
        res <- unlist(mclapply(modules, pred, mc.cores = mc.cores))
        # XXX: What about edge weight?
        module_imp[res, 1] <- module_imp[res, 1] + edge_importances[e]
    }

    # FIXME: paper actually claims normalization by edge count, not module size
    module_imp[, 1] <- module_imp[, 1] / sapply(modules, length)
    module_imp[, 2] <- module_imp[, 1] + module_imp[, 2]
    colnames(module_imp) <- c("edge", "total")

    return(module_imp)
}

unique_module_importance <- function(graph, modules, edge_importances,
                                     tree_importances, mc.cores = 1,
                                     collapse = mean) {
    module_imp <- module_importance(
        graph, modules,
        edge_importances, tree_importances,
        mc.cores = mc.cores
    )

    modules <- lapply(modules, function(m) unique(sort(m)))
    module_name <- function(mod) paste(mod, collapse = " ")
    module_names <- sapply(modules, module_name)
    by_module_name <- order(module_names)

    module_imp <- module_imp[by_module_name, ]
    tree_importances <- tree_importances[by_module_name]
    modules <- modules[by_module_name]

    unique_modules <- unique(modules)
    unique_module_imp <- matrix(0, nrow = length(unique_modules), ncol = 3)
    colnames(unique_module_imp) <- c("tree", "edge", "total")
    rownames(unique_module_imp) <- sapply(unique_modules, module_name)

    for (module in unique_modules) {
        select <- sapply(modules, function(m) identical(m, module))

        unique_module_imp[module_name(module), 1:2] <- c(
            collapse(tree_importances[select]),
            collapse(module_imp[select, 1])
        )
    }
    unique_module_imp[, 3] <- unique_module_imp[, 1] + unique_module_imp[, 2]

    return(list(table = unique_module_imp, modules = unique_modules))
}

DFNET_modules <- function(DFNET_graph, DFNET_object, DFNET_eImp) {
    N.trees <- length(DFNET_object$DFNET_MODULES)
    N.Nodes <- length(V(DFNET_graph$graph))

    N.ALL.TREES <- length(DFNET_object$DFNET_trees)

    # TREE IDS
    TREE_IDS <- (N.ALL.TREES - N.trees + 1):N.ALL.TREES

    SELECTED_NODES_X_OLD <- DFNET_object$DFNET_MODULES
    AUC_PER_TREE_OLD <- DFNET_object$DFNET_MODULES_AUC
    EDGE_IMP_FINAL <- DFNET_eImp


    EDGELIST <- as_edgelist(DFNET_graph$graph, names = TRUE)


    MODULE_MATRIX <- matrix(0, N.trees, N.Nodes)

    for (xx in 1:length(SELECTED_NODES_X_OLD)) {
        node.ids <- SELECTED_NODES_X_OLD[[xx]]
        MODULE_MATRIX[xx, node.ids] <- 1
    }

    # TREE IDS
    rownames(MODULE_MATRIX) <- TREE_IDS

    MODULE_MATRIX_unique <- unique(MODULE_MATRIX)

    # update TREE IDS
    TREE_IDS1 <- as.numeric(rownames(MODULE_MATRIX_unique))

    BEST_MODULES <- vector("list", dim(MODULE_MATRIX_unique)[1])
    names(BEST_MODULES) <- TREE_IDS1

    for (xx in 1:dim(MODULE_MATRIX_unique)[1]) {
        BEST_MODULES[[xx]] <- which(MODULE_MATRIX_unique[xx, ] == 1)
    }

    Module_IMP <- rep(0, length(BEST_MODULES))

    for (xx in 1:length(BEST_MODULES)) {
        module <- BEST_MODULES[[xx]]

        for (yy in 1:dim(EDGELIST)[1]) {
            edge <- EDGELIST[yy, ]
            ids <- match(edge, module)

            if (all(!is.na(ids))) {
                Module_IMP[xx] <- Module_IMP[xx] + EDGE_IMP_FINAL[yy]
            }
        }
    }

    Module_length <- sapply(BEST_MODULES, length)
    BEST_MODULES_IMP <- Module_IMP / Module_length

    names(BEST_MODULES_IMP) <- 1:length(BEST_MODULES_IMP)
    BEST_MODULES_IMP_SORTED <- sort(BEST_MODULES_IMP, decreasing = TRUE)
    ids <- as.numeric(names(BEST_MODULES_IMP_SORTED))

    BEST_MODULES_SORTED <- BEST_MODULES[ids]

    # BEST_MODULES_SORTED_AUC <- AUC_PER_TREE_OLD[ids]

    BEST_MODULES_SORTED

    # update TREE IDS
    TREE_IDS2 <- names(BEST_MODULES_SORTED)

    BEST_MODULES_IMP_SORTED

    # Get the AUC of these unique modules
    BEST_MODULES_AUC_IMP_SORTED <- vector("list", length(BEST_MODULES_IMP_SORTED))
    # TREE IDS
    names(BEST_MODULES_AUC_IMP_SORTED) <- TREE_IDS2

    for (xx in 1:length(BEST_MODULES_IMP_SORTED)) {
        m1 <- BEST_MODULES_SORTED[[xx]]

        for (yy in 1:length(SELECTED_NODES_X_OLD)) {
            m2 <- SELECTED_NODES_X_OLD[[yy]]

            # print(m1)
            # print(m2)

            if (all(m1 == unique(sort(m2)))) {
                BEST_MODULES_AUC_IMP_SORTED[[xx]] <- c(BEST_MODULES_AUC_IMP_SORTED[[xx]], AUC_PER_TREE_OLD[yy])
            }
        }
    }

    BEST_MODULES_AUC_IMP_SORTED <- sapply(BEST_MODULES_AUC_IMP_SORTED, mean)

    BEST_MODULES_SORTED_STRING <- sapply(BEST_MODULES_SORTED, function(x) {
        paste(x, collapse = " ")
    })


    RES <- data.frame(BEST_MODULES_SORTED_STRING, BEST_MODULES_IMP_SORTED, BEST_MODULES_AUC_IMP_SORTED)
    colnames(RES) <- c("Module", "#EDGE_IMP", "AUC")
    # rownames(RES) <- NULL


    IMP_total <- apply(RES[, c(2, 3)], 1, sum)
    RES2 <- cbind(RES, IMP_total)
    colnames(RES2) <- c("Module", "EDGE_IMP", "AUC", "IMP")
    RES2 <- RES2[order(RES2$IMP, decreasing = TRUE), ]


    # Fix an Issue
    # Order by AUC!
    # RES2  <- RES[order(RES$AUC, decreasing=TRUE),]
    # ids   <- which(RES2$AUC==RES2$AUC[1])
    # if(length(ids)>=2){
    # TEMP <- RES2[ids,]
    # TEMP <- TEMP[order(TEMP[,2], decreasing=TRUE),]
    # RES2[ids,] <- TEMP
    # }

    return(RES2)
}
