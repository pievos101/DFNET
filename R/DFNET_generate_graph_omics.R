#' Annotate a graph with multi-omic data.
#'
#' @param network
#' The graph to annotate as (weighted) edgelist.
#' @param features
#' The features (omics) to annotate the graph with as matrix.
#' @param target
#' The target feature as vector.
#' @param cut.off
#' The weight threshold below which edges within \code{network} are considered
#' insignificant.
#' @return a list of shape \code{(graph, feature.matrix, feature.names)},
#' wherein \code{graph} and \code{feature.matrix} are reduced variants of
#' \code{network} and \code{features} respectively, s.t. only nodes within the
#' network that correspond to features and only features that correspond to
#' nodes are kept.
DFNET_generate_graph_omics <- function(network, features, target, cut.off = NaN) {
    # @FIXME
    # This procedure performs multiple cleaning steps, which should probably be
    # extracted into their own functions.
    NET <- network
    OMICS <- features
    GENE_NAMES <- lapply(OMICS, colnames)

    # @FIXME
    # Something is wrong with read.table ...
    # So, clean data again

    # Harmonize Multi-Omics data
    # common_genes    <- intersect(colnames(omic1), colnames(omic2))
    common_genes <- Reduce(intersect, GENE_NAMES)

    for (xx in 1:length(OMICS)) {
        OMICS[[xx]] <- OMICS[[xx]][, common_genes]
    }

    # omic1  <- omic1[,common_genes]
    # omic2  <- omic2[,common_genes]

    # HARMONIZE WITH PPI (1)
    ids1 <- match(NET[, 1], common_genes)
    ids2 <- match(NET[, 2], common_genes)

    na.ids1 <- which(is.na(ids1))
    na.ids2 <- which(is.na(ids2))

    na.IDS <- unique(c(na.ids1, na.ids2))

    if (length(na.IDS) > 0) {
        NET <- NET[-na.IDS, ]
    }

    # HARMONIZE WITH THE PPI (2)
    PPI_genes <- unique(c(as.character(NET[, 1]), as.character(NET[, 2])))
    ids <- match(common_genes, PPI_genes)
    na.ids <- which(is.na(ids))

    if (length(na.ids) > 0) {
        for (xx in 1:length(OMICS)) {
            OMICS[[xx]] <- OMICS[[xx]][, -na.ids]
        }
    }

    if (!is.na(cut.off)) {
        thres <- quantile(NET[, 3], probs = cut.off)
        ids <- which(NET[, 3] > thres)

        # Renew data accordingly
        NET <- NET[ids, ]
        PPI_genes <- unique(c(as.character(NET[, 1]), as.character(NET[, 2])))
        ids <- match(colnames(OMICS[[1]]), PPI_genes)
        na.ids <- which(is.na(ids))
        for (xx in 1:length(OMICS)) {
            OMICS[[xx]] <- OMICS[[xx]][, -na.ids]
        }
    }

    # concert NET data to numeric values

    ids1 <- match(NET[, 1], colnames(OMICS[[1]]))
    ids2 <- match(NET[, 2], colnames(OMICS[[1]]))

    NET[, 1] <- ids1
    NET[, 2] <- ids2

    N.Nodes <- dim(OMICS[[1]])[2]

    for (xx in 1:length(OMICS)) {
        OMICS[[xx]] <- cbind(OMICS[[xx]], t(target))
    }


    feature.names <- colnames(OMICS[[1]])

    for (xx in 1:length(OMICS)) {
        colnames(OMICS[[xx]]) <- c(paste(LETTERS[xx], "N_", 1:N.Nodes, sep = ""), "target")
    }


    IN <- OMICS

    g <- graph_from_edgelist(as.matrix(NET[, 1:2]), directed = FALSE)


    return(list(graph = g, feature.matrix = IN, feature.names = feature.names))
}

DFNET_generate_graph_Omics <- DFNET_generate_graph_omics
