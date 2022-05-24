DFNET_generate_graph_omics <- function(PPI, FEATURES, TARGET, cut.off = NaN) {


    # omic1 <- FEATURES[[1]]
    # omic2 <- FEATURES[[2]]

    OMICS <- FEATURES
    GENE_NAMES <- lapply(OMICS, colnames)

    NET <- PPI

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

    # omic1       <- omic1[,-na.ids]
    # omic2       <- omic2[,-na.ids]

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
        # omic1  <- omic1[,-na.ids]
        # omic2  <- omic2[,-na.ids]
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

    # omic1           <- cbind(omic1, t(TARGET))
    # omic2           <- cbind(omic2, t(TARGET))

    for (xx in 1:length(OMICS)) {
        OMICS[[xx]] <- cbind(OMICS[[xx]], t(TARGET))
    }


    gene.names <- colnames(OMICS[[1]])

    # colnames(omic1) <- c(paste("AN_", 1:N.Nodes, sep=""),"target")
    # omic1           <- as.data.frame(omic1)
    # colnames(omic2) <- c(paste("BN_", 1:N.Nodes, sep=""),"target")
    # omic2           <- as.data.frame(omic2)

    for (xx in 1:length(OMICS)) {
        colnames(OMICS[[xx]]) <- c(paste(LETTERS[xx], "N_", 1:N.Nodes, sep = ""), "target")
    }


    IN <- OMICS # list(omic1, omic2)

    # remove NaNs from PPI --> no longer needed
    # ids1 <- which(is.na(NET[,1]))
    # ids2 <- which(is.na(NET[,2]))
    # NET  <- NET[-c(ids1,ids2),]

    g <- graph_from_edgelist(as.matrix(NET[, 1:2]), directed = FALSE)


    return(list(graph = g, Feature_Matrix = IN, gene.names = gene.names))
}

DFNET_generate_graph_Omics <- DFNET_generate_graph_omics
