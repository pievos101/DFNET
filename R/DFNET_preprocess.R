
DFNET_preprocess <- function(DFNET_graph) {
    # XXX: Express in terms of relat
    range01 <- function(x) {
        res <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
        if (is.infinite(res[1]) || is.na(res[1])) {
            return(x)
        } else {
            return(res)
        }
    }


    FEAT <- DFNET_graph$Feature_Matrix

    for (xx in 1:length(FEAT)) {
        FEAT[[xx]] <- as.data.frame(apply(FEAT[[xx]], 2, range01))
        colnames(FEAT[[xx]]) <- colnames(DFNET_graph$Feature_Matrix[[xx]])
    }

    DFNET_graph$Feature_Matrix <- FEAT
    return(DFNET_graph)
}
