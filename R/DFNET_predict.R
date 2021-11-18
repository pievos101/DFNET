# predict outcome based on forest
DFNET_predict <- function(DFNET_object, DFNET_graph){
  require(pROC)

  # retrieve full data
  dataset <- do.call(cbind, DFNET_graph[[2]])
  # leave only one target column
  dataset <- dataset[, -which(colnames(dataset) == "target")[-1]]

  pred <- data.frame()
  for (tree in DFNET_object$DFNET_trees) {
    pred <- rbind(pred, predict(tree, dataset)$predictions)
  }
  
  val <- apply(pred, 2, function(x){
    return(as.numeric(names(sort(table(x),decreasing=TRUE))[1]))
    })
  
  names(val) <- rownames(dataset)

return(val)  
}