# predict outcome based on forest
DFNET_predict <- function(DFNET_object, DFNET_graph, n.last.trees=NaN){
  #require(pROC)

  # retrieve full data
  dataset <- do.call(cbind, DFNET_graph_test[[2]])
  # leave only one target column
  dataset <- dataset[, -which(colnames(dataset) == "target")[-1]]

  pred <- data.frame()

  if(is.na(n.last.trees)){
    DECISION_TREES_ALL <- DFNET_object$DFNET_trees
  }else{
    DECISION_TREES_ALL <- tail(DFNET_object$DFNET_trees, n.last.trees)
  }

  for (tree in DECISION_TREES_ALL) {
    pred <- rbind(pred, predict(tree, dataset)$predictions)
  }
  
  val <- apply(pred, 2, function(x){
    return(as.numeric(names(sort(table(x),decreasing=TRUE))[1]))
    })
  
  names(val) <- rownames(dataset)

return(val)  
}