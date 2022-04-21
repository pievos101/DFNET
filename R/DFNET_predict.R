# predict outcome based on forest
DFNET_predict <- function(DFNET_object, DFNET_graph, n.last.trees=NaN, tree.ID=NaN){
  #require(pROC)

  # retrieve full data
  dataset <- do.call(cbind, DFNET_graph[[2]])
  # leave only one target column
  n.mod   <- length(DFNET_graph[[2]]) 
  dataset <- dataset[, -head(which(colnames(dataset) == "target"),n.mod-1)]

  pred <- data.frame()

  #
  if(is.na(n.last.trees) & is.na(tree.ID)){
    DECISION_TREES_ALL <- DFNET_object$DFNET_trees
  }

  #
  if(!is.na(n.last.trees) & is.na(tree.ID)){
    DECISION_TREES_ALL <- tail(DFNET_object$DFNET_trees, n.last.trees)
  }

  # 
  if(is.na(n.last.trees) & !is.na(tree.ID)){
    DECISION_TREES_ALL <- DFNET_object$DFNET_trees[tree.ID]
  }
  
  
  # Predict 
  pred  <- matrix(NaN, length(DECISION_TREES_ALL), dim(dataset)[1])
  count <- 1
  for (tree in DECISION_TREES_ALL) {
    #pred <- rbind(pred, predict(tree, dataset)$predictions)
    if (count %% 100 == 0){
      cat(count, " of ", length(DECISION_TREES_ALL) , " trees \n")
    }
    pred[count,] <- predict(tree, dataset)$predictions
    count <- count + 1
  }
  
  val <- apply(pred, 2, function(x){
    return(as.numeric(names(sort(table(x),decreasing=TRUE))[1]))
    })
  
  names(val) <- rownames(dataset)

return(val)  
}