DFNET_explain <- function(DFNET_object, DFNET_graph, n.last.trees=NaN, tree.ID=NaN){
	
	#DFNET_graph contains the data
	require(treeshap)
	# retrieve full data
	# retrieve full data
    dataset <- do.call(cbind, DFNET_graph[[2]])
    # leave only one target column
    n.mod   <- length(DFNET_graph[[2]]) 
    dataset <- dataset[, -head(which(colnames(dataset) == "target"),n.mod-1)]


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
  
	# Unify 
	forest <- unify_forest(DECISION_TREES_ALL, dataset)

	# calculate treeshap
	forest_shap <- treeshap(forest, dataset[, -dim(dataset)[2]])
	#sv <- forest_shap$shaps
	return(forest_shap)
}

unify_forest <- function(trees, data) {
  #:# create, fix, and concatenate unified forests (single trees)
  unified_model <- data.frame()
  i <- 0
  support <- 0
  for (tree in trees) {
    cat(i + 1, " of ", length(trees), " trees unified \n")
    unified_tree <- ranger.unify(tree, data)
    unified_tree$model$Tree <- i
    unified_tree$model$Yes <- unified_tree$model$Yes + support
    unified_tree$model$No <- unified_tree$model$No + support
    unified_tree$model$Missing <- unified_tree$model$Missing + support
    unified_model <- rbind(unified_model, unified_tree$model)
    i <- i + 1
    support <- support + nrow(unified_tree$model)
  }
  unified_model$Prediction <- unified_model$Prediction / length(trees)
  unified_forest <- unified_tree
  unified_forest$model <- unified_model
  unified_forest$data <- data
  unified_forest
}