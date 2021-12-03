DFNET_explain <- function(forest, treeID=NaN, plots=FALSE){
	
	require(treeshap)
	# Unify 

	# TreeSHAP


}

unify_forest <- function(trees, data) {
  #:# create, fix, and concatenate unified forests (single trees)
  unified_model <- data.frame()
  i <- 0
  support <- 0
  for (tree in trees) {
    cat(i + 1, " of ", length(trees), "\n")
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