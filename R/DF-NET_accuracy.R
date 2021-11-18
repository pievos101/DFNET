# RFNET calc Accuracy
range01        <- function(x){(x-min(x))/(max(x)-min(x))}

DFNET_accuracy <- function(DFNET_graph, DFNET_object, n.last.trees=NaN){

g        <- DFNET_graph[[1]]
IN       <- DFNET_graph[[2]]

MultiModalData <- !is.data.frame(IN)

if(MultiModalData){
	target  <- IN[[1]][,"target"] # MM DATA
}else{
	target  <- IN[,"target"] 
}

EDGELIST <- as_edgelist(g, names = TRUE)

if(is.na(n.last.trees)){
	DECISION_TREES_ALL <- DFNET_object$DFNET_trees
}else{
	DECISION_TREES_ALL <- tail(DFNET_object$DFNET_trees, n.last.trees)
}

# Calculate Classifier prediction
PRED            <- sapply(DECISION_TREES_ALL,function(x){x$predictions})
AUC_ALL_TREES   <- apply(PRED,2,function(x){auc(target,x, na.rm=TRUE, levels = c(0, 1), direction = "<")[1]})

PRED_VEC_X   <- apply(PRED,1, function(x){
				sum(x==1, na.rm=TRUE)-sum(x==0, na.rm=TRUE)
				})

PRED_VEC_X[PRED_VEC_X==0] <- NaN
PRED_VEC_X[PRED_VEC_X<0]  <- 0
PRED_VEC_X[PRED_VEC_X>0]  <- 1


#print(auc(target, PRED_VEC,  na.rm=TRUE, levels = c(0, 1), direction = "<")[1])
accuracy = auc(target, PRED_VEC_X,  na.rm=TRUE, levels = c(0, 1), direction = "<")[1]

return(accuracy)

}