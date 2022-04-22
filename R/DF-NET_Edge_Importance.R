# DFNET - Edge Importances
range01        <- function(x){(x-min(x))/(max(x)-min(x))}

DFNET_Edge_Importance <- function(DFNET_graph, DFNET_object, parallel=FALSE){

g        <- DFNET_graph[[1]]
IN       <- DFNET_graph[[2]]

MultiModalData <- !is.data.frame(IN)

if(MultiModalData){
	target  <- IN[[1]][,"target"] # MM DATA
}else{
	target  <- IN[,"target"] 
}


EDGELIST <- as_edgelist(g, names = TRUE)


DECISION_TREES_ALL <- DFNET_object$DFNET_trees

# Get the variables involved in each TREE
X_TREE_VARS <- lapply(DECISION_TREES_ALL, 
	function(x){names(x$variable.importance)})

TREE_VARS <- lapply(X_TREE_VARS, function(x){
				as.numeric(substring(x,4))})

# Get the out-of-bag performance of each TREE
TREE_AUCS <- sapply(DECISION_TREES_ALL, function(x){
				auc(target, x$predictions, 
				na.rm=TRUE, levels = c(0, 1), direction = "<")[1]})


# Calculate EDGE IMPORTANCE
EDGE_IMP <- numeric(dim(EDGELIST)[1])
EDGELIST_LIST = lapply(apply(EDGELIST,1,list),unlist)

if(parallel){
require(parallel)
for (xx in 1:length(TREE_VARS)){

  if (xx %% 100 == 0) cat(xx, " of ", length(TREE_VARS), " trees \n")

  vars <- TREE_VARS[[xx]]
  res  <- mclapply(EDGELIST_LIST,
  	function(x){all(is.element(x,vars))}, mc.cores=detectCores()-2)
  res <- unlist(res)

  EDGE_IMP[res] <- EDGE_IMP[res] + TREE_AUCS[xx] 
}
} # if parallel

if(!parallel){

# VERSION 1 ###################
for (xx in 1:length(TREE_VARS)){

  if (xx %% 100 == 0) cat(xx, " of ", length(TREE_VARS), " trees \n")

  vars <- TREE_VARS[[xx]]
  res  <- apply(EDGELIST,1,function(x){all(is.element(x,vars))})
  EDGE_IMP[res] <- EDGE_IMP[res] + TREE_AUCS[xx] 
}
################################

# VERSION 2 #####################
#for (xx in 1:length(EDGELIST)){

#	if (xx %% 1000 == 0) cat(xx, " of ", length(EDGELIST), " edges \n")
#	edge <- EDGELIST[[xx]]

#	res <- sapply(TREE_VARS, function(x){all(is.element(edge,x))})
#	EDGE_IMP[xx] <- sum(TREE_AUCS[res])
#}
#################################

}# if not parallel


EDGE_IMP_FINAL <- range01(EDGE_IMP)

return(EDGE_IMP_FINAL)
}


DFNET_Edge_Importance_old <- function(DFNET_graph, DFNET_object){

g        <- DFNET_graph[[1]]
IN       <- DFNET_graph[[2]]

MultiModalData <- !is.data.frame(IN)

if(MultiModalData){
	target  <- IN[[1]][,"target"] # MM DATA
}else{
	target  <- IN[,"target"] 
}


EDGELIST <- as_edgelist(g, names = TRUE)


DECISION_TREES_ALL <- DFNET_object$DFNET_trees

# Calculate EDGE IMPORTANCE
EDGE_IMP <- vector("list", dim(EDGELIST)[1])

for (xx in 1:dim(EDGELIST)[1]){

 if (xx %% 1000 == 0) cat(xx, " of ", dim(EDGELIST)[1], " edges \n")

 edge <- paste("N_", EDGELIST[xx,], sep="")

	for(yy in 1:length(DECISION_TREES_ALL)){

	     	Vimp <- names(DECISION_TREES_ALL[[yy]]$variable.importance)
	     	#print(Vimp)

	     	if(MultiModalData){
	     	 Vimp_names <- substring(Vimp,2)	
	     	}else{
	         Vimp_names <- Vimp    
	        }
	     	
	     	#print(edge)
	     	#print(Vimp_names)
	     	check <- match(edge,Vimp_names)
			
		if(all(!is.na(check))){
			EDGE_IMP[[xx]] <- c(EDGE_IMP[[xx]], auc(target, DECISION_TREES_ALL[[yy]]$predictions, na.rm=TRUE, levels = c(0, 1), direction = "<")[1])
		}
	}
} 

EDGE_IMP_FINAL <- sapply(EDGE_IMP, sum)
EDGE_IMP_FINAL <- range01(EDGE_IMP_FINAL)


return(EDGE_IMP_FINAL)


}