# Calculate Module Importance

DFNET_modules <- function(DFNET_graph, DFNET_object, DFNET_eImp){
	
N.trees <- length(DFNET_object$DFNET_MODULES)
N.Nodes <- length(V(DFNET_graph$graph))

N.ALL.TREES <- length(DFNET_object$DFNET_trees)

# TREE IDS
TREE_IDS    <- (N.ALL.TREES-N.trees+1):N.ALL.TREES

SELECTED_NODES_X_OLD <- DFNET_object$DFNET_MODULES
AUC_PER_TREE_OLD     <- DFNET_object$DFNET_MODULES_AUC
EDGE_IMP_FINAL       <- DFNET_eImp


EDGELIST <- as_edgelist(DFNET_graph$graph, names = TRUE)


MODULE_MATRIX <- matrix(0, N.trees, N.Nodes)

for (xx in 1:length(SELECTED_NODES_X_OLD)){
    node.ids <- SELECTED_NODES_X_OLD[[xx]]
	 MODULE_MATRIX[xx, node.ids] <- 1
}

# TREE IDS
rownames(MODULE_MATRIX) <- TREE_IDS

MODULE_MATRIX_unique <- unique(MODULE_MATRIX)

# update TREE IDS
TREE_IDS1     <- as.numeric(rownames(MODULE_MATRIX_unique))

BEST_MODULES <- vector("list", dim(MODULE_MATRIX_unique)[1])
names(BEST_MODULES) <- TREE_IDS1

for(xx in 1:dim(MODULE_MATRIX_unique)[1]){
	BEST_MODULES[[xx]] <- which(MODULE_MATRIX_unique[xx,]==1)
}

Module_IMP <- rep(0, length(BEST_MODULES))

for( xx in 1:length(BEST_MODULES)){

	module <- BEST_MODULES[[xx]]
 
 for(yy in 1:dim(EDGELIST)[1]){
    
    edge <- EDGELIST[yy,]
    ids  <- match(edge, module)

     if(all(!is.na(ids))){
     	Module_IMP[xx] <- Module_IMP[xx] + EDGE_IMP_FINAL[yy]
     }

 }  
}

Module_length    <- sapply(BEST_MODULES,length)
BEST_MODULES_IMP <- Module_IMP/Module_length 

names(BEST_MODULES_IMP) <- 1:length(BEST_MODULES_IMP)
BEST_MODULES_IMP_SORTED <- sort(BEST_MODULES_IMP, decreasing=TRUE)
ids <- as.numeric(names(BEST_MODULES_IMP_SORTED))

BEST_MODULES_SORTED     <- BEST_MODULES[ids]

#BEST_MODULES_SORTED_AUC <- AUC_PER_TREE_OLD[ids]

BEST_MODULES_SORTED

# update TREE IDS
TREE_IDS2 <- names(BEST_MODULES_SORTED)

BEST_MODULES_IMP_SORTED

# Get the AUC of these unique modules
BEST_MODULES_AUC_IMP_SORTED <- vector("list",length(BEST_MODULES_IMP_SORTED))
# TREE IDS
names(BEST_MODULES_AUC_IMP_SORTED) <- TREE_IDS2

for(xx in 1:length(BEST_MODULES_IMP_SORTED)){

 m1 <- BEST_MODULES_SORTED[[xx]]
		
	for(yy in 1:length(SELECTED_NODES_X_OLD)){
 
     m2 <- SELECTED_NODES_X_OLD[[yy]]

     #print(m1)
     #print(m2)

     if(all(m1==unique(sort(m2)))){
     	BEST_MODULES_AUC_IMP_SORTED[[xx]] <- c(BEST_MODULES_AUC_IMP_SORTED[[xx]], AUC_PER_TREE_OLD[yy]) 
     }

	}

}

BEST_MODULES_AUC_IMP_SORTED <- sapply(BEST_MODULES_AUC_IMP_SORTED, mean)

BEST_MODULES_SORTED_STRING <- sapply(BEST_MODULES_SORTED, function(x){
	paste(x,collapse=" ")
	})


RES <- data.frame(BEST_MODULES_SORTED_STRING,BEST_MODULES_IMP_SORTED, BEST_MODULES_AUC_IMP_SORTED)
colnames(RES) <- c("Module","#EDGE_IMP","AUC")
#rownames(RES) <- NULL


IMP_total <- apply(RES[,c(2,3)],1,sum)
RES2      <- cbind(RES, IMP_total)
colnames(RES2) <- c("Module","EDGE_IMP","AUC","IMP")
RES2      <- RES2[order(RES2$IMP, decreasing=TRUE),]


# Fix an Issue
# Order by AUC!
#RES2  <- RES[order(RES$AUC, decreasing=TRUE),]
#ids   <- which(RES2$AUC==RES2$AUC[1])
#if(length(ids)>=2){
#TEMP <- RES2[ids,]
#TEMP <- TEMP[order(TEMP[,2], decreasing=TRUE),]
#RES2[ids,] <- TEMP
#}

return(RES2)

} 