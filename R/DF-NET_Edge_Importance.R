# RFNET calc Accuracy
range01        <- function(x){(x-min(x))/(max(x)-min(x))}

DFNET_Edge_Importance <- function(DFNET_graph, DFNET_object){

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