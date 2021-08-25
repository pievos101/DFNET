
#id       <- which.max(DFNET_object$DFNET_MODULES_AUC)
#vertices <- DFNET_object$DFNET_MODULES[[id]]
#library(igraph)

DFNET_get_module <- function(Nodes, DFNET_graph, DFNET_Eimp){

vertices <- Nodes 
#vertices <- as.numeric(strsplit(DFNET_mod[1,1]," ")[[1]])

module.genes <- DFNET_graph$gene.names[vertices]

g <- induced.subgraph(DFNET_graph[[1]], vertices)

vertex_attr(g) <- list(name = module.genes)
                        #color = rep("yellow", gorder(g)))

#plot(g, vertex.shape="none")


# Extract Edge Importance for the edges
EDGELIST <- as_edgelist(DFNET_graph$graph, names = TRUE)

ids <- apply(EDGELIST,1,function(x){all(is.element(x, vertices))})
ids <- which(ids)

MODULE_EDGELIST_SAVE  <- EDGELIST[ids,]
MODULE_EDGELIST       <- EDGELIST[ids,]
MODULE_EDGELIST[,1]   <- DFNET_graph$gene.names[MODULE_EDGELIST_SAVE[,1]]
MODULE_EDGELIST[,2]   <- DFNET_graph$gene.names[MODULE_EDGELIST_SAVE[,2]]

MODULE_EDGE_IMP       <- cbind(MODULE_EDGELIST,DFNET_Eimp[ids])

FINAL_MODULE_EDGE_IMP <- as_edgelist(g)

FINAL_EDGE_IMP <- rep(NaN,dim(FINAL_MODULE_EDGE_IMP)[1])

for(xx in 1:dim(FINAL_MODULE_EDGE_IMP)[1]){

	for(yy in 1:dim(MODULE_EDGE_IMP)[1]){

		if(all(is.element(FINAL_MODULE_EDGE_IMP[xx,1:2],MODULE_EDGE_IMP[yy,1:2]))){
			FINAL_EDGE_IMP[xx] <- MODULE_EDGE_IMP[yy, 3]
		}
	}
}

RES <- cbind(FINAL_MODULE_EDGE_IMP, FINAL_EDGE_IMP ) 

# In case edges are dublicated !?
RES <- unique(RES)
colnames(RES) <- c("GENE","GENE","EDGE_IMP") 
return(as.data.frame(RES))

}

# Plot the edge Importance

#RES_EDGE_IMP <- apply(RES,1, function(x){paste(x[1],"-",x[2])})
#RES_EDGE_IMP <- cbind(RES_EDGE_IMP,RES[,3])
#RES_EDGE_IMP <- as.data.frame(RES_EDGE_IMP)
#colnames(RES_EDGE_IMP) <- c("EDGE","IMP")
#RES_EDGE_IMP$IMP <- as.numeric(RES_EDGE_IMP$IMP)
#RES_EDGE_IMP <- RES_EDGE_IMP[order(RES_EDGE_IMP$IMP, decreasing=TRUE),]

#edge_imp        <- RES_EDGE_IMP$IMP
#names(edge_imp) <- RES_EDGE_IMP$EDGE

#barplot(sort(edge_imp), horiz = TRUE, las=1, xlab="Edge Importance", 
#	cex.names=0.35, col="cadetblue")

#library(ggplot2)
#p <- ggplot(RES_EDGE_IMP, aes(y=IMP, x=EDGE)) + 
#    geom_bar(position="dodge", stat="identity") +
#    ylab("Edge Importance") +
#    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))