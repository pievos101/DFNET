#DFNET_calc_feature_importance 

vert1 <- paste("AN_", vertices,sep="")
vert2 <- paste("BN_", vertices,sep="")

SUMIMP1 <- numeric(length(vertices))
COUNT1  <- numeric(length(vertices))

SUMIMP2 <- numeric(length(vertices))
COUNT2  <- numeric(length(vertices))

for(xx in 1:length(DFNET_object$DFNET_trees)){

	VARIMP <- DFNET_object$DFNET_trees[[xx]]$variable.importance
	NN	   <- names(VARIMP)
	
	ids    <- match(vert1, NN)

	if(!all(is.na(ids))){
	COUNT1[!is.na(ids)]  <- COUNT1[!is.na(ids)] + 1
	SUMIMP1[!is.na(ids)] <- SUMIMP1[!is.na(ids)] + VARIMP[ids[!is.na(ids)]] 
	}

	ids    <- match(vert2, NN)

	if(!all(is.na(ids))){
	COUNT2[!is.na(ids)]  <- COUNT2[!is.na(ids)] + 1
	SUMIMP2[!is.na(ids)] <- SUMIMP2[!is.na(ids)] + VARIMP[ids[!is.na(ids)]] 
	}

}

FEATURE1_imp <- SUMIMP1/COUNT1
FEATURE2_imp <- SUMIMP2/COUNT2

names(FEATURE1_imp) <- DFNET_graph$gene.names[vertices]
names(FEATURE2_imp) <- DFNET_graph$gene.names[vertices]

## GGPLOT
library(ggplot2)
library(reshape)

RES1 <- cbind(names(FEATURE1_imp),FEATURE1_imp)
RES2 <- cbind(names(FEATURE2_imp),FEATURE2_imp)
RES1 <- cbind(RES1,"mRNA")
RES2 <- cbind(RES2,"Methylation")

RES <- rbind(RES1,RES2)
rownames(RES) <- NULL
colnames(RES) <- c("Gene","IMP","Type")
RES <- as.data.frame(RES)
RES$IMP <- as.numeric(RES$IMP)

p <- ggplot(RES, aes(fill=Type, y=IMP, x=Gene)) + 
    geom_bar(position="dodge", stat="identity") +
    ylab("Impurity Importance") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot(p)    