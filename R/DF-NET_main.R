# Peform the novel module selection algorithm
#library(ranger)
#library(pROC)

DFNET <- function(DFNET_graph, ntrees=100, niter=200, init.mtry=NaN){

g          <- DFNET_graph[[1]]
IN         <- DFNET_graph[[2]] #its a list # MM DATA

# MM -- START
MultiModal <- !is.data.frame(IN)
if(MultiModal){
    target  <- IN[[1]][,"target"] # MM DATA
}else{
    target  <- IN[,"target"] 
}
## MM -- END

Nodes   <- V(g)
N.Nodes <- length(Nodes)

SELECTED_NODES   <- list()
SELECTED_NODES_WEIGHTS <- list()

if(is.na(init.mtry)){
 n.nodes.per.tree <- ceiling(sqrt(N.Nodes))
}else{
 n.nodes.per.tree <- init.mtry   
}

N.trees          <- ntrees 

# Sample Start Nodes for Random Walk
Start.Nodes <- sample(Nodes, N.trees, replace=TRUE)

count <- 1
for (xx in Start.Nodes){

    SELECTED_NODES[[count]] <- as.numeric(random_walk(g, xx, n.nodes.per.tree))

    # Not optimal solution
    while(length(SELECTED_NODES[[count]]) < n.nodes.per.tree){
        print("retry ..")
        xx2 <- sample(Nodes, 1)
        SELECTED_NODES[[count]] <- as.numeric(random_walk(g, xx2, n.nodes.per.tree))
    }
  
    SELECTED_NODES_WEIGHTS[[count]] <- as.numeric(table(SELECTED_NODES[[count]]))
    
    count <- count + 1 
}

# Create the TREES
DECISION_TREES <- list()

for (xx in 1:length(SELECTED_NODES)){
  
    # feature weights are sorted
    UNIQUE_NODES <- sort(unique(SELECTED_NODES[[xx]])) 
    UNIQUE_NODES_WEIGHTS <- SELECTED_NODES_WEIGHTS[[xx]]
    WEIGHTS <- UNIQUE_NODES_WEIGHTS
    
    # START -- MM Data
    # Collect the features from the MM space
    if(MultiModal){

     MM_DATA  <- IN[[1]][,UNIQUE_NODES]
        #print("MM Data!")
     for(mm in 2:length(IN)){
      MM_DATA <- cbind(MM_DATA, IN[[mm]][,UNIQUE_NODES])    
      WEIGHTS <- c(WEIGHTS, UNIQUE_NODES_WEIGHTS)
     }

     MM_DATA <- cbind(MM_DATA, target)
     MM_DATA <- as.data.frame(MM_DATA)
    
    }else{
    
     MM_DATA   <- IN[,UNIQUE_NODES]
     MM_DATA   <- as.data.frame(cbind(MM_DATA, target))
    }
    # END -- MM Data

    #print(dim(MM_DATA))

    # Peform Feature Selection
    rf.sim  <-  ranger(dependent.variable.name = "target",
                    data=MM_DATA, # MM DATA     
                    #data=IN[,c(SELECTED_NODES[[xx]], N.Nodes+1)], 
                    split.select.weights=WEIGHTS / sum(WEIGHTS),
                    verbose = FALSE,
                    classification=TRUE, 
                    importance ="impurity", 
                    num.trees=1, 
                    mtry=dim(MM_DATA)[2] - 1) # at each split consider all variables
                    #replace = TRUE)

    DECISION_TREES[[xx]] <- rf.sim

}

# Calculate the Accuracy of the Graph Classifier
# Generate Outcome based on Majority Vote

PRED            <- sapply(DECISION_TREES,function(x){x$predictions})
AUC_PER_TREE    <- apply(PRED,2,function(x){auc(target,x, na.rm=TRUE, levels = c(0, 1), direction = "<")[1]})
range01         <- function(x){(x-min(x))/(max(x)-min(x))}
AUC_PER_TREE_01 <- range01(AUC_PER_TREE)


PRED_VEC   <- apply(PRED,1, function(x){
                sum(x==1, na.rm=TRUE)-sum(x==0, na.rm=TRUE)
                })

PRED_VEC[PRED_VEC==0] <- NaN
PRED_VEC[PRED_VEC<0]  <- 0
PRED_VEC[PRED_VEC>0]  <- 1

auc(target, PRED_VEC,  na.rm=TRUE, levels = c(0, 1), direction = "<")


##################################################################
# NOW ITERATE - GREEDY STEP
##################################################################

DECISION_TREES_ALL      <- list()
DECISION_TREES_ALL      <- DECISION_TREES
SELECTED_NODES_X        <- SELECTED_NODES
SELECTED_NODES_X_OLD    <- SELECTED_NODES_X
SELECTED_NODES_WEIGHTS  <- list()

WALK.DEPTH     <- rep(ceiling(n.nodes.per.tree), N.trees) 
WALK.DEPTH_OLD <- WALK.DEPTH
#iter  <- 20

range01             <- function(x){
    if(all(x==1)){return(x)}
    (x-min(x))/(max(x)-min(x))
}

AUC_PER_TREE_01_OLD <- AUC_PER_TREE_01
AUC_PER_TREE_OLD    <- AUC_PER_TREE 


#while ( any(WALK.DEPTH > 2) ){
ITER     <- niter
converge <- rep(FALSE, N.trees)

for(xx in 1:ITER){

cat(xx, " of ", ITER, "\n")

#if(ITER%%10==0){converge <- rep(FALSE, N.trees)}
#WALK.DEPTH_OLD <- WALK.DEPTH

# Sample
  #print(AUC_PER_TREE_01)
  ids_keep             <- sample(1:length(AUC_PER_TREE_OLD), prob = AUC_PER_TREE_OLD, replace = TRUE)

  SELECTED_NODES_X     <- SELECTED_NODES_X_OLD[ids_keep]

  WALK.DEPTH           <- WALK.DEPTH_OLD[ids_keep]

  #AUC_PER_TREE_01_OLD  <- AUC_PER_TREE_01_OLD[ids_keep] 
 
  #SELECTED_NODES_X_OLD <- SELECTED_NODES_X_OLD[ids_keep]

  converge <- converge[ids_keep]

  #print(AUC_PER_TREE_01)
  #print(SELECTED_NODES_X)
  #print(converge)
  #print(WALK.DEPTH)

# Sample Start Nodes for Random Walk
  Start.Nodes_X <- sapply(SELECTED_NODES_X,function(x){sample(x,1)}) 

count <- 1
for (xx in Start.Nodes_X){

    #print(xx)
    #if(converge[xx]){
        #SELECTED_NODES_X[[count]] <- SELECTED_NODES_X_OLD[[count]]
    #   count <- count + 1 
    #   next
    #}

    SELECTED_NODES_X[[count]] <- as.numeric(random_walk(g, xx, WALK.DEPTH[count]))
    SELECTED_NODES_WEIGHTS[[count]] <- as.numeric(table(SELECTED_NODES_X[[count]]))

    count <- count + 1 
}

# Create the TREES
DECISION_TREES <- list()

for (xx in 1:length(SELECTED_NODES_X)){

  # feature weights are sorted
  UNIQUE_NODES <- sort(unique(SELECTED_NODES_X[[xx]])) 
  UNIQUE_NODES_WEIGHTS <- SELECTED_NODES_WEIGHTS[[xx]]
  WEIGHTS <- UNIQUE_NODES_WEIGHTS 

    # if(converge[xx]){next}
    # Peform Feature Selection

    # START -- MM Data
    # Collect the features from the MM space
    if(MultiModal){

     MM_DATA  <- IN[[1]][,UNIQUE_NODES]
        #print("MM Data!")
     for(mm in 2:length(IN)){
      MM_DATA <- cbind(MM_DATA, IN[[mm]][,UNIQUE_NODES])    
      WEIGHTS <- c(WEIGHTS, UNIQUE_NODES_WEIGHTS)
     }

     MM_DATA <- cbind(MM_DATA, target)
     MM_DATA <- as.data.frame(MM_DATA)
    
    }else{
    
     MM_DATA   <- IN[,UNIQUE_NODES]
     MM_DATA   <- as.data.frame(cbind(MM_DATA, target))

    }
    # END -- MM Data
    rf.sim  <-  ranger(dependent.variable.name = "target",
                    data = MM_DATA, # MM DATA
                    #data = IN[,c(SELECTED_NODES_X[[xx]], N.Nodes+1)], 
                    split.select.weights=WEIGHTS / sum(WEIGHTS),
                    verbose = FALSE,
                    classification = TRUE, 
                    importance = "impurity", 
                    num.trees = 1, 
                    mtry = dim(MM_DATA)[2] - 1, # at each split consider all variables
                    replace = TRUE) # crucial parameter!

    DECISION_TREES[[xx]] <- rf.sim

}

PRED            <- sapply(DECISION_TREES,function(x){x$predictions})
AUC_PER_TREE    <- apply(PRED,2,function(x){auc(target,x, na.rm=TRUE, levels = c(0, 1), direction = "<")[1]})
#AUC_PER_TREE   <- range01(AUC_PER_TREE) #@FIXME

#print(AUC_PER_TREE)

ids_schrink     <- which(AUC_PER_TREE>=AUC_PER_TREE_OLD)
ids_schrink_not <- which(AUC_PER_TREE<AUC_PER_TREE_OLD)

# save the run

#AUC_PER_TREE_01_OLD <- AUC_PER_TREE_01
if(length(ids_schrink)>0){
 AUC_PER_TREE_OLD[ids_schrink]        <- AUC_PER_TREE[ids_schrink]
 SELECTED_NODES_X_OLD[ids_schrink]    <- SELECTED_NODES_X[ids_schrink] 
 WALK.DEPTH_OLD[ids_schrink]          <- WALK.DEPTH[ids_schrink] - 1
}
if(length(ids_schrink_not)>0){
 # ?
 AUC_PER_TREE_OLD[ids_schrink_not]     <- AUC_PER_TREE_OLD[ids_schrink_not]
 SELECTED_NODES_X_OLD[ids_schrink_not] <- SELECTED_NODES_X_OLD[ids_schrink_not] 
 WALK.DEPTH_OLD[ids_schrink_not]       <- WALK.DEPTH_OLD[ids_schrink_not]
 
 DECISION_TREES[ids_schrink_not]       <- NULL 
}

DECISION_TREES_ALL <- c(DECISION_TREES_ALL, DECISION_TREES)

#print(WALK.DEPTH_OLD)

#converge[ids_schrink]     <- FALSE 
converge[ids_schrink_not] <- TRUE 

WALK.DEPTH_OLD[WALK.DEPTH_OLD<2]  <- 2

#print(SELECTED_NODES_X_OLD)

#print(converge)

#if(all(converge==TRUE) & (ITER%%10==0)){
#    print("reset")
#   converge <- rep(FALSE, N.trees)
#}

#print(converge)

#if(all(WALK.DEPTH_OLD==WALK.DEPTH)){
#   print("EARLY STOP!")
#   break
#}

#if(all(converge)){
#   print("EARLY STOP!")
#   break
#}


#WALK.DEPTH_OLD <- WALK.DEPTH
#print(WALK.DEPTH)

}# End of for loop

return(list(DFNET_graph=DFNET_graph, DFNET_trees=DECISION_TREES_ALL, 
            DFNET_MODULES=SELECTED_NODES_X_OLD, DFNET_MODULES_AUC=AUC_PER_TREE_OLD))

}# End of function
