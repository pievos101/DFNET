# 
# Random Forests on Graphs
#
library(DFNET)
#source("~/GitLab/rf-feature-selection-on-graphs/DF-NET_main.R")
source("~/GitLab/rf-feature-selection-on-graphs/DF-NET_generate_Graph_MM.R")
#source("~/GitLab/rf-feature-selection-on-graphs/DF-NET_accuracy.R")
#source("~/GitLab/rf-feature-selection-on-graphs/DF-NET_modules.R")
#source("~/GitLab/rf-feature-selection-on-graphs/DF-NET_Edge_Importance.R")
source("~/GitLab/rf-feature-selection-on-graphs/DF-NET_load_Graph_MM.R")


N.Nodes        <- 50	#25 #50 #75 #100 

DFNET_graph    <- NULL
while(length(DFNET_graph)==0){
 DFNET_graph    <- DFNET_generate_graph_MM(N.Nodes=N.Nodes)
}

#load("Graph.RData")

#DFNET_graph    <- DFNET_load_graph_MM()


LOOPS          <- 50
#NITER          <- c(2, 5, 10, 20, 50, 100, 200, 300)
NITER          <- c(100, 200, 300)


COVERAGE       <- matrix(NaN, LOOPS, length(NITER))
AUC            <- matrix(NaN, LOOPS, length(NITER))
N_MODULES      <- matrix(NaN, LOOPS, length(NITER))

for (xx in 1:length(NITER)){ 

cat(xx, " of ", length(NITER), "done! ---------------- \n")

 for(yy in 1:LOOPS){

  #
  DFNET_graph    <- NULL
  while(length(DFNET_graph)==0){
  DFNET_graph    <- DFNET_generate_graph_MM(N.Nodes=N.Nodes)
  }
  #

  DFNET_object   <- DFNET(DFNET_graph, ntrees=100, niter=NITER[xx])

  #DFNET_acc      <- DFNET_accuracy(DFNET_graph, DFNET_object)

  DFNET_Eimp     <- DFNET_Edge_Importance(DFNET_graph, DFNET_object)

  DFNET_mod      <- DFNET_modules(DFNET_graph, DFNET_object, DFNET_Eimp)
  DFNET_mod

  bm             <- as.numeric(strsplit(DFNET_mod[1,1]," ")[[1]])
  sm             <- DFNET_graph$Selected_Module
  
  print(bm)
  print(sm)
 
  print(DFNET_mod)

  COVERAGE[yy,xx]    <- as.numeric(identical(bm,sm))
  #AUC[yy,xx]		 <- DFNET_acc
  N_MODULES[yy,xx]   <- dim(DFNET_mod)[1] 

  print(COVERAGE)
  #print(AUC)
  print(N_MODULES)

 }
}

# Save the entire R Workspace
#save.image("RESULTS.RData")
#save(DFNET_graph, file="Graph.RData")


# Do some Plots
par(mfrow=c(2,3))
barplot(apply(COVERAGE,2,sum)/N.Nodes, ylab="Coverage", xlab="Iterations", 
	names.arg=NITER, col="cadetblue", ylim=c(0,1), las=2)

boxplot(N_MODULES, ylab="Number of unique modules detected", xlab="Iterations", 
	names=NITER, col="cadetblue", las=2)

boxplot(AUC, ylab="AUC", xlab="Iterations", 
	names=NITER, col="cadetblue", ylim=c(0,1), las=2)

