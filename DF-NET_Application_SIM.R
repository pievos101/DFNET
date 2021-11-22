# RFNET Application (linekd Omics)
library(ranger)
library(igraph)
library(pROC)
library(DFNET)

## PPI
#PPI      <- read.table("~/LinkedOmics/BRCA/BRCA_PPI.txt")
PPI      <- read.table("~/LinkedOmics/KIRC/KIDNEY_PPI.txt")

## Features
#mRNA     <- read.table("~/LinkedOmics/BRCA/BRCA_mRNA_FEATURES.txt")
mRNA     <- read.table("~/LinkedOmics/KIRC/KIDNEY_mRNA_FEATURES.txt")

#Methy    <- read.table("~/LinkedOmics/BRCA/BRCA_Methy_FEATURES.txt")
Methy    <- read.table("~/LinkedOmics/KIRC/KIDNEY_Methy_FEATURES.txt")

#Mut      <- read.table("~/LinkedOmics/BRCA/BRCA_Mut_FEATURES.txt")

# Outcome class
#TARGET   <- read.table("~/LinkedOmics/BRCA/BRCA_SURVIVAL.txt")
TARGET   <- read.table("~/LinkedOmics/KIRC/KIDNEY_SURVIVAL.txt")

#@FIXME -- UGLY
# Replace NANs with mean
na.ids <- which(apply(Methy,2,function(x){any(is.na(x))}))

for (xx in na.ids){

	ids <- which(is.na(Methy[,xx]))
	Methy[ids,xx] <- mean(Methy[,xx], na.rm=TRUE)

}
#-----------------------------

# Read Data
DFNET_graph  <- DFNET_generate_graph_Omics(PPI, list(mRNA, Methy), TARGET, 0.90)

# Make data balanced -------------------------------------------- #
    TT        <- table(unlist(TARGET))
    id        <- which.min(TT)
    down_samp <- TT[id]
    class     <- as.numeric(names(TT[id]))
    down_ids  <- sample(which(unlist(TARGET)!=class), down_samp)

    ids1 <- down_ids
    ids2 <- which(unlist(TARGET)==class)

    TARGET2 <- unlist(TARGET)[c(ids1,ids2)]

    for (xx in 1:length(DFNET_graph$Feature_Matrix)){
        DFNET_graph$Feature_Matrix[[xx]] <- DFNET_graph$Feature_Matrix[[xx]][c(ids1,ids2),]
    }
    # ---------------------------------------- #


# Normalize 
DFNET_graph  <- DFNET_preprocess(DFNET_graph)


# DFNET - START SIMULATIONS ----------------------------- #
N.SIM  <- 20

DFNET_RESULT <- matrix(NaN, N.SIM, 5)
colnames(DFNET_RESULT) <- c("Sensitivity","Specificity",
                        "Precision","Recall","Accuracy")
RF_RESULT <- matrix(NaN, N.SIM, 5)
colnames(RF_RESULT) <- c("Sensitivity","Specificity",
                        "Precision","Recall","Accuracy")
for(sim in 1:N.SIM){

    
    # Create TRAIN set ----------------------------------- #
    DFNET_graph_train <- DFNET_graph
    ## 80% of the sample size
    smp_size  <- floor(0.80 * nrow(DFNET_graph$Feature_Matrix[[1]]))
    train_ids <- sample(seq_len(nrow(DFNET_graph$Feature_Matrix[[1]])), size = smp_size)
    for(xx in 1:length(DFNET_graph_train$Feature_Matrix)){
        DFNET_graph_train$Feature_Matrix[[xx]] <- DFNET_graph$Feature_Matrix[[xx]][train_ids,]
    }
    table(TARGET2[train_ids])

    # Create TEST set ------------------------------------ #
    DFNET_graph_test  <- DFNET_graph
    test_ids <- (1:nrow(DFNET_graph$Feature_Matrix[[1]]))[-train_ids]
    for(xx in 1:length(DFNET_graph_test$Feature_Matrix)){
        DFNET_graph_test$Feature_Matrix[[xx]] <- DFNET_graph$Feature_Matrix[[xx]][test_ids,]
    }
    table(TARGET2[test_ids])


    # DFNET ------------------------------------ #
    DFNET_object <- DFNET(DFNET_graph_train, ntrees=500, niter=0, init.mtry=15)

    # PREDICTION
    DFNET_pred   <- DFNET_predict(DFNET_object, DFNET_graph_test)

    # PERFORMANCE
    DFNET_perf   <- DFNET_performance(DFNET_pred, TARGET2[test_ids])

    #DFNET_perf$byClass
    DFNET_RESULT[sim,1] <- DFNET_perf$byClass["Sensitivity"]
    DFNET_RESULT[sim,2] <- DFNET_perf$byClass["Specificity"]
    DFNET_RESULT[sim,3] <- DFNET_perf$byClass["Precision"]
    DFNET_RESULT[sim,4] <- DFNET_perf$byClass["Recall"]
    DFNET_RESULT[sim,5] <- DFNET_perf$byClass["Balanced Accuracy"]

print("DFNET RESULT")
print(DFNET_RESULT)
    
    # Vanilla RF
    # TRAIN
    dataset1 <- DFNET_graph_train$Feature_Matrix[[1]]
    dataset2 <- DFNET_graph_train$Feature_Matrix[[2]]
    DATASETX <- cbind(dataset1, dataset2)
    TRAIN_DATASET  <- DATASETX[, -which(colnames(DATASETX) == "target")[-1]]
    
    vanilla_rf  <-  ranger(dependent.variable.name = "target",
                    data = TRAIN_DATASET, # MM DATA
                    classification = TRUE, 
                    importance = "impurity", 
                    num.trees = 500, 
                    mtry = 15,
                    replace = TRUE) 

    
    # Vanilla RF
    # TEST
    dataset1 <- DFNET_graph_test$Feature_Matrix[[1]]
    dataset2 <- DFNET_graph_test$Feature_Matrix[[2]]
    DATASETX <- cbind(dataset1, dataset2)
    TEST_DATASET  <- DATASETX[, -which(colnames(DATASETX) == "target")[-1]]
    
    pp    <- predict(vanilla_rf, TEST_DATASET)
    pred  <- pp$predictions

    RF_perf <- DFNET_performance(as.factor(pred), as.factor(DFNET_graph_test$Feature_Matrix[[1]][,"target"]))
    RF_RESULT[sim,1] <- RF_perf$byClass["Sensitivity"]
    RF_RESULT[sim,2] <- RF_perf$byClass["Specificity"]
    RF_RESULT[sim,3] <- RF_perf$byClass["Precision"]
    RF_RESULT[sim,4] <- RF_perf$byClass["Recall"]
    RF_RESULT[sim,5] <- RF_perf$byClass["Balanced Accuracy"]

print("VANILLA RF RESULT")
print(RF_RESULT)

}

dev.off()
par(mfrow=c(1,2))
boxplot(RF_RESULT, ylim=c(0,1), las=2, main="VANILLA RF")
boxplot(DFNET_RESULT, ylim=c(0,1), las=2, main="DFNET")

