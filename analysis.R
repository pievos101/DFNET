library(ranger)
library(igraph)
library(pROC)


PPI      <- read.table("./data/KIDNEY_PPI.txt")
mRNA     <- read.table("./data/KIDNEY_mRNA_FEATURES.txt")
Methy    <- read.table("./data/KIDNEY_Methy_FEATURES.txt")
TARGET   <- read.table("./data/KIDNEY_SURVIVAL.txt")

#@FIXME -- UGLY
# Replace NANs with mean
na.ids <- which(apply(Methy,2,function(x){any(is.na(x))}))

for (xx in na.ids){
  
  ids <- which(is.na(Methy[,xx]))
  Methy[ids,xx] <- mean(Methy[,xx], na.rm=TRUE)
  
}


#-----------------------------
source("source.R")


graph  <- DFNET_generate_graph_Omics(PPI, list(mRNA, Methy), TARGET, 0.95)
saveRDS(graph, "graph.rds")

graph <- readRDS("graph.rds")

# retrieve full data
dataset <- do.call(cbind, graph[[2]])
# leave only one target column
dataset <- dataset[, -which(colnames(dataset) == "target")[-1]]
dim(dataset)

set.seed(123, kind="L'Ecuyer-CMRG")
ntrees <- 100

dfnet <- DFNET(graph, ntrees=ntrees, niter=300, init.mtry=15)
saveRDS(dfnet, "dfnet.rds")

dfnet <- readRDS("dfnet.rds")

eimp <- DFNET_Edge_Importance(graph, dfnet)
saveRDS(eimp, "eimp.rds")

eimp <- readRDS("eimp.rds")

modules <- DFNET_modules(graph, dfnet, eimp)
saveRDS(modules, "modules.rds")

modules <- readRDS("modules.rds")

nodes <- unlist(stringr::str_split(modules$Module, " "))
length(unique(nodes))
sum(table(nodes) > 1)
# frequency of nodes
hist(table(nodes))

# retrieve all trees
trees <- dfnet$DFNET_trees
length(trees)

sapply(trees, function(x) length(x$variable.importance))

single_tree <- trees[[101]]
single_tree
y_hat <- predict(single_tree, dataset)
y_hat$predictions

oob_error <- mean(sapply(trees, function(x) x$prediction.error))
oob_error # 0.66 out-of-bag accuracy

# predict outcome based on forest
predict_dfnet <- function(model, newdata) {
  pred <- data.frame()
  for (tree in model$DFNET_trees) {
    pred <- rbind(pred, predict(tree, newdata)$predictions)
  }
  as.vector(apply(pred, 2, mean))
}
y_hat <- predict_dfnet(dfnet, dataset)
mean(dataset$target==(y_hat>0.5)) # to be expected since evaluated on training data


#-----------------------------
# install.packages("devtools")
# devtools::install_github("ModelOriented/treeshap")
library(treeshap)

test <- ranger.unify(single_tree, dataset[, names(single_tree$variable.importance)])
test$model

# create a forest from single trees
# it will be a unified (data.frame) representation with treeshap
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
forest <- unify_forest(trees[(length(trees)-ntrees+1):length(trees)], dataset)
saveRDS(forest, "forest.rds")

forest <- readRDS("forest.rds")

# calculate treeshap
forest_shap <- treeshap(forest, dataset[, -dim(dataset)[2]])
sv <- forest_shap$shaps


# plot
library(ggplot2)

global_sv <- colMeans(abs(sv))
df <- data.frame(variable = factor(names(global_sv)), importance = as.vector(global_sv))
df$variable <- reorder(df$variable, df$importance)
df <- df[order(df$importance, decreasing = TRUE)[1:20], ]
p <- ggplot(df, aes(x = variable, y = importance)) +
  geom_bar(stat = "identity", fill = colors_discrete_drwhy(1)) +
  coord_flip() +
  theme_drwhy_vertical() +
  ylab("mean(|SHAP|)") + xlab("") +
  labs(title = "Feature Importance") +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position = "none")
p



# local explanation
treeshap::plot_contribution(forest_shap, 1)
tmp <- dfnet
tmp$DFNET_trees <- trees[(length(trees)-ntrees+1):length(trees)]
predict_dfnet(tmp, dataset[300, ])
treeshap:::predict.model_unified(forest, dataset[300, ])
treeshap::plot_contribution(forest_shap, 300)







### module clustering experiment

f <- forest$model
variable_count <- (dim(dataset)[2] - 1) / 2
first <- TRUE
for (tree in split(f, f$Tree)) {
  nodes <- as.numeric(stringr::str_split(
    as.character(na.omit(unique(tree$Feature))),
    "_",
    simplify = TRUE
  )[, 2])
  row <- ifelse(1:variable_count %in% nodes, 1, 0)
  if (!first) {
    module_node_matrix <- rbind(module_node_matrix, row)
  } else {
    module_node_matrix <- matrix(row, nrow = 1)
    first <- FALSE
  }
}
rownames(module_node_matrix) <- 1:dim(module_node_matrix)[1]
colnames(module_node_matrix) <- 1:dim(module_node_matrix)[2]
dim(module_node_matrix)
save_module_node_matrix <- module_node_matrix

table(apply(save_module_node_matrix, 2, function(x) sum(x)))
module_node_matrix <- save_module_node_matrix[, apply(save_module_node_matrix,
                                                      2,
                                                      function(x) sum(x) > 0)]

module_node_dist <- dist(module_node_matrix, method = "manhattan")
clust_modules <- hclust(module_node_dist)
K <- 3
clust_modules_labels <- cutree(tree = clust_modules, k = K)
plot(x = clust_modules, main = "Clustering modules by them having similar variables",
     labels =  row.names(clust_modules), cex = 0.5)
rect.hclust(tree = clust_modules, k = K, border = 1:K, cluster = clust_modules_labels)


library(genieclust) # genieclust algorithm penalizes small clusters
genie_modules <- gclust(module_node_dist, 0.3, distance = "manhattan")
genie_modules_labels <- cutree(tree = genie_modules, k = K)
table(genie_modules_labels)
plot(genie_modules)


### variable clustering experiment

node_module_matrix <- t(module_node_matrix)
dim(node_module_matrix)
table(apply(node_module_matrix, 1, function(x) sum(x)))
node_module_matrix <- node_module_matrix[apply(node_module_matrix, 1, function(x) sum(x) > 1),]
dim(node_module_matrix)
node_module_dist <- dist(node_module_matrix, method = "manhattan")

clust_nodes <- hclust(node_module_dist)
clust_nodes_labels <- cutree(tree = clust_nodes, k = K)
table(clust_nodes_labels)
plot(x = clust_nodes, main = "Clustering variables by them being in similar modules",
     labels =  row.names(clust_nodes), cex = 0.5)


genie_nodes <- gclust(node_module_dist, 0.3, distance = "manhattan")
genie_nodes_labels <- cutree(tree = genie_nodes, k = K)
table(genie_nodes_labels)
plot(genie_nodes)