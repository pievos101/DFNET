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

# retrieve full data
dataset <- do.call(cbind, graph[[2]])
# leave only one target column
dataset <- dataset[, -which(colnames(dataset) == "target")[-1]]
dim(dataset)

set.seed(123, kind="L'Ecuyer-CMRG")
dfnet <- DFNET(graph, ntrees=100, niter=300, init.mtry=15)

# retrieve all trees
trees <- dfnet$DFNET_trees
length(trees)

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
  as.vector(apply(pred, 2, median))
}
y_hat <- predict_dfnet(dfnet, dataset)
mean(dataset$target==y_hat) # to be expected since evaluated on training data


#-----------------------------
# install.packages("devtools")
# devtools::install_github("ModelOriented/treeshap")
library(treeshap)

# create forest from single trees
unify_forest <- function(trees, data) {
  unified_model <- data.frame()
  i <- 0
  for (tree in trees) {
    cat(i + 1, " of ", length(trees), "\n")
    unified_tree <- ranger.unify(tree, data)
    unified_tree$model$Tree <- i
    unified_model <- rbind(unified_model, unified_tree$model)
    i <- i + 1
  }
  # unified_model$Feature <- as.numeric(lapply(
  #     strsplit(forest$model$Feature, "_"), function(x) ifelse(length(x[-1]) == 0, NA, x[-1])
  # ))
  unified_forest <- unified_tree
  unified_forest$model <- unified_model
  unified_forest
}
forest <- unify_forest(dfnet$DFNET_trees, dataset)

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


# join mm features into edge importance ? 
variable_count <- dim(sv)[2]/2
sv_join <- sv[,1:variable_count] + sv[,(variable_count+1):(2*variable_count)]
colnames(sv_join) <- as.numeric(lapply(
  strsplit(colnames(sv_join), "_"), function(x) ifelse(length(x[-1]) == 0, NA, x[-1])
))
global_sv_joined <- colMeans(abs(sv_join))
df_joined <- data.frame(
  variable = factor(names(global_sv_joined)),
  importance = as.vector(global_sv_joined)
)
df_joined$variable <- reorder(df_joined$variable, df_joined$importance)
df_joined <- df_joined[order(df_joined$importance, decreasing = TRUE)[1:20], ]

p_joined <- ggplot(df_joined, aes(x = variable, y = importance)) +
  geom_bar(stat = "identity", fill = colors_discrete_drwhy(1)) +
  coord_flip() +
  theme_drwhy_vertical() +
  ylab("mean(|SHAP_A + SHAP_B|)") + xlab("") +
  labs(title = "Edge Importance") +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position = "none")
p_joined
