# treeSHAP plots
# forest_shap <- DFNET_explain()

sv <- forest_shap$shaps

# plot
library(ggplot2)

global_sv <- colMeans(abs(sv))

# convert to gene names
NN <- names(global_sv)
gene_names <- DFNET_graph_test$gene.names
NN2 <- strsplit(NN, "_")
NN3 <- sapply(NN2, function(x) {
    return(x[1])
})
NN4 <- gsub("AN", "mRNA", NN3)
NN4 <- gsub("BN", "Methy", NN4)
NN5 <- paste(gene_names, NN4, sep = "_")
names(global_sv) <- NN5

df <- data.frame(variable = factor(names(global_sv)), importance = as.vector(global_sv))
df$variable <- reorder(df$variable, df$importance)
df <- df[order(df$importance, decreasing = TRUE)[1:15], ]
p <- ggplot(df, aes(x = variable, y = importance)) +
    geom_bar(stat = "identity", fill = colors_discrete_drwhy(1)) +
    coord_flip() +
    theme_drwhy_vertical() +
    ylab("mean(|SHAP|)") +
    xlab("") +
    labs(title = "Feature Importance") +
    scale_y_continuous(labels = scales::comma) +
    theme(legend.position = "none")
p



# join mm features into node importance
variable_count <- dim(sv)[2] / 2
sv_join <- sv[, 1:variable_count] + sv[, (variable_count + 1):(2 * variable_count)]
colnames(sv_join) <- as.numeric(lapply(
    strsplit(colnames(sv_join), "_"), function(x) ifelse(length(x[-1]) == 0, NA, x[-1])
))
global_sv_joined <- colMeans(abs(sv_join))
names(global_sv_joined) <- DFNET_graph_test$gene.names[-length(DFNET_graph_test$gene.names)]
df_joined <- data.frame(
    variable = factor(names(global_sv_joined)),
    importance = as.vector(global_sv_joined)
)
df_joined$variable <- reorder(df_joined$variable, df_joined$importance)
df_joined <- df_joined[order(df_joined$importance, decreasing = TRUE)[1:15], ]

p_joined <- ggplot(df_joined, aes(x = variable, y = importance)) +
    geom_bar(stat = "identity", fill = colors_discrete_drwhy(1)) +
    coord_flip() +
    theme_drwhy_vertical() +
    ylab("mean(|SHAP_A + SHAP_B|)") +
    xlab("") +
    labs(title = "Node Importance") +
    scale_y_continuous(labels = scales::comma) +
    theme(legend.position = "none")
p_joined


# local explanation
treeshap::plot_contribution(forest_shap, 20)
