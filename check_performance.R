# Performance measures

DFNET_performance <- function(pred, y){

require(caret)
require(e1071)

res <- confusionMatrix(as.factor(pree), as.factor(y), mode = "prec_recall", positive="1")

return(res)

}