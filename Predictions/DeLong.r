library(pROC)

df = read.csv('/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/plt/TaskI/predictions/y_test_pred_all_models.csv', header=TRUE, sep=',',row.names=1)
y_test = factor(df[['y_test']])
copan_pred <- as.vector(df[['copan']])
print(y_test)
roc_copan <- roc(y_test, copan_pred)
models <- c('mOTUs','humann','clinical')
for (m in models) {
    print(m)
    roc_foe <- roc(y_test, as.vector(df[[m]]))
    print(roc.test(roc_foe,roc_copan,method=c('delong')))
}